set.seed(12345)
library(caret)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(microbiome)
library(doParallel)
library(Boruta)
library(DMwR)
#library(doRNG)

args<-commandArgs(TRUE)

infile <- as.character(args[1]) #     infile <- "./Test_data/infiles/phylo_training_set.rds"
Variable <- as.character(args[2]) #     Variable <- "condition"
PrevCutoff <- as.numeric(as.character(args[3])) # PrevCutoff <-  0.1
Feature_selection <- as.logical(as.character(args[4])) # Feature_selection <- T
ncores  <- as.numeric(as.character(args[5])) # ncores <-  15
simple_training_params = F
print(infile)
print(Variable)
print(PrevCutoff)
print(Feature_selection)
print(ncores)

#######################################################################################################################################
##########################################              Load the functions                       ######################################
#######################################################################################################################################

#path_functions <- "/home/luna.kuleuven.be/u0141268/github_shared_code_and_publications/Liver_Disease_Microbiome_ML/Functions"
path_functions  <- as.character(args[6]) 
source(paste0(path_functions,"/Machine_learning_functions.R"))

#######################################################################################################################################
##########################################             ML running parameters                     ######################################
#######################################################################################################################################

#### Parameters classifier ####
n_K <- 10 ### Number of CV in the outer loop
#n_K <- 2 ### Number of CV in the outer loop
cv_folds = 5  # cv_folds = 2   
cv_repeats = 5 #  cv_repeats = 1  

if(simple_training_params == T){
	n_K <- 10 ### Number of CV in the outer loop
	cv_folds = 2  # cv_folds = 2   
	cv_repeats = 1 #  cv_repeats = 1  
}

#######################################################################################################################################
##########################################                  Read the data                        ######################################
#######################################################################################################################################

if(grepl(".rds$",infile)){
	in.phylo <- readRDS(infile)
}

if( grepl(".RData$",infile) | grepl(".Rdata$",infile) ){
	in.phylo <- load_RData(file=infile)
}

#### Re-format the data into a data.frame format ####
MEspDF <-  data_format(in.obj = in.phylo, PrevCutoff = PrevCutoff, Variable=Variable )

if(is.character(MEspDF$Variable)){
	# ### Set the reference and the postive class, the positve class is always the second level of the factor
	MEspDF$Variable <- paste0(Variable,".",MEspDF$Variable)	
	MEspDF$Variable <- factor(MEspDF$Variable)
}

cat("\n","N samples = ",nrow(MEspDF), "\n" )

#######################################################################################################################################
##########################################            Check if the data is balanced              ######################################
#######################################################################################################################################
## https://thesai.org/Downloads/Volume13No6/Paper_27-Survey_on_Highly_Imbalanced_Multi_class_Data.pdf
## According to Hamid et al. (2022), datasets become “moderately to highly imbalanced” when the majority-to-minority ratio exceeds 80:20 and reach “high imbalance” at 50:1 or more. Extreme imbalance is observed at ratios of 1,000:1 or higher. My rule—flagging imbalance when a minority class is less than one-fifth the size of the majority—automatically captures cases starting at an 80:20 imbalance, which aligns well with standard definitions for moderate imbalance.
##  Then implement the 1:5 imbalance rule
#Frequency table
if(class(MEspDF$Variable) == "factor" | class(MEspDF$Variable) == "character" ){

	#Frequency table
	tbl <- table(MEspDF$Variable)

	# total samples
	n <- sum(tbl)

	# check if the smallest class proportion < 1/5 of largest
	is_imbalanced <- (min(tbl) / max(tbl)) < (1/5)

	if (is_imbalanced == T) {
		cat("\nData is imbalanced\n")
		cat("\nmin/max class ratio:",  (min(tbl) / max(tbl)),"\n" )
	} else {
		cat("\nData is balanced\n")
		cat("\nmin/max class ratio:",  (min(tbl) / max(tbl)),"\n" )
	}

	check_and_report_balance(tbl=tbl, output_file = "data_balance_report.txt",verbose = TRUE  )
	
}

#######################################################################################################################################
##########################################            Create the outfiles directory              ######################################
#######################################################################################################################################

out_dir <- paste0("RandomForest_FS_",Feature_selection)
dir.create( out_dir )
setwd(out_dir)

#######################################################################################################################################
##########################################             Create the outer loop                     ######################################
#######################################################################################################################################

outer_loop <- createFolds(MEspDF$Variable, k = n_K, list = TRUE, returnTrain = T)
#save(file="outer_loop_cv_sample_index.RData",outer_loop)
saveRDS(file="outer_loop_cv_sample_index.rds",outer_loop)

#######################################################################################################################################
##########################################             Set the multithread enviroment                  ################################
#######################################################################################################################################

cl <- parallel::makeForkCluster(ncores)
doParallel::registerDoParallel(cl)
#registerDoRNG(12345) 

#######################################################################################################################################
##########################################             Create the inner loop                     ######################################
#######################################################################################################################################

fitControl <- trainControl(## 5-fold CV for the inner loop
                           method = "repeatedcv",
                           number = cv_folds,
                           repeats = cv_repeats,
                           allowParallel = T)

#######################################################################################################################################
##########################################              Run the inner loop                       ######################################
#######################################################################################################################################

#### List and data.frames to save the results ####
list_CV_models <- list()
stats_model <- data.frame()
sum <- 1

#### Stats and plot varaibles ####
if( class(MEspDF$Variable) == "factor" | class(MEspDF$Variable) == "character"  ){
	list_prediction <- list()
	list_response <- list()
	list_ggroc_auc <- c()
	list_ggroc_curves <- list()
	pooled_outer_loop_cv_predictions <- data.frame()
}
importance_df <- data.frame()

#for(outer in names(outer_loop)){
for(i in seq_along(names(outer_loop))){

	# outer <- "Fold01"  Fold1  i<-7
	outer <- names(outer_loop)[i]	
	cat("Training",outer," Feature selection ",Feature_selection,"\n")
	
	ind_train <- outer_loop[[outer]]
	Train <- MEspDF[ind_train,]
	Test <- MEspDF[-ind_train,]


	if(Feature_selection == T){
		set.seed(12345 + i * 100)  # Unique seed
		# Boruta feature selection
		print("START Boruta")	
		boruta_model <- Boruta(Variable ~ ., data = Train, doTrace = 0)
		final_vars <- getSelectedAttributes(boruta_model, withTentative = TRUE)
		cat("N FEATURES: ", length(final_vars), "\n")
		print("END Boruta")

		### Handle zero features selected ###
		if (length(final_vars) == 0) {
			warning("Fold ", i, ": Boruta selected 0 features. Using all features instead.")
			final_vars <- colnames(Train)
			final_vars <- final_vars[!grepl("Variable", final_vars)]
		}

	}else{
		final_vars <- colnames(Train)
		final_vars <- final_vars[!grepl("Variable",final_vars)]	
	}

	final_vars_sub <- c("Variable", final_vars)
	Train_selected <- Train[, match(final_vars_sub, colnames(Train))]
	Test_selected <- Test[, match(final_vars_sub, colnames(Test))]

	###########################
	#####  Apply SMOTE   ######
	###########################
	if (is_imbalanced == T) {
		cat("Applying SMOTE\n")
		Train_selected$Variable <- as.factor(Train_selected$Variable)
		set.seed(12345 + i * 1000 + i * 10 + 1)   # Unique seed
		Train_selected <- smart_smote(data = Train_selected)          
	} 
	
	###########################
	##### Train the model #####
	set.seed(12345 + i * 1000 + i * 10 + 2)  
	Fit_model <- train(Variable ~ ., data = Train_selected, method = "rf", trControl =    fitControl )
	# Save the model
	list_CV_models[[sum]] <- Fit_model
	
	###########################
	##### Test the model #####
	tempStats <- stats_model_func(Fit.Model = Fit_model, Test.df = Test_selected, RetDF = T, df_complete = MEspDF)
	tempStats <- cbind( Fold = outer , tempStats)
	print(tempStats)
	
	############################################################
	##### Save the probs for the ROC and performance plots #####
	if( class(MEspDF$Variable) == "factor" | class(MEspDF$Variable) == "character"  ){

		PredRF_VAR <- predict(Fit_model, newdata =  Test_selected )
		PredRFProb <- predict(Fit_model, newdata = Test_selected,  type = "prob")
		list_prediction[[sum]] <- PredRFProb[, levels(Train_selected$Variable)[2] ]
		list_response[[sum]] <- PredRF_VAR
		
		### stats for all the predictions in the outer loop, will help to check inconsistencies
		pooled_outer_loop_cv_predictions <- rbind( pooled_outer_loop_cv_predictions,
						data.frame( Real = Test_selected$Variable, Predictec = PredRF_VAR , PredRFProb ) )
		
		if(length(unique(MEspDF$Variable)) == 2){		
			rocobj2 <- pROC::roc(Test$Variable, PredRFProb[, levels(Train_selected$Variable)[2]] )
			list_ggroc_auc[sum] <- as.numeric(auc(rocobj2))
			list_ggroc_curves[[sum]]  <- rocobj2	
		}else{
			ret_list <- ROC_AUC_func(Test_df = Test, pred_probs = PredRFProb )
			list_ggroc_auc[[sum]] <- ret_list$auc_list
			list_ggroc_curves[[sum]] <- ret_list$roc_list
		}
	}

	##############################################
	##### Save the importance for each model #####
	
	# Check number of features before calculating variable importance
	# varImp() fails with a single feature (returns NaN)
	n_features <- length(Fit_model$coefnames)

	if(n_features > 1){
	    # Multiple features - calculate importance normally
	    importance_df <- rbind(importance_df, data.frame(outer, Feature = rownames(varImp(Fit_model)$importance), varImp(Fit_model)$importance))
	} else {
	    # Single feature - assign 100% importance by default
	    importance_df <- rbind(importance_df, data.frame(outer, Feature = Fit_model$coefnames, Overall = 100))
	}

	##########################
	##### Save the stats #####	
	stats_model <- rbind(stats_model,tempStats)
	sum <- 1 + sum
	
	rm(Fit_model, ind_train, Train, Test, final_vars_sub, final_vars,Train_selected, Test_selected,tempStats)
}

names(list_CV_models) <- names(outer_loop)
#save(file=paste0("list_CV_models",".FeatureSelection_",Feature_selection,".RData"),list_CV_models)
#save(file= "list_CV_models.RData",list_CV_models)
saveRDS(file= "list_CV_models.rds",list_CV_models)

saveRDS(file= "stats_model.rds",stats_model)
#save(file= "stats_model.RData",stats_model)
write.table(stats_model,"stats_model.tsv",col.names=T,row.names = F,quote=FALSE,sep = "\t")

saveRDS(file= "pooled_outer_loop_cv_predictions.rds",pooled_outer_loop_cv_predictions)

### Save the stats and plots for the pooled apprach
if( class(MEspDF$Variable) == "factor" | class(MEspDF$Variable) == "character"  ){

	# Detect number of classes and select probability columns dynamically
	n_classes <- length(unique(pooled_outer_loop_cv_predictions$Real))
	prob_cols <- 3:(2 + n_classes)
	
	# Fix column names to match class levels
	colnames(pooled_outer_loop_cv_predictions)[prob_cols] <- levels(factor(pooled_outer_loop_cv_predictions$Real))
	
	# Calculate stats
	stats_pooled_model <- stats_model_func_discrete(
		PredProbs = pooled_outer_loop_cv_predictions[, prob_cols], 
		PredClass = pooled_outer_loop_cv_predictions$Predictec, 
		RealClass = pooled_outer_loop_cv_predictions$Real, 
		return_data.frame = TRUE
	)
	write.table(stats_pooled_model, "stats_pooled_model.tsv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
	saveRDS(stats_pooled_model, file = "stats_pooled_model.rds")
	
	p_roc <- plot_roc_simple_curve(
		PredProbs = pooled_outer_loop_cv_predictions[, prob_cols], 
		PredClass = pooled_outer_loop_cv_predictions$Predictec, 
		RealClass = pooled_outer_loop_cv_predictions$Real, 
		mcc_value = stats_pooled_model$MCC
	)

	p_roc <- p_roc + labs(title = "Pooled outer loop cv ROC Curve")
	ggsave("pooled_outer_loop_roc_curve.pdf", p_roc, width = 12, height = 7)
	saveRDS(file="pooled_outer_loop_roc_curve.rds",p_roc)
}

###################################################################################################################################
##############################              Average performance metrics in the outer-loop               ###########################
###################################################################################################################################

Summary_Statistics_Nested_CV <- Stats_CV_function(Stats_CV = stats_model )

# Save to file
write.table(Summary_Statistics_Nested_CV, "Summary_Statistics_Nested_CV.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Plot the stats
stats_names <- colnames(stats_model)[-1]
plot_list <- vector("list", length(colnames(stats_model)[-1]))
for(i in 1:length(stats_names)){ plot_list[[i]] <- create_metric_boxplot(data=stats_model, metric_name=stats_names[i]) }
# Merge
p_combined <- ggarrange(plotlist = plot_list)
# Save
ggsave("Metrics_CV_Boxplots.pdf", p_combined, width = 10, height = 8)

#################################################################################
################    Plot ROC curves in case of classification    ################
if( class(MEspDF$Variable) == "factor" | class(MEspDF$Variable) == "character"  ){
	
	plot_auc_function(list_ggroc_curves=list_ggroc_curves, list_ggroc_auc = list_ggroc_auc , 
		Variable = Variable,  Categories_variable = unique(as.character(MEspDF$Variable)))
}

#########################################################
#############        Importance plot        #############
#########################################################

if( class(MEspDF$Variable) == "factor" | class(MEspDF$Variable) == "character"  ){
	Imp_plot_discrete_function(importance_df = importance_df,  MEspDF = MEspDF)
}

##################################################################################################
#############        Performance of the models along all the outer loops plot        #############
##################################################################################################


##### IMPORTANT ######
# Exact run-to-run reproducibility is not the primary objective.
# Due to the stochastic nature of the pipeline (Random Forest, SMOTE,
# and feature selection), small performance differences across runs
# are expected and reflect model instability.
# Model selection is therefore based on mean performance and variability
# across outer cross-validation folds.


# Retrain nested CV models on outer loops and evaluate
#ctrl <- trainControl(method = "none")  # No CV, just train
#ctrl <- fitControl
Stats_CV_models_outer_loops <- data.frame()

#for(outer in names(outer_loop)){
for(i in seq_along(names(outer_loop))){
	# outer <- "Fold01"  # i<-1
	outer <- names(outer_loop)[i]
	print(outer)
	
	#### Data partitions #####
	ind_train <- outer_loop[[outer]]
	Train <- MEspDF[ind_train,]
	Test <- MEspDF[-ind_train,]
	
	for(j in seq_along(list_CV_models)){
		# Extract previously trained model
		modelList <- names(list_CV_models)[j]
		Fit_model <- list_CV_models[[modelList]]
		
		# Select same features used by this model
		Train_selected <- Train[, match(c("Variable", Fit_model$coefnames), colnames(Train))]
		Test_selected <- Test[, match(c("Variable", Fit_model$coefnames), colnames(Test))]
		
		###########################
		#####  Apply SMOTE   ######
		###########################
		if (is_imbalanced == T) {
			cat("Applying SMOTE\n")
			Train_selected$Variable <- as.factor(Train_selected$Variable)
			set.seed(12345 + i * 1000 + j * 10 + 1)
			Train_selected <- smart_smote(data = Train_selected)          
		} 


		set.seed(12345 + i * 1000 + j * 10 + 2)  
		Fit_model_RT <- train(Variable ~ ., data = Train_selected, method = "rf", trControl = fitControl, tuneGrid = Fit_model$bestTune )
		
		# Calculate performance stats on test set
		tempStats <- stats_model_func(
			Fit.Model = Fit_model_RT, 
			Test.df = Test_selected, 
			RetDF = T, 
			df_complete = MEspDF
		)
		
		# Add hyperparameters to stats
		tempStats <- cbind(tempStats, Fit_model_RT$bestTune)
		tempStats <- data.frame(Model_Fold = modelList, Test_Fold = outer, tempStats)
		
		# Accumulate results
		Stats_CV_models_outer_loops <- rbind(Stats_CV_models_outer_loops, tempStats)
		
		# Clean up
		rm(tempStats, Fit_model, Fit_model_RT, Train_selected, Test_selected)
	}
	
	rm(ind_train, Train, Test)
}

#save(file = "Stats_CV_models_outer_loops.RData", Stats_CV_models_outer_loops)
saveRDS(file = "Stats_CV_models_outer_loops.rds", Stats_CV_models_outer_loops)

################################################################
### Check the diferences bewteen the first and second model ####
stats_model2 <- Stats_CV_models_outer_loops[Stats_CV_models_outer_loops$Model_Fold == Stats_CV_models_outer_loops$Test_Fold,]
sink("Different_models.txt")
cat("Mean MCC diff = ",  mean(stats_model$MCC - stats_model2$MCC, na.rm = TRUE), "; sd = ",  sd(stats_model$MCC - stats_model2$MCC, na.rm = TRUE), "\n" )
cat("Mean AUC diff = ",  mean(stats_model$AUC - stats_model2$AUC, na.rm = TRUE), "; sd = ",  sd(stats_model$AUC - stats_model2$AUC, na.rm = TRUE), "\n" )
sink()


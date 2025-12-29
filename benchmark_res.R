set.seed(12345)
library(caret)
library(phyloseq)
library(microbiome)

path_functions <- "/home/luna.kuleuven.be/u0141268/github_shared_code_and_publications/Liver_Disease_Microbiome_ML/Functions"
source(paste0(path_functions,"/Machine_learning_functions.R"))

Variable <- "normal_vs_cirrhosis"
PrevCutoff <- 0
#######################################################################################################################################
##########################################                  Read the data                        ######################################
#######################################################################################################################################

###########################################
###### phylo_holdout_validation_set  ######
infile <- "/home/luna.kuleuven.be/u0141268/github_shared_code_and_publications/Liver_Disease_Microbiome_ML/Results/infiles/phylo_holdout_validation_set.rds"

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

cat("\n","Hould out N samples = ",nrow(MEspDF), "\n" )

table(MEspDF$Variable)

#################################
###### phylo_training_set  ######
infile <- "/home/luna.kuleuven.be/u0141268/github_shared_code_and_publications/Liver_Disease_Microbiome_ML/Results/infiles/phylo_training_set.rds"

if(grepl(".rds$",infile)){
	in.phylo <- readRDS(infile)
}

if( grepl(".RData$",infile) | grepl(".Rdata$",infile) ){
	in.phylo <- load_RData(file=infile)
}

#### Re-format the data into a data.frame format ####
MEspDF_training <-  data_format(in.obj = in.phylo, PrevCutoff = PrevCutoff, Variable=Variable )

if(is.character(MEspDF_training$Variable)){
	# ### Set the reference and the postive class, the positve class is always the second level of the factor
	MEspDF_training$Variable <- paste0(Variable,".",MEspDF_training$Variable)	
	MEspDF_training$Variable <- factor(MEspDF_training$Variable)
}

cat("\n","Training set N samples = ",nrow(MEspDF_training), "\n" )

table(MEspDF_training$Variable)

##########################################################
################      Average stats       ################

#### Pooled AUC
pooled_outer_loop_cv_predictions <- readRDS("/home/luna.kuleuven.be/u0141268/github_shared_code_and_publications/Liver_Disease_Microbiome_ML/Results/RandomForest_FS_TRUE/pooled_outer_loop_cv_predictions.rds")
roc_all <- pROC::roc(pooled_outer_loop_cv_predictions$Real, pooled_outer_loop_cv_predictions$normal_vs_cirrhosis.cirrhosis)

cat("\n","Pooled AUC = ",roc_all$auc, "\n" )


### Mean AUC
stats_model <- readRDS("/home/luna.kuleuven.be/u0141268/github_shared_code_and_publications/Liver_Disease_Microbiome_ML/Results/RandomForest_FS_TRUE/stats_model.rds")

cat("\n","Mean AUC = ",mean(stats_model$AUC), "\n" )



#################################################################
################      Try different models       ################


BEST_MODEL_FULL <- readRDS("/home/luna.kuleuven.be/u0141268/github_shared_code_and_publications/Liver_Disease_Microbiome_ML/Results/RandomForest_FS_TRUE/BEST_MODEL_FULL.rds")

bmstats_df <- stats_model_func(Fit.Model = BEST_MODEL_FULL, Test.df = MEspDF, RetDF = T, df_complete = MEspDF)
bmstats_df$AUC


load("/home/luna.kuleuven.be/u0141268/github_shared_code_and_publications/Liver_Disease_Microbiome_ML/Results/RandomForest_FS_TRUE/list_CV_models.RData")

print("Prediction CV over the outer loop")
for(i in list_CV_models){
	tempStats <- stats_model_func(Fit.Model = i, Test.df = MEspDF, RetDF = T, df_complete = MEspDF)$AUC
	print(tempStats)
}


print("Prediction CV over the outer loop retraining the model")
ctrl <- trainControl(method = "none")

for(model in list_CV_models){
	Best_features <- colnames(model$trainingData)[-1]
	DF_TRAIN <- MEspDF_training[,match( c( "Variable" , Best_features ) , colnames(MEspDF_training) )]
	
	BEST_MODEL_FULL <- train(Variable ~ ., data = DF_TRAIN, method = "rf", trControl = ctrl, tuneGrid = model$bestTune)	
	
	tempStats <- stats_model_func(Fit.Model = BEST_MODEL_FULL, Test.df = MEspDF, RetDF = T, df_complete = MEspDF)
	print(tempStats$AUC)	

}



set.seed(12345)
library(caret)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(microbiome)
library(doParallel)
library(Boruta)
library(DMwR)

args<-commandArgs(TRUE)

infile <- as.character(args[1]) #     infile <- "/home/luna.kuleuven.be/u0141268/Dropbox/Jiyeon_Jorge/research/research_liver_diseases/Meta_analysis/ML_predictions/infiles_case_controls/cirrhosis_species_abd_adj_prop_phyloseq.rds"  infile <- "/home/luna.kuleuven.be/u0141268/github_shared_code_and_publications/Liver_Disease_Microbiome_ML/Results/infiles/phylo_training_set.rds"
Variable <- as.character(args[2]) #     Variable <- "normal_vs_cirrhosis"
PrevCutoff <- as.numeric(args[3]) # PrevCutoff <-  0.2
ncores  <- as.numeric(args[4]) # ncores <-  15

print(infile)
print(Variable)
print(PrevCutoff)
print(ncores)

# path_functions <- "/home/luna.kuleuven.be/u0141268/github_shared_code_and_publications/Liver_Disease_Microbiome_ML/Functions"
path_functions  <- as.character(args[5]) 
source(paste0(path_functions,"/Machine_learning_functions.R"))

### Indir: Directory with the results ###
WorkingDir  <- as.character(args[6])  # WorkingDir <- "./RandomForest_FS_TRUE"

### Machine laerning method ###
MethodML <- as.character(args[7]) # MethodML <- "rf"

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

##########################################################################################################################################
##########################################                  Change to the WD                        ######################################
##########################################################################################################################################

setwd(WorkingDir)

#######################################################################################################################################
##########################################            Check if the data is balanced              ######################################
#######################################################################################################################################
## https://thesai.org/Downloads/Volume13No6/Paper_27-Survey_on_Highly_Imbalanced_Multi_class_Data.pdf
## According to Hamid et al. (2022), datasets become “moderately to highly imbalanced” when the majority-to-minority ratio exceeds 80:20 and reach “high imbalance” at 50:1 or more. Extreme imbalance is observed at ratios of 1,000:1 or higher. My rule—flagging imbalance when a minority class is less than one-fifth the size of the majority—automatically captures cases starting at an 80:20 imbalance, which aligns well with standard definitions for moderate imbalance.
##  Then implement the 1:5 imbalance rule

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
	
	Task <- "classification"
}

##################################################################################################################################################
################################                  Read the train data partitions and models                       ################################
##################################################################################################################################################

# Read the NCV models
list_CV_models <- readRDS("list_CV_models.rds")

# Read the data partitions
outer_loop_cv_sample_index <- readRDS("outer_loop_cv_sample_index.rds")

# Read the stats for the outer loop
#load("Stats_CV_models_outer_loops.RData")
Stats_CV_models_outer_loops <- readRDS("Stats_CV_models_outer_loops.rds")

#######################################################################################################################################
##########################################             Set the multithread enviroment                  ################################
#######################################################################################################################################

cl <- parallel::makeForkCluster(ncores)
doParallel::registerDoParallel(cl)

#############################################################################################################################################
################################                  Best parameters for the final models                       ################################
#############################################################################################################################################

####################################################################################
#### Criterion 1: MODE hyperparameter values across models and Feature selected ####
stats_model <- readRDS("stats_model.rds")
df_best_hyperparameter <- best_hyperparameter_selection(List_models = list_CV_models , Performance = stats_model, tie_threshold = 0.01, verbose = TRUE  )
Best_features <- feature_selection( List_models = list_CV_models, cutoff = 0.6 )

# Train the model
BEST_MODEL_NestCV_common_features <- train_final_model(Best_features = names(Best_features), df_best_hyperparameter = df_best_hyperparameter, 
	in_df = MEspDF, ML_method = MethodML, is_imbalanced = is_imbalanced, cv_folds = 5, cv_repeats = 5 )

# Save the model
saveRDS(file = "BEST_MODEL_NestCV_common_features.rds",BEST_MODEL_NestCV_common_features)

##########################################################################
#### Criterion 2: Hyperparameters of the single best-performing model ####

BestModel <- select_best_model( stats_model = stats_model)
Best_NCV_Model <- list_CV_models[[BestModel$Fold]]

Best_features_BNCV <- colnames(Best_NCV_Model$trainingData) 
Best_features_BNCV <- Best_features_BNCV[!Best_features_BNCV %in% ".outcome"]

# Train the model
BEST_MODEL_NestCV <- train_final_model(Best_features = Best_features_BNCV, 
			df_best_hyperparameter = Best_NCV_Model$bestTune, 
			in_df = MEspDF, ML_method = MethodML, is_imbalanced = is_imbalanced )

# Save the model
saveRDS(file = "BEST_MODEL_NestCV.rds",BEST_MODEL_NestCV)

#######################################################################################################################################
#### Criterion 3: Hyperparameters of the model that, on average, achieved the best performance across outer cross-validation folds ####

# Read the stats for the outer loop
#load("Stats_CV_models_outer_loops.RData")
#Stats_CV_model_FIT <- Stats_CV_models_outer_loops
Stats_CV_model_FIT <- readRDS("Stats_CV_models_outer_loops.rds")


Best_outerloop_Model <- best_performances_outer_loop_function( Stats_CV_model_FIT = Stats_CV_model_FIT, Task = Task, list_CV_models = list_CV_models )
Best_features_outerloop_performance <- colnames(Best_outerloop_Model$trainingData) 
Best_features_outerloop_performance <- Best_features_outerloop_performance[!Best_features_outerloop_performance %in% ".outcome"]

# Train the model
BEST_MODEL_outerloop_performance <- train_final_model(Best_features = Best_features_outerloop_performance, 
	df_best_hyperparameter = Best_outerloop_Model$bestTune, 
	in_df = MEspDF, ML_method = MethodML, 
	is_imbalanced = is_imbalanced )

# Save the model
saveRDS(file = "BEST_MODEL_outerloop_performance.rds",BEST_MODEL_outerloop_performance)








set.seed(12345)
library(caret)
library(phyloseq)
library(microbiome)
library(pROC)
library(ggplot2)

args<-commandArgs(TRUE)

infile <- as.character(args[1]) #     infile <- "./infiles/phylo_holdout_validation_set.rds"
Variable <- as.character(args[2]) #   Variable <- "condition"
Model_file <- as.character(args[3]) #  Model_file <- "./RandomForest_FS_TRUE/BEST_MODEL_outerloop_performance.rds"
PrevCutoff <- 0
OutDir <- as.character(args[4])  # OutDir <- "BEST_MODEL_outerloop_performance_stats"
#######################################################################################################################################
##########################################              Load the functions                       ######################################
#######################################################################################################################################

# path_functions <- "/home/luna.kuleuven.be/u0141268/github_shared_code_and_publications/Liver_Disease_Microbiome_ML/Functions"
# path_functions <- "/raeslab/scratch/jorvaz/github_shared_code_and_publications/Liver_Disease_Microbiome_ML/Functions" 
path_functions  <- as.character(args[5]) 
source(paste0(path_functions,"/Machine_learning_functions.R"))

#######################################################################################################################################
##########################################                Load the model                         ######################################
#######################################################################################################################################

Fit_model <- readRDS(Model_file)

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

##############################################################################################################################################
##########################################                  Change to the outdir                        ######################################
##############################################################################################################################################
dir.create(OutDir)
setwd(OutDir)

#######################################################################################################################################
##########################################                  Performance stats                    ######################################
#######################################################################################################################################

tempStats <- stats_model_func(Fit.Model = Fit_model, Test.df = MEspDF, RetDF = T, df_complete = MEspDF)
cat("\n","Stats", "\n" )
print(tempStats)
cat("\n")

Ret_stats <- data.frame( Stats_name = rownames(t(tempStats))  , Stats = c(t(tempStats)) )
write.table(Ret_stats,"Validation_stats.tsv",col.names=T,row.names = F,quote=FALSE,sep = "\t")

#######################################################################################################################################
##########################################            Predictions and ROC curve                  ######################################
#######################################################################################################################################

# Predictions
Pred_VAR <- predict(Fit_model, newdata = MEspDF)
names(Pred_VAR) <- rownames(MEspDF)
Pred_Prob <- predict(Fit_model, newdata = MEspDF, type = "prob")
all(names(Pred_VAR) == rownames(Pred_Prob))

# ROC curve (works for 2 or more classes)
n_classes <- length(levels(factor(MEspDF$Variable)))

if (n_classes == 2) {
	# Binary ROC
	roc_curve <- pROC::roc(MEspDF$Variable, Pred_Prob[, 2], levels = levels(MEspDF$Variable), quiet = TRUE)
	auc_value <- as.numeric(pROC::auc(roc_curve))
	mcc_value <- tempStats$MCC
	
	p_roc <- pROC::ggroc(roc_curve, size = 1.2, color = "steelblue") +
		geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "gray50") +
		annotate("text", x = 0.2, y = 0.15, label = paste0("AUC = ", round(auc_value, 3)), size = 5) +
		annotate("text", x = 0.2, y = 0.10, label = paste0("MCC = ", round(mcc_value, 3)), size = 5) +
		theme_bw() +
		labs(title = "ROC Curve", x = "Specificity", y = "Sensitivity")
	
} else {
	# Multiclass ROC (One-vs-Rest)
	classes <- levels(factor(MEspDF$Variable))
	colnames(Pred_Prob) <- classes
	
	roc_list <- list()
	auc_list <- c()
	mcc_list <- c()
	
	for (this_class in classes) {
		# ROC and AUC
		binary_response <- ifelse(MEspDF$Variable == this_class, 1, 0)
		roc_list[[this_class]] <- pROC::roc(binary_response, Pred_Prob[, this_class], quiet = TRUE)
		auc_list[this_class] <- as.numeric(pROC::auc(roc_list[[this_class]]))
		
		# MCC per class (One-vs-Rest)
		binary_pred <- ifelse(Pred_VAR == this_class, 1, 0)
		TP <- sum(binary_response == 1 & binary_pred == 1)
		TN <- sum(binary_response == 0 & binary_pred == 0)
		FP <- sum(binary_response == 0 & binary_pred == 1)
		FN <- sum(binary_response == 1 & binary_pred == 0)
		
		denom_parts <- c((TP + FP), (TP + FN), (TN + FP), (TN + FN))
		if (any(denom_parts == 0)) {
			mcc_list[this_class] <- 0
		} else {
			mcc_list[this_class] <- ((TP * TN) - (FP * FN)) / sqrt(prod(denom_parts))
		}
	}
	
	median_auc <- median(auc_list)
	median_mcc <- median(mcc_list)
	
	# Legend labels with AUC and MCC per class
	legend_labels <- paste0(
		gsub("Enterotype\\.", "", classes),
		" (AUC=", round(auc_list, 2), ", MCC=", round(mcc_list, 2), ")"
	)
	names(legend_labels) <- classes
	
	p_roc <- pROC::ggroc(roc_list, size = 1) +
		geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "gray50") +
		annotate("text", x = 0.3, y = 0.15, 
			label = paste0("Median AUC = ", round(median_auc, 3), "\nMedian MCC = ", round(median_mcc, 3)), 
			size = 4, hjust = 0) +
		theme_bw() +
		labs(title = "ROC Curve (One-vs-Rest)", x = "Specificity", y = "Sensitivity", color = "Class") +
		theme(legend.position = "right") +
		scale_color_brewer(palette = "Set1", labels = legend_labels)
}

ggsave("ROC_curve_validation.pdf", p_roc, width = 12, height = 9)
cat("ROC curve saved to: ROC_curve_validation.pdf\n")


#######################################################################################################################################
##########################################              Confusion Matrix                         ######################################
#######################################################################################################################################

# Create confusion matrix
conf_matrix <- confusionMatrix(Pred_VAR, MEspDF$Variable)

# Extract table
conf_table <- as.data.frame.matrix(conf_matrix$table)

# Add row and column totals
conf_table$Total <- rowSums(conf_table)
conf_table <- rbind(conf_table, Total = colSums(conf_table))

# Save confusion matrix
write.table(conf_table, "Confusion_Matrix_validation_dataset.tsv", 
            col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")

cat("\n=== CONFUSION MATRIX ===\n")
print(conf_matrix$table)
cat("\nConfusion matrix saved to: Confusion_Matrix_validation_dataset.tsv\n\n")

##############################################################################################################################################
##########################################              Save a list with all the stats              ##########################################
##############################################################################################################################################

List_stats_validation_dataset <- list(conf_matrix = conf_matrix, Stats=tempStats, ROC_curve= p_roc, Pred_VAR = Pred_VAR, Pred_Prob = Pred_Prob)
saveRDS(file="List_stats_validation_dataset.rds",List_stats_validation_dataset)

	


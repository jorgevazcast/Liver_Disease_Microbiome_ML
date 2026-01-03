set.seed(12345)
library(caret)
library(phyloseq)
library(microbiome)
library(pROC)
library(ggplot2)

args<-commandArgs(TRUE)

infile <- as.character(args[1]) #     infile <- "./infiles/phylo_holdout_validation_set.rds"
Variable <- as.character(args[2]) #     Variable <- "normal_vs_cirrhosis"   Variable <- "condition"
Model_file <- as.character(args[3]) #     Model_file <- "./RandomForest_FS_TRUE/BEST_MODEL_outerloop_performance.rds"
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

all( names(Pred_VAR) == rownames(Pred_Prob)  )

# ROC curve
roc_curve <- pROC::roc(MEspDF$Variable, Pred_Prob[, 2],  levels = levels(MEspDF$Variable))
auc_value <- as.numeric(pROC::auc(roc_curve))

# Extract MCC from stats
mcc_value <- tempStats$MCC

# Prepare ROC data for ggplot
roc_data <- data.frame(
    sensitivity = roc_curve$sensitivities,
    specificity = roc_curve$specificities
)

# Create ROC plot
p_roc <- ggroc(roc_curve, size = 1.2, color = "steelblue") +
	geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "gray50") +
    annotate("text", x = 0.2, y = 0.15, label = paste0("AUC = ", round(auc_value, 3)), size = 5 ) +
    annotate("text", x = 0.2, y = 0.1, label = paste0("MCC = ", round(mcc_value, 3)), size = 5 ) +
    theme_bw(base_size = 12) +
    labs(title = "ROC Curve - Validation Set",x = "Specificity",y = "Sensitivity") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))



ggsave("ROC_curve_validation.pdf", p_roc, width = 7, height = 7)
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

	


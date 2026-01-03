set.seed(12345)
library(ggplot2)
library(ggpubr)
library(caret)
library(shapviz)
library(fastshap)
library(randomForest)

args<-commandArgs(TRUE)
#library(phyloseq)

# Load your phyloseq object and model
#inphylo <- readRDS("/home/luna.kuleuven.be/u0141268/Dropbox/Jiyeon_Jorge/research/research_liver_diseases/Predcition_ML/infiles/in_phylo_ASV_prop.rds")
outdir <- as.character(args[1]) # indir <- "/home/luna.kuleuven.be/u0141268/Dropbox/Jiyeon_Jorge/research/research_liver_diseases/Predcition_ML/PREDICTIONS//shotgun_Species/normal_vs_liver_cancer/Ret_abundance/RandomForest_FS_TRUE"
infileModel <- as.character(args[2])  #infileModel <- "./RandomForest_FS_TRUE/BEST_MODEL_outerloop_performance.rds"
Bmodel <- readRDS(infileModel)

print(outdir)

#dir.create(outdir)
setwd(outdir)
levels2use <- Bmodel$levels
if(length(levels2use) == 2){

 if(any(grepl("Other",levels2use))){
	Levels <- levels2use[!grepl("Other",levels2use)]
 }else{
	Levels <- levels(factor(levels2use))[1]
 }
}

# -------------------------------------------
# 1. Extract predictors from the model training data
# caret's trainingData contains predictors + ".outcome"
X <- Bmodel$trainingData
y <- X$.outcome  # target variable
X <- X[, !(names(X) %in% ".outcome"), drop = FALSE]

# -------------------------------------------
# 2. Define prediction wrapper for fastshap
#    For classification, we can focus on a specific class (e.g., "Disease")
#    Here, I take the first column of the probability matrix
pred_fun <- function(object, newdata, Cat) {
  predict(object, newdata = newdata, type = "prob")[, Cat]
}

# -------------------------------------------
# 3. Compute SHAP values
# nsim = number of Monte Carlo samples (increase for more precision)
shap_values <- fastshap::explain(
  object = Bmodel$finalModel,
  X = X,
  pred_wrapper = function(obj, newdata) pred_fun(obj, newdata, Cat = Levels),
  nsim = 1000
)

# -------------------------------------------
# 4. Create shapviz object
# X = human-readable feature dataframe
sv <- shapviz(shap_values, X = as.data.frame(X))

# -------------------------------------------
# 5. Visualizations
# Global importance (mean |SHAP|)
Barplot <- sv_importance(sv, kind = "bar")

# Beeswarm plot (shows direction and distribution)
sv_importance(sv, kind = "beeswarm")

# Dependence plot for one feature
sv_dependence(sv, v = names(X)[1])


saveRDS(file="SHAP_object.rds",sv)

if( dim(Barplot$data)[1] < 20){
	ggsave(file = "SHAP_Imp.pdf", Barplot, width = 10, height = 7)
} else if( dim(Barplot$data)[1] > 20){
	ggsave(file = "SHAP_Imp.pdf", Barplot, width = 10, height = 10)
} else if( dim(Barplot$data)[1] > 100){
	ggsave(file = "SHAP_Imp.pdf", p, width = 10, height = 17)
} else if( dim(Barplot$data)[1] > 160){
	ggsave(file = "SHAP_Imp.pdf", p, width = 10, height = 20)
}




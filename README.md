# MLPredictR
Scripts for Implementing different Machine Learning Algorithms in R using the caret Library. The repository is currently under development, and additional algorithms may be included in the future. At present, it features implementations of Random Forest, GLMNet, and XGBoost, along with scripts for running these algorithms on Zeus (specifically for XGBoost).

## Requirements

- [ ] [caret](https://cran.r-project.org/web/packages/caret/index.html)
- [ ] [phyloseq](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html)
- [ ] [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html) 
- [ ] [Boruta](https://cran.r-project.org/web/packages/Boruta/index.html) 
- [ ] [xgboost](https://cran.r-project.org/web/packages/xgboost/index.html) 
- [ ] [pROC](https://cran.r-project.org/web/packages/pROC/index.html) 
- [ ] [microbiome](https://www.bioconductor.org/packages/release/bioc/html/microbiome.html) 
- [ ] [cvAUC](https://cran.r-project.org/web/packages/cvAUC/index.html) 
- [ ] [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) 
- [ ] [ggpubr](https://cran.r-project.org/web/packages/ggpubr/index.html) 
- [ ] [DMwR](https://cran.r-project.org/src/contrib/Archive/DMwR/) 

# Classification
The best model is selected based on the one the maximize the AUC
```
########################################
### Directory of the data functions ####

DIR_FUNCTIONS=$HOME"/github_projects/mlpredictr/Functions" ## MODIFY THIS VARIABLE WHERE YOU HAVE THE REPO

##########################################
### Importance plots and stats script ####

Plot_AUC_IMP=$HOME"/github_projects/mlpredictr/Scripts/Plot_AUC_importance.R"

#################################################
###       Infiles and filter parameters      ####
infile=$HOME"/github_projects/mlpredictr/Examples/Example_phyloseq.RData"  
PrevCutoff=0.2
ncores=15

Variable="gender"

#################################################################################
#######################      Classification examples       ######################
#################################################################################

mkdir Classification
cd Classification

###################
##### XGboost #####

### Model ###
ML_method=$HOME"/github_projects/mlpredictr/Scripts/XGboost.R"
Feature_selection=F ### For XGboost it is not really needed

GridMatrix="xgbGrid_small"  ### "xgbGrid_large",GridMatrix <- "xgbGrid_mid" "xgbGrid_small" #### Optimal xgbGrid_mid but takes a lot of time!!!!!
Rscript --vanilla $ML_method $infile $Variable $PrevCutoff $Feature_selection $ncores $GridMatrix $DIR_FUNCTIONS

### Stats ###
Dir_Results="XGboost_FS_FALSE"
Rscript --vanilla $Plot_AUC_IMP $infile $Variable $PrevCutoff $Dir_Results $DIR_FUNCTIONS


########################
##### RandomForest #####

### Model ###
ML_method=$HOME"/github_projects/mlpredictr/Scripts/RandomForest.R"
Feature_selection=T

Rscript --vanilla $ML_method $infile $Variable $PrevCutoff $Feature_selection $ncores $DIR_FUNCTIONS

### Stats ###
Dir_Results="./RandomForest_FS_TRUE/"
Rscript --vanilla $Plot_AUC_IMP $infile $Variable $PrevCutoff $Dir_Results $DIR_FUNCTIONS

##################
##### glmnet #####

### Model ###
ML_method=$HOME"/github_projects/mlpredictr/Scripts/glmnet.R"
Feature_selection=T ## GLM has its own feature selection methods

Rscript --vanilla $ML_method $infile $Variable $PrevCutoff $Feature_selection $ncores $DIR_FUNCTIONS

### Stats ###
Dir_Results="./glmnet_FS_TRUE/"
Rscript --vanilla $Plot_AUC_IMP $infile $Variable $PrevCutoff $Dir_Results $DIR_FUNCTIONS


```

# Regression
The best model is selected based on the one the minimizes the RSEM

```
########################################
### Directory of the data functions ####

DIR_FUNCTIONS=$HOME"/github_projects/mlpredictr/Functions" ## MODIFY THIS VARIABLE WHERE YOU HAVE THE REPO

##########################################
### Importance plots and stats script ####
Plot_RMSE_importance=$HOME"/github_projects/mlpredictr/Scripts/Plot_RMSE_importance.R"

#################################################
###       Infiles and filter parameters      ####
infile=$HOME"/github_projects/mlpredictr/Examples/Example_phyloseq.RData"  
PrevCutoff=0.2
ncores=15

#################################################################################
#######################      Classification examples       ######################
#################################################################################

mkdir Regression
cd Regression

Variable="moisture"

###################
##### XGboost #####

### Model ###
ML_method=$HOME"/github_projects/mlpredictr/Scripts/XGboost.R"
Feature_selection=F ### For XGboost it is not really needed

GridMatrix="xgbGrid_small"  ### "xgbGrid_large",GridMatrix <- "xgbGrid_mid" "xgbGrid_small" #### Optimal xgbGrid_mid but takes a lot of time!!!!!
Rscript --vanilla $ML_method $infile $Variable $PrevCutoff $Feature_selection $ncores $GridMatrix $DIR_FUNCTIONS

### Stats ###
Dir_Results="XGboost_FS_FALSE"
Rscript --vanilla $Plot_RMSE_importance $infile $Variable $PrevCutoff $Dir_Results $DIR_FUNCTIONS


########################
##### RandomForest #####

### Model ###
ML_method=$HOME"/github_projects/mlpredictr/Scripts/RandomForest.R"
Feature_selection=T

Rscript --vanilla $ML_method $infile $Variable $PrevCutoff $Feature_selection $ncores $DIR_FUNCTIONS

### Stats ###
Dir_Results="./RandomForest_FS_TRUE/"
Rscript --vanilla $Plot_RMSE_importance $infile $Variable $PrevCutoff $Dir_Results $DIR_FUNCTIONS

##################
##### glmnet #####

### Model ###
ML_method=$HOME"/github_projects/mlpredictr/Scripts/glmnet.R"
Feature_selection=T ## GLM has its own feature selection methods

Rscript --vanilla $ML_method $infile $Variable $PrevCutoff $Feature_selection $ncores $DIR_FUNCTIONS

### Stats ###
Dir_Results="./glmnet_FS_TRUE/"
Rscript --vanilla $Plot_RMSE_importance $infile $Variable $PrevCutoff $Dir_Results $DIR_FUNCTIONS

```


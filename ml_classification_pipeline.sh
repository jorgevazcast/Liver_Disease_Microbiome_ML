#!/bin/bash


################################################################################
#                    MACHINE LEARNING PIPELINE - RANDOM FOREST                 #
#                                                                              #
# Description: Trains a Random Forest classifier using nested cross-validation #
#              with feature selection and SMOTE balancing for microbiome data  #
#                                                                              #
# Author: Jorge Francisco Vázquez Castellanos                                  #
# Date: [Date]                                                                 #
################################################################################

### sh ml_classification_pipeline.sh

####################################################################################
#                           GLOBAL VARIABLES                                       #
####################################################################################
mkdir Results
cd Results

# Home directory - use $HOME for current user
HOMEC=$HOME

# Path to the Functions directory containing Machine_learning_functions.R
DIR_FUNCTIONS=$HOMEC"/github_shared_code_and_publications/Liver_Disease_Microbiome_ML/Functions"
# NOTE: MODIFY THIS PATH to match your local repository location

# Prevalence cutoff for feature filtering (0-1)
# Features present in fewer than this proportion of samples will be removed
PrevCutoff=0.2  # 20% prevalence threshold

# Number of CPU cores for parallel processing
ncores=15

# Size for the hold-out partition
hold_out_size=0.1

# Feature selection flag
# TRUE = Apply Boruta feature selection within each CV fold
# FALSE = Use all features (after prevalence filtering)
Feature_selection=T

####################################################################################
#                           SCRIPT PATHS                                           #
####################################################################################

# Path to the file for the file partitions
infile_partition=$HOMEC"/github_shared_code_and_publications/Liver_Disease_Microbiome_ML/Scripts/infile_partition.R"

# Path to Random Forest training script
ML_method_rf=$HOMEC"/github_shared_code_and_publications/Liver_Disease_Microbiome_ML/Scripts/RandomForest.R"

# Path to the script for the best model estimation
Best_model=$HOMEC"/github_shared_code_and_publications/Liver_Disease_Microbiome_ML/Scripts/Best_model.R"

# Path to prediction script (for external validation)
Predict_ML=$HOMEC"/github_shared_code_and_publications/Liver_Disease_Microbiome_ML/Scripts/Predict_ML.R"

####################################################################################
#                           INPUT DATA                                             #
####################################################################################

# Input phyloseq object (.rds format)
# Should contain microbiome abundance data and metadata with target variable
infile=/home/luna.kuleuven.be/u0141268/Dropbox/Jiyeon_Jorge/research/research_liver_diseases/Meta_analysis/ML_predictions/infiles_case_controls/cirrhosis_species_abd_adj_prop_phyloseq.rds

# Target variable for classification
# Must match a column name in the phyloseq metadata
Variable="normal_vs_cirrhosis"

####################################################################################
#                     Parition dataset                                             #
####################################################################################

Rscript --vanilla $infile_partition $infile $Variable $hold_out_size

####################################################################################
#                        MODEL TRAINING - RANDOM FOREST                            #
####################################################################################

# Train Random Forest model with nested cross-validation
# 
# Arguments:
#   1. $infile             - Input phyloseq object (.rds or .RData)
#   2. $Variable           - Target variable name (from metadata)
#   3. $PrevCutoff         - Prevalence threshold (0-1)
#   4. $Feature_selection  - Enable Boruta feature selection (TRUE/FALSE)
#   5. $ncores             - Number of cores for parallel processing
#   6. $DIR_FUNCTIONS      - Path to Functions directory
#
# Outputs:
#   - list_CV_models.rds              - All trained CV models
#   - stats_model.tsv                 - Performance metrics per fold
#   - Summary_Statistics_Nested_CV.tsv - Aggregated CV performance
#   - BEST_MODEL_FULL.rds             - Final model (no SMOTE)
#   - BEST_MODEL_FULL_BALANCED.rds    - Final model (with SMOTE if imbalanced)
#   - Metrics_CV_Boxplots.pdf         - Performance visualization
#   - data_balance_report.txt         - Class balance information

echo "================================"
echo "Starting Random Forest Training"
echo "================================"
echo "Input file: $infile"
echo "Variable: $Variable"
echo "Prevalence cutoff: $PrevCutoff"
echo "Feature selection: $Feature_selection"
echo "Cores: $ncores"
echo "================================"

Rscript --vanilla $ML_method_rf \
    ./infiles/phylo_training_set.rds \
    $Variable \
    $PrevCutoff \
    $Feature_selection \
    $ncores \
    $DIR_FUNCTIONS


###########################################################################
#                     Train the best model                                #
###########################################################################

echo "================================"
echo "Train the best model"
echo "================================"
echo "Input file: $infile"
echo "Variable: $Variable"
echo "Prevalence cutoff: $PrevCutoff"
echo "Feature selection: $Feature_selection"
echo "Cores: $ncores"
echo "================================"

Rscript --vanilla $Best_model \
    ./infiles/phylo_training_set.rds \
    $Variable \
    $PrevCutoff \
    $ncores \
    $DIR_FUNCTIONS \
    ./RandomForest_FS_TRUE \
    "rf"
    


####################################################################################
#                     OPTIONAL: EXTERNAL VALIDATION                                #
####################################################################################

external_dataset="./infiles/phylo_holdout_validation_set.rds"
# 
Rscript --vanilla $Predict_ML \
     $external_dataset \
     $Variable \
     "./RandomForest_FS_TRUE/BEST_MODEL_outerloop_performance.rds" \
     "BEST_MODEL_outerloop_performance_stats"\
     $DIR_FUNCTIONS

Rscript --vanilla $Predict_ML \
     $external_dataset \
     $Variable \
     "./RandomForest_FS_TRUE/BEST_MODEL_NestCV.rds" \
     "BEST_MODEL_NestCV_stats"\
     $DIR_FUNCTIONS
     
Rscript --vanilla $Predict_ML \
     $external_dataset \
     $Variable \
     "./RandomForest_FS_TRUE/BEST_MODEL_NestCV_common_features.rds" \
     "BEST_MODEL_NestCV_common_features_stats"\
     $DIR_FUNCTIONS
          

################################################################################
#                              END OF SCRIPT                                    #
################################################################################






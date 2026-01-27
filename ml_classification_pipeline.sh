#!/bin/bash
set -e  # Exit on error

################################################################################
#                    MACHINE LEARNING PIPELINE - RANDOM FOREST                 #
#                                                                              #
# Trains a Random Forest classifier using nested cross-validation with         #
# Boruta feature selection and SMOTE balancing for microbiome data.            #
#                                                                              #
# Usage: bash ml_classification_pipeline.sh                                    #
# Author: Jorge Francisco Vázquez Castellanos                                  #
################################################################################

####################################################################################
#                              SETUP & PARAMETERS                                  #
####################################################################################

# Working directory
HOMEC=$(pwd)

# Scripts paths
infile_partition="$HOMEC/Scripts/infile_partition.R"
ML_method_rf="$HOMEC/Scripts/RandomForest.R"
Best_model="$HOMEC/Scripts/Best_model.R"
Predict_ML="$HOMEC/Scripts/Predict_ML.R"
SHAP_biomarker_importance="$HOMEC/Scripts/SHAP_biomarker_importance.R"
DIR_FUNCTIONS="$HOMEC/Functions"

# Verify required files exist
for script in "$infile_partition" "$ML_method_rf" "$Best_model" "$Predict_ML" "$SHAP_biomarker_importance"; do
    [ -f "$script" ] || { echo "Error: $script not found"; exit 1; }
done
[ -f "$DIR_FUNCTIONS/Machine_learning_functions.R" ] || { echo "Error: Machine_learning_functions.R not found"; exit 1; }

# Input data
infile="$HOMEC/Test_data/in_phylo.rds"
[ -f "$infile" ] || { echo "Error: $infile not found"; exit 1; }

# Model parameters
Variable="condition"          # Target variable (column in metadata)
PrevCutoff=0.2                # Prevalence filter (20%)
ncores=15                     # CPU cores for model training
hold_out_size=0.3             # Hold-out set proportion (30%)
Feature_selection=TRUE        # Boruta feature selection (TRUE/FALSE)

# Output directories
mkdir -p Results
cd Results
MODEL_DIR="./RandomForest_FS_${Feature_selection}"

####################################################################################
#                         1. DATA PARTITION                                        #
####################################################################################

echo "=== 1. Partitioning data ==="
Rscript --vanilla "$infile_partition" "$infile" "$Variable" "$hold_out_size"

####################################################################################
#                         2. MODEL TRAINING                                        #
####################################################################################

echo "=== 2. Training Random Forest ==="
Rscript --vanilla "$ML_method_rf" \
    ./infiles/phylo_training_set.rds \
    "$Variable" \
    "$PrevCutoff" \
    "$Feature_selection" \
    "$ncores" \
    "$DIR_FUNCTIONS"

####################################################################################
#                         3. BEST MODEL SELECTION                                  #
####################################################################################

echo "=== 3. Selecting best model ==="
Rscript --vanilla "$Best_model" \
    ./infiles/phylo_training_set.rds \
    "$Variable" \
    "$PrevCutoff" \
    "$ncores" \
    "$DIR_FUNCTIONS" \
    "$MODEL_DIR" \
    "rf"

####################################################################################
#                         4. EXTERNAL VALIDATION                                   #
####################################################################################

echo "=== 4. External validation ==="
external_dataset="./infiles/phylo_holdout_validation_set.rds"
models=("BEST_MODEL_outerloop_performance" "BEST_MODEL_NestCV" "BEST_MODEL_NestCV_common_features")

for model in "${models[@]}"; do
    Rscript --vanilla "$Predict_ML" \
        "$external_dataset" \
        "$Variable" \
        "${MODEL_DIR}/${model}.rds" \
        "${model}_stats" \
        "$DIR_FUNCTIONS"
done

####################################################################################
#                         5. SHAP ANALYSIS (PARALLEL)                              #
####################################################################################

echo "=== 5. SHAP feature importance (parallel) ==="
Ncores_shap=3

# Build combinations: outdir|inmodel
combinations=()
for model in "${models[@]}"; do
    combinations+=("./${model}_stats|${MODEL_DIR}/${model}.rds")
done

# Run in parallel
printf "%s\n" "${combinations[@]}" | xargs -P "$Ncores_shap" -I {} bash -c '
    outdir="${1%|*}"
    inmodel="${1#*|}"
    echo "Processing: $inmodel"
    Rscript --vanilla '"$SHAP_biomarker_importance"' "$outdir" "$inmodel"
' _ {}

echo "=== Pipeline completed ==="

# MLPredictR
Scripts for implementing Random Forest classification in R using the caret library with nested cross-validation, Boruta feature selection, and SHAP-based feature importance.

## Pipeline workflow
![ML Classification Pipeline](images/Raw_pipeline.png)
*Figure 1: Nested cross-validation framework for microbiome classification*

## Requirements

- [caret](https://cran.r-project.org/web/packages/caret/index.html)
- [coin](https://cran.r-project.org/web/packages/coin/index.html/) 
- [phyloseq](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html)
- [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html) 
- [Boruta](https://cran.r-project.org/web/packages/Boruta/index.html) 
- [pROC](https://cran.r-project.org/web/packages/pROC/index.html) 
- [microbiome](https://www.bioconductor.org/packages/release/bioc/html/microbiome.html) 
- [cvAUC](https://cran.r-project.org/web/packages/cvAUC/index.html) 
- [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) 
- [ggpubr](https://cran.r-project.org/web/packages/ggpubr/index.html) 
- [DMwR](https://cran.r-project.org/src/contrib/Archive/DMwR/) 
- [fastshap](https://cran.r-project.org/web/packages/fastshap/index.html/) 
- [rstatix](https://cran.r-project.org/web/packages/rstatix/index.html/) 
- [shapviz](https://cran.r-project.org/web/packages/shapviz/index.html/) 

## Quick Start

### 1. Clone the repository
```bash
git clone https://github.com/jorgevazcast/Liver_Disease_Microbiome_ML.git
cd Liver_Disease_Microbiome_ML
```

### 2. Configure the pipeline

Edit `ml_classification_pipeline.sh` and modify these parameters:
```bash
# Input data
infile="$HOMEC/Test_data/in_phylo.rds"  # Your phyloseq object

# Model parameters
Variable="condition"          # Target variable (column in metadata)
PrevCutoff=0.2                # Prevalence filter (20%)
ncores=15                     # CPU cores for model training
hold_out_size=0.3             # Hold-out set proportion (30%)
Feature_selection=TRUE        # Boruta feature selection (TRUE/FALSE)
```

### 3. Run the pipeline
```bash
bash ml_classification_pipeline.sh
```

## Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 1 | `infile_partition.R` | Splits data into training and hold-out sets |
| 2 | `RandomForest.R` | Trains RF with nested CV + optional Boruta feature selection |
| 3 | `Best_model.R` | Selects and retrains best model on full training set |
| 4 | `Predict_ML.R` | Validates models on hold-out set |
| 5 | `SHAP_biomarker_importance.R` | Computes SHAP values for feature importance |

## Output Structure
```
Results/
├── infiles/
│   ├── phylo_training_set.rds
│   └── phylo_holdout_validation_set.rds
├── RandomForest_FS_TRUE/
│   ├── list_CV_models.rds                    # All CV fold models
│   ├── stats_model.tsv                       # Per-fold metrics
│   ├── Summary_Statistics_Nested_CV.tsv      # Aggregated CV stats
│   ├── BEST_MODEL_NestCV.rds                 # Best model from nested CV
│   ├── BEST_MODEL_NestCV_common_features.rds # Best model with common features
│   ├── BEST_MODEL_outerloop_performance.rds  # Best model from outer loop
│   ├── Metrics_CV_Boxplots.pdf               # Performance plots
│   └── data_balance_report.txt               # Class balance info
├── BEST_MODEL_*_stats/
│   └── SHAP_*.pdf                            # SHAP importance plots
└── *.tsv                                     # Validation results
```

## Input Data Format

The input must be a **phyloseq object** (.rds) containing:

- `otu_table`: Feature abundance matrix
- `sample_data`: Metadata with the target variable
- `tax_table`: (optional) Taxonomy information

## Performance Metrics

- AUC, Accuracy, Kappa, MCC, F1-score
- Sensitivity, Specificity, Precision
- PPV, NPV

## Citation

If you use this pipeline, please cite:
```
[Add citation here]
```

## License

MIT License - Copyright (c) 2025 Jorge Francisco Vázquez Castellanos

See [LICENSE](LICENSE) for details.

## Contact

Jorge Francisco Vázquez Castellanos

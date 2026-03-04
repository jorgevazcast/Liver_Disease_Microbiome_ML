set.seed(12345)
library(caret)
library(phyloseq)
library(microbiome)

args<-commandArgs(TRUE)

infile <- as.character(args[1])   # infile   <- "/Liver_Disease_Microbiome_ML/Test_data/in_phylo_LOCO_CV.rds"
Variable <- as.character(args[2]) # Variable <- "condition"
Cohort <-  as.character(args[3])  # Cohort <-  "cohort_1"


# Load phyloseq and extract metadata
all_phylo <- readRDS(infile)
Metadata <- sample_data(all_phylo)
Metadata$Variable <- Metadata[[Variable]]

# Remove samples with missing target variable
Metadata <- Metadata[!is.na(Metadata$Variable), ]
all_phylo <- prune_samples(rownames(Metadata), all_phylo)

# Remove taxa with zero abundance across all samples
all_taxa <- taxa_sums(all_phylo)
all_phylo <- prune_taxa(names(all_taxa[all_taxa != 0]), all_phylo)

# Subset the cohort
training_phylo <- subset_samples(all_phylo, cohort!=Cohort)
holdout_phylo <- subset_samples(all_phylo, cohort==Cohort)

# Remove taxa with zero abundance in each set separately
# (some taxa may become zero after split)
all_taxa <- taxa_sums(training_phylo)
#training_phylo <- prune_taxa(names(all_taxa[all_taxa != 0]), training_phylo)

all_taxa <- taxa_sums(holdout_phylo)
#holdout_phylo <- prune_taxa(names(all_taxa[all_taxa != 0]), holdout_phylo)

# Create output directory and save
dir.create("infiles", showWarnings = FALSE)
setwd("infiles")

saveRDS(training_phylo, "phylo_training_set.rds")
saveRDS(holdout_phylo, "phylo_holdout_validation_set.rds")

cat("Training set:", nsamples(training_phylo), "samples\n")
cat("Hold-out set:", nsamples(holdout_phylo), "samples\n")


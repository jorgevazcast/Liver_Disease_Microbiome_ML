set.seed(12345)
library(caret)
library(phyloseq)
library(microbiome)

args<-commandArgs(TRUE)

infile <- as.character(args[1]) #     infile <- "/home/luna.kuleuven.be/u0141268/Dropbox/Jiyeon_Jorge/research/research_liver_diseases/Meta_analysis/ML_predictions/infiles_case_controls/cirrhosis_species_abd_adj_prop_phyloseq.rds"
Variable <- as.character(args[2]) #     Variable <- "normal_vs_cirrhosis"
hold_out_size <- as.numeric(as.character(args[3])) # hold_out_size <- 0.2



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

# Create stratified split
holdout_index <- createDataPartition(
    y = Metadata$Variable, 
    p = hold_out_size, 
    list = FALSE
)

training_phylo <- prune_samples(rownames(Metadata)[-holdout_index], all_phylo)
holdout_phylo <- prune_samples(rownames(Metadata)[holdout_index], all_phylo)

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


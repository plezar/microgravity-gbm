library(data.table)
library(MOFA2)
library(tidyverse)
source("util.R")
set.seed(123)

#### 1. Parameters ####
# Samples per condition (can differ by modality)
n_rna_A  <- 20
n_rna_B  <- 20
n_prot_A <- 15
n_prot_B <- 15

# Features per modality
n_genes <- 1000
n_prots <- 200

# Number of truly affected features (condition effect)
n_de_genes <- 100
n_de_prots <- 40

#### 2. Sample metadata (note: non-overlapping sample IDs) ####
rna_samples <- data.frame(
  sample_id = c(paste0("RNA_A_", seq_len(n_rna_A)),
                paste0("RNA_B_", seq_len(n_rna_B))),
  condition = rep(c("A", "B"), times = c(n_rna_A, n_rna_B))
)

prot_samples <- data.frame(
  sample_id = c(paste0("PROT_A_", seq_len(n_prot_A)),
                paste0("PROT_B_", seq_len(n_prot_B))),
  condition = rep(c("A", "B"), times = c(n_prot_A, n_prot_B))
)

#### 3. Simulate RNA-seq–like data ####
# Base expression matrix: genes x samples
rna_mat <- matrix(
  rnorm(n_genes * nrow(rna_samples), mean = 0, sd = 1),
  nrow = n_genes,
  ncol = nrow(rna_samples)
)
rownames(rna_mat) <- paste0("Gene", seq_len(n_genes))
colnames(rna_mat) <- rna_samples$sample_id

# Add a condition effect to first n_de_genes in condition B
de_genes <- seq_len(n_de_genes)
is_B_rna <- rna_samples$condition == "B"
rna_mat[de_genes, is_B_rna] <- rna_mat[de_genes, is_B_rna] + 1.5  # log2 fold-change-ish

#### 4. Simulate proteomics data ####
prot_mat <- matrix(
  rnorm(n_prots * nrow(prot_samples), mean = 0, sd = 1),
  nrow = n_prots,
  ncol = nrow(prot_samples)
)
rownames(prot_mat) <- paste0("Prot", seq_len(n_prots))
colnames(prot_mat) <- prot_samples$sample_id

# Add a condition effect to first n_de_prots in condition B
de_prots <- seq_len(n_de_prots)
is_B_prot <- prot_samples$condition == "B"
prot_mat[de_prots, is_B_prot] <- prot_mat[de_prots, is_B_prot] + 1.0

#### 5. Wrap into a simple list structure ####
multi_omics <- list(
  rna = list(
    X = rna_mat,
    samples = rna_samples
  ),
  proteomics = list(
    X = prot_mat,
    samples = prot_samples
  )
)

# ====MOFA====

MOFAobject <- create_mofa(dt)

#plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 5

train_opts <- get_default_training_options(MOFAobject)

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

outfile = file.path(tempdir(),"model.hdf5")

library(reticulate)
py_config()  # just to confirm you're using the uv env you showed
py_install("mofapy2", pip = TRUE)

MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk=F)


# ===downstream analysis=====

## ====add metadata======


samples_metadata(MOFAobject.trained) <- meta_dta %>% rename(sample = Sample) #, group = Patient.code
head(MOFAobject.trained@samples_metadata, n=3)

head(MOFAobject.trained@cache$variance_explained$r2_tota) # group 1

head(MOFAobject.trained@cache$variance_explained$r2_per_factor) # group 1

plot_variance_explained(MOFAobject.trained, x="view", y="factor")

colnames(MOFAobject.trained@samples_metadata)

plot_factors(MOFAobject.trained, 
             factors = 1:3,
             color_by = "IDH.status"
)

plot_factors(MOFAobject.trained, 
             factors = 1:3,
             color_by = "gsea.subtype.call"
)

## ====umap=====

library(uwot)

factors <- get_factors(MOFAobject.trained, factors = "all")



X <- factors$single_group   # samples × factors matrix

set.seed(123)
umap_out <- umap(X, n_neighbors = 5, min_dist = 0.3, metric = "euclidean")

# Convert to data.frame
umap_df <- as.data.frame(umap_out)
colnames(umap_df) <- c("UMAP1", "UMAP2")

ggplot(umap_df, aes(UMAP1, UMAP2)) +
  geom_point(size = 1) +
  theme_classic()


## ====clustering====

factors <- get_factors(MOFAobject.trained, factors = "all")
X <- factors$single_group

library(ComplexHeatmap)

library(ConsensusClusterPlus)

results <- ConsensusClusterPlus(t(X),
                                maxK=10,
                                reps=1000,
                                distance="euclidean",
                                title="/Users/mzarodniuk/Documents/Scripts/MREOmics/omics/results/ConsensusClusterPlus",
                                seed=42,
                                plot="pdf")

clustering <- results[[4]][["consensusClass"]]
Heatmap(X)


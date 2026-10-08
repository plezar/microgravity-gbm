library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(GSVA)

# GBM ---------

# build a query to get gene expression data for entire cohort
query_gbm_all = GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor",
  access = "open")

# download data
TCGA_DIR <- "01_bulk-mrnaseq/data/TCGA"
GDCdownload(query_gbm_all, directory=TCGA_DIR)

# get counts
tcga_gbm_data <- GDCprepare(query_gbm_all, directory = TCGA_DIR, summarizedExperiment = TRUE)
gbm_matrix <- assay(tcga_gbm_data, "unstranded")
dim(gbm_matrix)

rownames(gbm_matrix) <- gene_metadata$gene_name
gbm_matrix <- rowsum(gbm_matrix, group = rownames(gbm_matrix), reorder = F)

# extract gene and sample metadata from summarizedExperiment object
gene_metadata <- as.data.frame(rowData(tcga_gbm_data))
coldata <- as.data.frame(colData(tcga_gbm_data))

# add survival stats from cBioPortal
cbio_meta <- read_tsv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/data/TCGA/data_clinical_patient.txt", comment = "#")

coldata <- coldata %>% left_join(
  cbio_meta, by = c("patient"="PATIENT_ID")
)


## Normalization ----------

# Setting up countData object   
dds <- DESeqDataSetFromMatrix(countData = gbm_matrix,
                              colData = coldata,
                              design = ~ 1)

#rowData(dds) <- gene_metadata
saveRDS(dds, "01_bulk-mrnaseq/data/TCGA/TCGA-GBM.rds")

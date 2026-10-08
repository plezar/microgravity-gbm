#!/usr/bin/env Rscript

# example use 
# Rscript RNA/02_Enrichment_Analysis.R -i results/RNA/deseq2_wt_stiff_vs_soft.csv -o results/RNA/Enrichment/WT/StiffSoft/
# Rscript RNA/02_Enrichment_Analysis.R -i results/RNA/deseq2_wt_stiff_vs_base.csv -o results/RNA/Enrichment/WT/StiffBase/
# Rscript RNA/02_Enrichment_Analysis.R -i results/RNA/deseq2_wt_soft_vs_base.csv -o results/RNA/Enrichment/WT/SoftBase/

# Rscript RNA/02_Enrichment_Analysis.R -i results/RNA/deseq2_mut_stiff_vs_soft.csv -o results/RNA/Enrichment/MUT/StiffSoft/
# Rscript RNA/02_Enrichment_Analysis.R -i results/RNA/deseq2_mut_stiff_vs_base.csv -o results/RNA/Enrichment/MUT/StiffBase/
# Rscript RNA/02_Enrichment_Analysis.R -i results/RNA/deseq2_mut_soft_vs_base.csv -o results/RNA/Enrichment/MUT/SoftBase/

# Libraries -----------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(msigdbr)
  library(tidyverse)
})

# Cmd line options -----------------------------------
option_list <- list(
  make_option(c("-i", "--input"),
              type    = "character",
              help    = "Path to input CSV (with columns: baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, symbol)",
              metavar = "file"),
  make_option(c("-o", "--outdir"),
              type    = "character",
              help    = "Directory to write ORA/GSEA results into",
              metavar = "dir"),
  make_option(c("-p", "--padj"),
              type    = "double",
              default = 0.05,
              help    = "Adjusted-p cutoff for selecting DE genes [%default]",
              metavar = "double")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$outdir)) {
  print_help(opt_parser)
  stop("Both --input and --outdir must be supplied.\n", call. = FALSE)
}

input_csv <- opt$input
out_dir    <- opt$outdir
padj_cut   <- opt$padj

## DEBUG
#input_csv <- "/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/results/deseq2_uG_vs_KSC_U87.csv"
#out_dir <- "/Users/mzarodniuk/Documents/Scripts/MREomics/results"

# Read data -----------------------------------
df <- read.csv(input_csv, row.names = 1, stringsAsFactors = FALSE)
# expecting columns: baseMean, log2FoldChange, lfcSE, stat, pvalue, padj

# Convert rownames (SYMBOL) to Entrez
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = rownames(df),
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first"
)

# Add to dataframe
df$ENTREZID <- entrez_ids
df <- df[!is.na(df$ENTREZID),]

# DEFINE GENE LISTS AND UNIVERSE
# all tested genes (background)
count_mtx <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/data/count_mtx.csv", row.names = 1)
gene_universe <- mapIds(
  org.Hs.eg.db,
  keys = rownames(count_mtx),
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first"
)

# genes of interest (here: those with padj < 0.05; adjust threshold as needed)
gene_entrez   <- df %>% filter(padj < 0.05) %>% pull("ENTREZID")
gene_entrez_up   <- df %>% filter(padj < 0.05, log2FoldChange > 0) %>% pull("ENTREZID")
gene_entrez_down   <- df %>% filter(padj < 0.05, log2FoldChange < 0) %>% pull("ENTREZID")

# prepare ranked list for GSEA
gene_list <- df$log2FoldChange
names(gene_list) <- df$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# PREP OUTPUT DIRS
ora_dir  <- file.path(out_dir, "ORA")
gsea_dir <- file.path(out_dir, "GSEA")
dir.create(ora_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)

# ORA: GO upregulated (BP/CC/MF) -----------------------------------
go_ora_up <- enrichGO(gene          = gene_entrez_up,
                      universe      = gene_universe,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2)
saveRDS(go_ora_up, file.path(ora_dir, "GO_ORA_UP_ALL.rds"))
write.csv(as.data.frame(go_ora_up),
          file = file.path(ora_dir, "GO_ORA_UP_ALL.csv"),
          row.names = FALSE)

# ORA: GO downregulated (BP/CC/MF) -----------------------------------
go_ora_down <- enrichGO(gene          = gene_entrez_down,
                        universe      = gene_universe,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = "ENTREZID",
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.2)
saveRDS(go_ora_down, file.path(ora_dir, "GO_ORA_DOWN_ALL.rds"))
write.csv(as.data.frame(go_ora_down),
          file = file.path(ora_dir, "GO_ORA_DOWN_ALL.csv"),
          row.names = FALSE)

# ORA: KEGG -----------------------------------
kegg_ora <- enrichKEGG(gene          = gene_entrez,
                       universe      = gene_universe,
                       organism      = "hsa",
                       keyType       = "ncbi-geneid",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.2)
write.csv(as.data.frame(kegg_ora),
          file = file.path(ora_dir, "KEGG_ORA.csv"),
          row.names = FALSE)

# ORA: Reactome -----------------------------------
react_ora_up <- enrichPathway(gene          = gene_entrez_up,
                           universe      = gene_universe,
                           organism      = "human",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.2,
                           readable      = TRUE)
saveRDS(go_ora_up, file.path(ora_dir, "Reactome_ORA_UP_ALL.rds"))
write.csv(as.data.frame(react_ora_up),
          file = file.path(ora_dir, "Reactome_ORA_UP_ALL.csv"),
          row.names = FALSE)

react_ora_down <- enrichPathway(gene          = gene_entrez_down,
                              universe      = gene_universe,
                              organism      = "human",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05,
                              qvalueCutoff  = 0.2,
                              readable      = TRUE)
saveRDS(go_ora_down, file.path(ora_dir, "Reactome_ORA_DOWN_ALL.rds"))
write.csv(as.data.frame(react_ora_down),
          file = file.path(ora_dir, "Reactome_ORA_DOWN_ALL.csv"),
          row.names = FALSE)

# ORA: Hallmark (MSigDB H) -----------------------------------
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, entrez_gene) %>%
  distinct()

hallmark_ora <- enricher(gene          = gene_entrez,
                         universe      = gene_universe,
                         TERM2GENE     = hallmark_sets,
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.2)
write.csv(as.data.frame(hallmark_ora),
          file = file.path(ora_dir, "Hallmark_ORA.csv"),
          row.names = FALSE)

# GSEA: GO-BP -----------------------------------
go_gsea <- gseGO(geneList      = gene_list,
                 OrgDb         = org.Hs.eg.db,
                 ont           = "ALL",
                 keyType       = "ENTREZID",
                 nPerm         = 1000,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1)
saveRDS(go_gsea, file.path(gsea_dir, "GO_GSEA_ALL.rds"))
write.csv(as.data.frame(go_gsea),
          file = file.path(gsea_dir, "GO_GSEA_ALL.csv"),
          row.names = FALSE)

# GSEA: Hallmark -----------------------------------
hallmark_gsea <- GSEA(geneList      = gene_list,
                      TERM2GENE     = hallmark_sets,
                      pAdjustMethod = "BH",
                      nPerm         = 1000,
                      pvalueCutoff  = 0.05)
saveRDS(hallmark_gsea, file.path(gsea_dir, "Hallmark_GSEA.rds"))
write.csv(as.data.frame(hallmark_gsea),
          file = file.path(gsea_dir, "Hallmark_GSEA.csv"),
          row.names = FALSE)

message("All analyses complete. Results saved in:\n ",
        ora_dir, "\n ",
        gsea_dir, "\n")


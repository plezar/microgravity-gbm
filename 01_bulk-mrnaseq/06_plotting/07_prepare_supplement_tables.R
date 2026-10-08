#!/usr/bin/env Rscript

# ============================================================
# Write 4 DESeq2 result tables into 1 Excel file (4 worksheets)
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(writexl)
})

# ---- Input files ----
base_dir <- "/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/results"

files <- list(
  "U87THP_vs_U87_KSC"     = file.path(base_dir, "deseq2_U87THP_vs_U87_KSC.csv"),
  "U87THP_vs_U87_uG"      = file.path(base_dir, "deseq2_U87THP_vs_U87_uG.csv"),
  "uG_vs_KSC_U87"         = file.path(base_dir, "deseq2_uG_vs_KSC_U87.csv"),
  "uG_vs_KSC_U87THP1"     = file.path(base_dir, "deseq2_uG_vs_KSC_U87THP1.csv")
)

# ---- Read all tables ----
tables <- lapply(files, function(f) {
  if (!file.exists(f)) stop("File not found: ", f)
  read_csv(f, show_col_types = FALSE)
})

# ---- Output path ----
out_xlsx <- file.path(base_dir, "DESeq2_results_all_contrasts.xlsx")

# ---- Write workbook ----
write_xlsx(tables, path = out_xlsx)

message("Wrote Excel workbook: ", out_xlsx)
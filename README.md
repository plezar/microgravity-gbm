# Microgravity GBM analyses

Research scripts and notebooks for analyzing microgravity-associated changes in U87 glioblastoma spheroids and U87–THP-1 cocultures. Analyses cover bulk RNA sequencing, Xenium spatial transcriptomics, Olink proteomics, and image-based mechanical inference.

## Repository layout

| Directory | Contents |
| --- | --- |
| `01_bulk-mrnaseq/` | UMI-based RNA-seq preprocessing, count matrix assembly, DESeq2 differential expression, enrichment and HOMER motif analyses, plotting, and TCGA analyses. |
| `02_st-xenium/` | Xenium-to-Zarr conversion, SpatialData processing, consensus NMF gene expression programs (GEPs), enrichment, RCTD cell-type deconvolution, and spatial statistics. |
| `02_st-xenium/biomechanics/` | Segmentation preparation and TensionMap inference of cell pressure, junction tension, and stress. |
| `03_olink/` | Missing-value imputation, limma differential protein abundance analysis, and volcano plots. |
| `utils/` | Shared R functions for enrichment analysis and plotting. |

## Analysis entry points

### Bulk RNA sequencing

1. Follow the scripts in [`00_bulk_preprocess/`](01_bulk-mrnaseq/00_bulk_preprocess/). The [preprocessing README](01_bulk-mrnaseq/00_bulk_preprocess/readme.md) describes UMI extraction, trimming, alignment, deduplication, and gene counting.
2. Assemble the count matrix with `01_make_count_mtx.R`.
3. Use `03_deseq2.R` or `03_deseq2_collapse_replicates.R` for differential expression. The design includes gravity, cell composition, and their interaction.
4. Explore enrichment with `04_enrichment_analysis.R` and motif analysis with `05_homer/`.
5. Use `06_plotting/` for figures and supplementary tables, and `07_tcga/` for TCGA preparation, BayesPrism deconvolution, and survival analysis.

### Xenium spatial transcriptomics

1. Convert Xenium output to Zarr using `00_convert_zarrr.py` and process SpatialData objects in `01_merge_sdata.ipynb`.
2. Identify GEPs with `02_nmf.ipynb`, then examine enrichment in `03_gep_enrichment.R` and spatial organization in `04_gep_spatial.ipynb`.
3. Prepare inputs with `05_prepare_rctd.ipynb`, run RCTD with `05_run_rctd.R`, and inspect the reference and results with `05_rctd_ref.R` and `05_plot_rctd.R`.
4. The other `05_*` scripts cover GEP scoring, radial spatial analyses, linear models, generalized additive models, and Moran's I.
5. `06_homer/` contains motif analyses; `06_plot_olink_proteins.ipynb` and `07_olink_hits_analysis.R` connect protein hits to spatial expression.

The `biomechanics/` notebooks prepare segmentation masks and run TensionMap. Backend requirements vary; the MATLAB optimizer requires MATLAB and a valid license.

### Olink proteomics

The numbered scripts `01_impute_missing.R`, `02_limma.R`, and `03_plot_volcano.R` cover missing-value imputation, differential abundance, and visualization, respectively.

## Setup and data

Open `microgravity-gbm.Rproj` in RStudio for the R analyses. Install the packages used by the relevant scripts; these include DESeq2, tidyverse, limma, missForest, Seurat, and spacexr. Bulk preprocessing additionally uses UMI-tools, Trimmomatic, STAR, samtools, and featureCounts; motif analyses use HOMER.

The Python environment is recorded in [`xenium_env.yml`](xenium_env.yml). It is a platform-specific Conda export with pinned builds and a local environment prefix, so review it for your platform before recreating it. The preprocessing environment is recorded separately in [`umitools.yml`](01_bulk-mrnaseq/00_bulk_preprocess/umitools.yml).

These are study-specific analysis scripts and exploratory notebooks. Before running an analysis, review its input paths, sample metadata, parameters, required packages, and output directories. Several scripts use absolute paths or depend on objects created in earlier steps; numbering provides an organizational guide rather than an automated end-to-end pipeline.

Raw data, generated results, figures, local caches, and downloaded model artifacts are excluded from version control. Supply the required data and reference files locally in the locations expected by each analysis. Archived scripts are retained in `archive/` directories for reference.

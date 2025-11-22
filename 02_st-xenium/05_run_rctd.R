library(spacexr)
library(Seurat)

# ====Construct Ref. Object========

ref_object <- readRDS("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/data/U87_THP_combined_ref.rds")

ref_object <- JoinLayers(ref_object)
counts <- GetAssayData(ref_object, assay = "RNA", slot = "counts")

cluster <- as.factor(ref_object$cell_type)
names(cluster) <- colnames(ref_object)

nUMI <- ref_object$nCount_RNA
names(nUMI) <- colnames(ref_object)
nUMI <- colSums(counts)

reference <- Reference(counts, cluster, nUMI)


# ====Construct Query Object========

counts <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/data/filtered_xenium_counts.csv", row.names = 1)
counts <- t(counts)

coords <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/data/filtered_xenium_coords.csv", row.names = 1)

query <- SpatialRNA(coords, counts, colSums(counts))

RCTD <- create.RCTD(query, reference, max_cores = 10)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
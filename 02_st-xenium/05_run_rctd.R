library(spacexr)
library(Seurat)

# ====Construct Ref. Object========

ref_object <- readRDS("/users/mzarodn2/afs/Private/microgravity-gbm/02_st-xenium/data/U87_THP_combined_ref.rds")

ref_object <- JoinLayers(ref_object)
counts <- GetAssayData(ref_object, assay = "RNA", slot = "counts")

cluster <- as.factor(ref_object$cell_type)
names(cluster) <- colnames(ref_object)

nUMI <- ref_object$nCount_RNA
names(nUMI) <- colnames(ref_object)
nUMI <- colSums(counts)

reference <- Reference(counts, cluster, nUMI)


# ====Construct Query Object========

counts <- read.csv("/users/mzarodn2/afs/Private/microgravity-gbm/02_st-xenium/data/filtered_xenium_counts.csv", row.names = 1)
counts <- t(counts)

coords <- read.csv("/users/mzarodn2/afs/Private/microgravity-gbm/02_st-xenium/data/filtered_xenium_coords.csv", row.names = 1)

query <- SpatialRNA(coords, counts, colSums(counts))

RCTD <- create.RCTD(query, reference, max_cores = 11)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

annotations.df <- RCTD@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)

table(annotations)
#annotations
#THP1  U87 
#   2 5236

saveRDS(RCTD, "/users/mzarodn2/afs/Private/microgravity-gbm/02_st-xenium/data/RCTD_obj.rds")
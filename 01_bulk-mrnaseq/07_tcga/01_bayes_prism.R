library(BayesPrism)
library(Seurat)

# ===== scRNA-seq ======

seurat_obj <- readRDS("/Volumes/TOSHIBA/Patzke_project/RNA-seq/GSE84465/seurat_obj.rds")
cell.type.labels <- seurat_obj$cell.type
cell.state.labels <- seurat_obj$cell.state
sc.dat <- as.matrix(seurat_obj@assays$RNA$counts)

sc.stat <- plot.scRNA.outlier(
  input=t(sc.dat),
  cell.type.labels=cell.type.labels,
  species="hs",
  return.raw=TRUE
)

filter <- (sc.stat$Rb | sc.stat$Mrp | sc.stat$other_Rb |  sc.stat$chrM | sc.stat$MALAT1 | sc.stat$chrX | sc.stat$chrY | sc.stat$act | sc.stat$hb)
sc.dat.filtered <- sc.dat[!filter,]
sc.dat.filtered.pc <-  select.gene.type (t(sc.dat.filtered), gene.type = "protein_coding")
sc.dat <- t(sc.dat)

diff.exp.stat <- get.exp.stat(sc.dat=sc.dat[,colSums(sc.dat>0)>3],
                              cell.type.labels=as.character(cell.type.labels),
                              cell.state.labels=as.character(cell.state.labels),
                              pseudo.count=0.1,
                              cell.count.cutoff=50,
                              n.cores=9
)

sc.dat.filtered.pc.sig <- select.marker(sc.dat=sc.dat.filtered.pc,
                                        stat=diff.exp.stat,
                                        pval.max=0.01,
                                        lfc.min=2)


# ==== Bulk RNA-seq ====

dds <- readRDS("01_bulk-mrnaseq/data/TCGA/TCGA-GBM.rds")
keep <- rowSums(assay(dds) >= 10) >= 0.05*ncol(dds)
dds <- dds[keep,]



myPrism <- new.prism(
  reference=sc.dat.filtered.pc, 
  mixture=t(assay(dds)),
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.state.labels,
  key="Glioma",
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

bp.res <- run.prism(prism = myPrism, n.cores=8)
saveRDS(bp.res, "/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/results/tcga_bp.rds")


malignant_exp <- get.exp(bp=bp.res,
        state.or.type="type",
        cell.name="Glioma")

malignant_exp_norm <- vst(round(t(malignant_exp)))
saveRDS(malignant_exp_norm, "/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/results/tcga_malignant_expr.rds")

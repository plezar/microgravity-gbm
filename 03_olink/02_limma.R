library(PerformanceAnalytics)
library(tidyverse)
library(limma)
library(edgeR)
library(UpSetR)

# ==== Prepare Data =====

pct_missing <- read.csv("03_olink/data/olink_all_data_imputed_pct_missing.csv", row.names = 1)
olink_all_data_imputed <- read.csv("03_olink/data/olink_all_data_imputed.csv", row.names = 1)

pdf("03_olink/figures/ftr_correlation.pdf", width = 25, height = 15)
chart.Correlation(olink_all_data_imputed)
dev.off()

pdf("03_olink/figures/ftr_correlation_log.pdf", width = 25, height = 15)
chart.Correlation(log10(olink_all_data_imputed))
dev.off()

meta_data <- data.frame(
  Sample = rownames(olink_all_data_imputed),
  Perturbation = str_split(rownames(olink_all_data_imputed), "-", Inf, T)[,1],
  Composition_lvl1 = str_split(rownames(olink_all_data_imputed), "-", Inf, T)[,2],
  Agarose = str_split(rownames(olink_all_data_imputed), "-", Inf, T)[,3]
)

meta_data$Composition_lvl1 <- gsub("250", "", meta_data$Composition_lvl1)
meta_data$Composition_lvl2 <- if_else(meta_data$Agarose=="A", paste0(meta_data$Composition_lvl1, "-", meta_data$Agarose), meta_data$Composition_lvl1)
meta_data$Agarose <- NULL
rownames(meta_data) <- meta_data$Sample

meta_data <- cbind(meta_data, pct_missing)
meta_data$high_imputed <- if_else(meta_data$pct_missing > 0.1, "10pct", "no")


# ===== Limma =====

x <- DGEList(counts = t(olink_all_data_imputed), samples = meta_data)

# legal comparisons only
x <- x[,x$samples$Composition_lvl2!="UM1-A"]

# make sure these are factors
x$samples$Composition_lvl2 <- factor(x$samples$Composition_lvl2)
x$samples$Perturbation     <- factor(x$samples$Perturbation,
                                     levels = c("KSC", "uG"))  # baseline KSC

## ==== MDS =====
par(mfrow=c(2,2))

col.perturb<- x$samples$Perturbation
col.perturb <- rainbow(length(unique(col.perturb)))[match(col.perturb, unique(col.perturb))]

col.comp_lvl1 <- x$samples$Composition_lvl1
col.comp_lvl1 <- rainbow(length(unique(col.comp_lvl1)))[match(col.comp_lvl1, unique(col.comp_lvl1))]

col.comp_lvl2 <- x$samples$Composition_lvl2
col.comp_lvl2 <- rainbow(length(unique(col.comp_lvl2)))[match(col.comp_lvl2, unique(col.comp_lvl2))]

col.high_imputed <- x$samples$high_imputed
col.high_imputed <- rainbow(length(unique(col.high_imputed)))[match(col.high_imputed, unique(col.high_imputed))]

plotMDS(x, labels=x$samples$Perturbation, col=col.perturb)
title(main="A. Microgravity")
plotMDS(x, labels=x$samples$Composition_lvl1, col=col.comp_lvl1)
title(main="B. Composition 1")
plotMDS(x, labels=x$samples$Composition_lvl2, col=col.comp_lvl2)
title(main="C. Composition 2")
plotMDS(x, labels=x$samples$high_imputed, col=col.high_imputed)
title(main="C. Composition 2")

# removing the majority of UM1-A
table(x$samples$high_imputed, x$samples$Composition_lvl2)

# removing samples with high imputed fraction
x <- x[,x$samples$high_imputed=="no"]

## ==== MDS =====
par(mfrow=c(2,2))

col.perturb<- x$samples$Perturbation
col.perturb <- rainbow(length(unique(col.perturb)))[match(col.perturb, unique(col.perturb))]

col.comp_lvl1 <- x$samples$Composition_lvl1
col.comp_lvl1 <- rainbow(length(unique(col.comp_lvl1)))[match(col.comp_lvl1, unique(col.comp_lvl1))]

col.comp_lvl2 <- x$samples$Composition_lvl2
col.comp_lvl2 <- rainbow(length(unique(col.comp_lvl2)))[match(col.comp_lvl2, unique(col.comp_lvl2))]

col.high_imputed <- x$samples$high_imputed
col.high_imputed <- rainbow(length(unique(col.high_imputed)))[match(col.high_imputed, unique(col.high_imputed))]

plotMDS(x, labels=x$samples$Perturbation, col=col.perturb)
title(main="A. Microgravity")
plotMDS(x, labels=x$samples$Composition_lvl1, col=col.comp_lvl1)
title(main="B. Composition 1")
plotMDS(x, labels=x$samples$Composition_lvl2, col=col.comp_lvl2)
title(main="C. Composition 2")
plotMDS(x, labels=x$samples$high_imputed, col=col.high_imputed)
title(main="C. Composition 2")


# ===== All Groups =====

mm <- model.matrix(~ x$samples$Composition_lvl2 + x$samples$Perturbation)
mm_colnames <- colnames(mm)
mm_colnames <- gsub("-", "_", mm_colnames)
mm_colnames <- gsub("x\\$samples\\$Composition_lvl2", "", mm_colnames)
mm_colnames <- gsub("x\\$samples\\$Perturbation", "", mm_colnames)
colnames(mm) <- mm_colnames

y <- voom(x, mm, plot = T)
cont <- makeContrasts(uG_vs_KSC = uG, levels = mm)
vfit <- lmFit(y, mm)
vfit <- contrasts.fit(vfit, cont)
efit <- eBayes(vfit, trend = F)

uG_vs_KSC_all <- topTable(efit,
                          coef      = "uG_vs_KSC",
                          number    = Inf,
                          adjust.method = "BH",
                          sort.by   = "P")

uG_vs_KSC_all


# ===== Each Group Separately =====

get_da_prots <- function(x) {
  mm <- model.matrix(~ x$samples$Perturbation)
  mm_colnames <- colnames(mm)
  mm_colnames <- gsub("-", "_", mm_colnames)
  mm_colnames <- gsub("\\(", "", mm_colnames)
  mm_colnames <- gsub("\\)", "", mm_colnames)
  mm_colnames <- gsub("x\\$samples\\$Perturbation", "", mm_colnames)
  colnames(mm) <- mm_colnames
  y <- voom(x, mm, plot = F)
  cont <- makeContrasts(uG_vs_KSC = uG, levels = mm)
  vfit <- lmFit(y, mm)
  vfit <- contrasts.fit(vfit, cont)
  efit <- eBayes(vfit, trend = F)
  uG_vs_KSC_all <- topTable(efit,
                            coef      = "uG_vs_KSC",
                            number    = Inf,
                            adjust.method = "BH",
                            sort.by   = "P")
  return(uG_vs_KSC_all)
}

# check which comparisons are legal
table(x$samples$Composition_lvl2, x$samples$Perturbation)

da_res_list <- lapply(unique(x$samples$Composition_lvl2), function(c) {
  get_da_prots(x[, x$samples$Composition_lvl2 == c])
})

names(da_res_list) <- unique(x$samples$Composition_lvl2)

# upset plots
up_lists <- lapply(da_res_list, function(x) {
  rownames(x)[x$logFC > 0 & x$adj.P.Val < 0.05]
})
up_input <- UpSetR::fromList(up_lists)

upset(
  up_input,
  nsets = length(up_lists),
  nintersects = 50,
  order.by = "freq"
)

saveRDS(da_res_list, "03_olink/results/da_res_list.rds")

da_res_list <- lapply(da_res_list, function(x) {x$gene <- rownames(x); return(x)})
da_res_list <- bind_rows(da_res_list)
de_proteins <- unique(da_res_list[da_res_list$adj.P.Val<0.05,]$gene)
writeLines(de_proteins, "/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/03_olink/results/olink_da_proteins.txt")

# ==== Interaction model =====

mm <- model.matrix(~ 0 + Composition_lvl2 * Perturbation,
                       data = x$samples)
mm_colnames <- colnames(mm)
mm_colnames <- gsub("-", "_", mm_colnames)
mm_colnames <- gsub(":", "int", mm_colnames)
mm_colnames <- gsub("Composition_lvl2", "", mm_colnames)
mm_colnames <- gsub("Perturbation", "", mm_colnames)
colnames(mm) <- mm_colnames

y    <- voom(x, mm, plot = FALSE)
fit  <- lmFit(y, mm)
efit <- eBayes(fit)

# Extract composition levels
comps <- levels(x$samples$Composition_lvl2)
comps <- gsub("-", "_", comps)

# Create contrast strings
contrast_list <- lapply(comps, function(cmp) {
  # interaction column name
  interaction_term <- paste0(cmp, "intuG")
  # main effect
  main <- "uG"
  
  # if this composition is the baseline, the interaction won't exist
  if (!interaction_term %in% colnames(mm)) {
    return(setNames(main, cmp))
  }
  
  # otherwise: main effect + interaction
  paste0(main, " + ", interaction_term) |> 
    setNames(cmp)
})

contrast_list <- unlist(contrast_list)

cont <- makeContrasts(contrasts = contrast_list, levels = mm)

fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2)

da_res_list_2 <- lapply(contrast_list, function(nm) {
  topTable(fit2,
           coef = nm,
           number = Inf,
           adjust.method = "BH",
           sort.by = "P")
})
names(da_res_list_2) <- names(contrast_list)

up_lists <- lapply(da_res_list_2, function(x) {
  rownames(x)[x$logFC > 0 & x$adj.P.Val < 0.05]
})



# ===== Saturated combined factor model =====

# combined factor: one level per (Composition, Perturbation) combo
x$samples$CompPert <- interaction(
  x$samples$Composition_lvl2,
  x$samples$Perturbation,
  sep = "_"
)

table(x$samples$CompPert)

design <- model.matrix(~ 0 + CompPert, data = x$samples)
colnames(design) <- gsub("-", ".", colnames(design))
colnames(design)


y    <- voom(x, design, plot = FALSE)
fit  <- lmFit(y, design)
fit  <- eBayes(fit)

cont <- makeContrasts(
  UTU_A = CompPertUTU.A_uG - CompPertUTU.A_KSC,
  UTU   = CompPertUTU_uG   - CompPertUTU_KSC,
  UTD_A = CompPertUTD.A_uG - CompPertUTD.A_KSC,
  UM2_A = CompPertUM2.A_uG - CompPertUM2.A_KSC,
  U_A   = CompPertU.A_uG   - CompPertU.A_KSC,
  U     = CompPertU_uG     - CompPertU_KSC,
  levels = design
)

fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2)

da_res_list_3 <- lapply(colnames(cont), function(nm) {
  topTable(fit2,
           coef = nm,
           number = Inf,
           adjust.method = "BH",
           sort.by = "P")
})
names(da_res_list_3) <- colnames(cont)

up_lists <- lapply(da_res_list_3, function(x) {
  rownames(x)[x$logFC > 0 & x$adj.P.Val < 0.05]
})

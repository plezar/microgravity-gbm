library(GSVA)
library(org.Hs.eg.db)
library(tidyverse)
library(survival)
library(patchwork)
library(survminer)
library(cowplot)
library(DESeq2)

# ==== Normalize data ====
dds <- readRDS("01_bulk-mrnaseq/data/TCGA/TCGA-GBM.rds")

#keep <- rowSums(assay(dds) >= 10) >= 0.05*ncol(dds)
#dds <- dds[keep,]
#vsd <- vst(dds)
#vsd <- assay(vsd)
vsd <- readRDS("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/results/tcga_malignant_expr.rds")

# ==== Prepare gene sets ====
gene_set <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/xenium_cNMF_k4_topgenes.csv", row.names = 1)
#gene_set <- as.list(top_genes_nmf)
names(gene_set) <- paste0("GEP", 1:4)

u87_degs <- read.csv("01_bulk-mrnaseq/results/deseq2_uG_vs_KSC_U87.csv", row.names = 1)
u87_degs <- u87_degs %>% arrange(desc(stat)) %>% head(50) %>% rownames()

u87thp_degs <- read.csv("01_bulk-mrnaseq/results/deseq2_uG_vs_KSC_U87THP1.csv", row.names = 1)
u87thp_degs <- u87thp_degs %>% arrange(desc(stat)) %>% head(50) %>% rownames()

gene_set[["GEP_U87"]] <- u87_degs
gene_set[["GEP_U87THP1"]] <- u87thp_degs

# ==== ssGSEA ====
param    <- ssgseaParam(vsd, gene_set)
gsva_out <- gsva(param)
gsva_out <- t(scale(t(gsva_out)))

clinical_cols <- c("patient", "paper_IDH.status", "paper_MGMT.promoter.status", "SEX", "AGE", "OS_STATUS", "OS_MONTHS", "DFS_STATUS", "DFS_MONTHS")
gsva_out <- as.data.frame(t(gsva_out))
gsva_out <- cbind(gsva_out, colData(dds)[,clinical_cols])

patient_id <- "patient"
gsva_patient <- gsva_out %>%
  group_by(.data[[patient_id]]) %>%
  summarise(
    across(starts_with("GEP"), ~mean(.x, na.rm = TRUE), .names = "{.col}"),
    # carry clinical vars (take first non-NA)
    across(all_of(setdiff(clinical_cols, patient_id)), ~dplyr::first(na.omit(.x))),
    .groups = "drop"
  )

gsva_patient$OS_STATUS <- factor(gsva_patient$OS_STATUS, levels = c("0:LIVING", "1:DECEASED"), labels = c(0, 1)) %>% as.vector() %>% as.numeric()
gsva_patient$DFS_STATUS <- factor(gsva_patient$DFS_STATUS, levels = c("0:DiseaseFree", "1:Recurred/Progressed"), labels = c(0, 1)) %>% as.vector() %>% as.numeric()

gsva_patient$GEP1_strata <- ifelse(gsva_patient$GEP1 >= quantile(gsva_patient$GEP1, 0.5, na.rm=TRUE), "HIGH", "LOW")
gsva_patient$GEP2_strata <- ifelse(gsva_patient$GEP2 >= quantile(gsva_patient$GEP2, 0.5, na.rm=TRUE), "HIGH", "LOW")
gsva_patient$GEP3_strata <- ifelse(gsva_patient$GEP3 >= quantile(gsva_patient$GEP3, 0.5, na.rm=TRUE), "HIGH", "LOW")
gsva_patient$GEP4_strata <- ifelse(gsva_patient$GEP4 >= quantile(gsva_patient$GEP4, 0.5, na.rm=TRUE), "HIGH", "LOW")
gsva_patient$GEP_U87_strata <- ifelse(gsva_patient$GEP_U87 >= quantile(gsva_patient$GEP_U87, 0.5, na.rm=TRUE), "HIGH", "LOW")
gsva_patient$GEP_U87THP1_strata <- ifelse(gsva_patient$GEP_U87THP1 >= quantile(gsva_patient$GEP_U87THP1, 0.5, na.rm=TRUE), "HIGH", "LOW")


# ==== KM survival =====

## === OS ====

fit1 <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ GEP1_strata, data = gsva_patient)
p1   <- ggsurvplot(fit1, data = gsva_patient, pval = TRUE, risk.table = TRUE)

fit2 <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ GEP2_strata, data = gsva_patient)
p2   <- ggsurvplot(fit2, data = gsva_patient, pval = TRUE, risk.table = TRUE)

fit3 <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ GEP3_strata, data = gsva_patient)
p3   <- ggsurvplot(fit3, data = gsva_patient, pval = TRUE, risk.table = TRUE)

fit4 <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ GEP4_strata, data = gsva_patient)
p4   <- ggsurvplot(fit4, data = gsva_patient, pval = TRUE, risk.table = TRUE)

km_row <- (p1$plot | p2$plot | p3$plot | p4$plot)
rt_row <- (p1$table | p2$table | p3$table | p4$table)
p_OS <- km_row / rt_row + plot_layout(heights = c(3, 1))

# RNA-seq gene signatures
fit1 <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ GEP_U87_strata, data = gsva_patient)
p1   <- ggsurvplot(fit1, data = gsva_patient, pval = TRUE, risk.table = TRUE)

fit2 <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ GEP_U87THP1_strata, data = gsva_patient)
p2   <- ggsurvplot(fit2, data = gsva_patient, pval = TRUE, risk.table = TRUE)


## === PFS ====

fit1 <- survfit(Surv(DFS_MONTHS, DFS_STATUS) ~ GEP1_strata, data = gsva_patient)
p1   <- ggsurvplot(fit1, data = gsva_patient, pval = TRUE, risk.table = TRUE)

fit2 <- survfit(Surv(DFS_MONTHS, DFS_STATUS) ~ GEP2_strata, data = gsva_patient)
p2   <- ggsurvplot(fit2, data = gsva_patient, pval = TRUE, risk.table = TRUE)

fit3 <- survfit(Surv(DFS_MONTHS, DFS_STATUS) ~ GEP3_strata, data = gsva_patient)
p3   <- ggsurvplot(fit3, data = gsva_patient, pval = TRUE, risk.table = TRUE)

fit4 <- survfit(Surv(DFS_MONTHS, DFS_STATUS) ~ GEP4_strata, data = gsva_patient)
p4   <- ggsurvplot(fit4, data = gsva_patient, pval = TRUE, risk.table = TRUE)

km_row <- (p1$plot | p2$plot | p3$plot | p4$plot)
rt_row <- (p1$table | p2$table | p3$table | p4$table)
p_PFS <- km_row / rt_row + plot_layout(heights = c(3, 1))




# RNA-seq gene signatures
colnames(gsva_patient)[colnames(gsva_patient)=="GEP_U87_strata"] <- "U87_uG"
fit1 <- survfit(Surv(DFS_MONTHS, DFS_STATUS) ~ U87_uG, data = gsva_patient)
p1   <- ggsurvplot(fit1, data = gsva_patient,
                   pval = TRUE,
                   risk.table = TRUE,
                   surv.median.line="hv", conf.int=T,
                   palette = c("black", "grey"),
                   pval.coord = c(50, 0.75))

km_row <- p1$plot + xlab("Time (months)") + theme_cowplot(14) + theme(axis.line = element_blank(),
                                                  axis.title = element_blank(),
                                                  axis.text.x = element_blank(),
                                                  panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
rt_row <- p1$table + xlab("Time (months)") + ggtitle(NULL) + theme_cowplot(14) + theme(axis.line = element_blank(),
                                                                                       axis.title = element_text(size = 12),
                                                                                       panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

p_PFS <- km_row / rt_row + plot_layout(heights = c(3, 1))
ggsave(plot = p_PFS,
       path="01_bulk-mrnaseq/results/figures/",
       filename="KM_PFS_TCGA.svg",
       device = svglite::svglite,
       width=8, height=4, dpi=700)


fit2 <- survfit(Surv(DFS_MONTHS, DFS_STATUS) ~ GEP_U87THP1_strata, data = gsva_patient)
p2   <- ggsurvplot(fit2, data = gsva_patient, pval = TRUE, risk.table = TRUE)



# ==== CoxPH ====

## ==== OS =====

## --- GEP1 ---
fit_cox_gep1 <- coxph(
  Surv(OS_MONTHS, OS_STATUS) ~ GEP1_strata + AGE + SEX + paper_MGMT.promoter.status + paper_IDH.status,
  data = gsva_patient
)
df_gep1 <- droplevels(model.frame(fit_cox_gep1))
p_forest_gep1 <- ggforest(fit_cox_gep1, data = df_gep1, fontsize = 0.9)

## --- GEP2 ---
fit_cox_gep2 <- coxph(
  Surv(OS_MONTHS, OS_STATUS) ~ GEP2_strata + AGE + SEX + paper_MGMT.promoter.status + paper_IDH.status,
  data = gsva_patient
)
df_gep2 <- droplevels(model.frame(fit_cox_gep2))
p_forest_gep2 <- ggforest(fit_cox_gep2, data = df_gep2, fontsize = 0.9)

## --- GEP3 ---
fit_cox_gep3 <- coxph(
  Surv(OS_MONTHS, OS_STATUS) ~ GEP3_strata + AGE + SEX + paper_MGMT.promoter.status + paper_IDH.status,
  data = gsva_patient
)
df_gep3 <- droplevels(model.frame(fit_cox_gep3))
p_forest_gep3 <- ggforest(fit_cox_gep3, data = df_gep3, fontsize = 0.9)

## --- GEP4 ---
fit_cox_gep4 <- coxph(
  Surv(OS_MONTHS, OS_STATUS) ~ GEP4_strata + AGE + SEX + paper_MGMT.promoter.status + paper_IDH.status,
  data = gsva_patient
)
df_gep4 <- droplevels(model.frame(fit_cox_gep4))
p_forest_gep4 <- ggforest(fit_cox_gep4, data = df_gep4, fontsize = 0.9)

## --- GBM ---
fit_cox_gep5 <- coxph(
  Surv(OS_MONTHS, OS_STATUS) ~ GEP_U87_strata + AGE + SEX + paper_MGMT.promoter.status + paper_IDH.status,
  data = gsva_patient
)
df_gep5 <- droplevels(model.frame(fit_cox_gep5))
p_forest_gep5 <- ggforest(fit_cox_gep5, data = df_gep5, fontsize = 0.9)


## ==== PFS =====

## --- GEP1 ---
fit_cox_gep1 <- coxph(
  Surv(DFS_MONTHS, DFS_STATUS) ~ GEP1_strata + AGE + SEX + paper_MGMT.promoter.status + paper_IDH.status,
  data = gsva_patient
)
df_gep1 <- droplevels(model.frame(fit_cox_gep1))
p_forest_gep1 <- ggforest(fit_cox_gep1, data = df_gep1, fontsize = 0.9)

## --- GEP2 ---
fit_cox_gep2 <- coxph(
  Surv(DFS_MONTHS, DFS_STATUS) ~ GEP2_strata + AGE + SEX + paper_MGMT.promoter.status + paper_IDH.status,
  data = gsva_patient
)
df_gep2 <- droplevels(model.frame(fit_cox_gep2))
p_forest_gep2 <- ggforest(fit_cox_gep2, data = df_gep2, fontsize = 0.9)

## --- GEP3 ---
fit_cox_gep3 <- coxph(
  Surv(DFS_MONTHS, DFS_STATUS) ~ GEP3_strata + AGE + SEX + paper_MGMT.promoter.status + paper_IDH.status,
  data = gsva_patient
)
df_gep3 <- droplevels(model.frame(fit_cox_gep3))
p_forest_gep3 <- ggforest(fit_cox_gep3, data = df_gep3, fontsize = 0.9)

## --- GEP4 ---
fit_cox_gep4 <- coxph(
  Surv(DFS_MONTHS, DFS_STATUS) ~ GEP4_strata + AGE + SEX + paper_MGMT.promoter.status + paper_IDH.status,
  data = gsva_patient
)
df_gep4 <- droplevels(model.frame(fit_cox_gep4))
p_forest_gep4 <- ggforest(fit_cox_gep4, data = df_gep4, fontsize = 0.9)

## --- GBM ---
colnames(gsva_patient)[colnames(gsva_patient)=="AGE"] <- "Age"
colnames(gsva_patient)[colnames(gsva_patient)=="SEX"] <- "Sex"
colnames(gsva_patient)[colnames(gsva_patient)=="paper_MGMT.promoter.status"] <- "MGMT"
colnames(gsva_patient)[colnames(gsva_patient)=="paper_IDH.status"] <- "IDH"

fit_cox_gep5 <- coxph(
  Surv(DFS_MONTHS, DFS_STATUS) ~ U87_uG + Age + Sex + MGMT + IDH,
  data = gsva_patient
)
df_gep5 <- droplevels(model.frame(fit_cox_gep5))
p_forest_gep5 <- ggforest(fit_cox_gep5, data = df_gep5, fontsize = 0.9)

ggsave(plot = p_forest_gep5,
       path="01_bulk-mrnaseq/results/figures/",
       filename="CoxPH_PFS_TCGA.svg",
       device = svglite::svglite,
       width=6, height=5, dpi=700)

library(tidyverse)
library(ggpubr)
library(rstatix)
library(emmeans)
source("utils/utils.R")

dt <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/spatial_stats.csv")
radial_layer_areas <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/radial_layer_areas.csv")

dt <- dt %>%
  mutate(Composition_lvl1 = if_else(grepl("U87", region), "U87", "U87_THP-1"),
         Composition_lvl2 = str_split(region, "_", Inf, T)[,3],
         Perturbation = if_else(grepl("uG", region), "uG", "KSC"))

## ===== GEPs Moran I =====

moran_dt <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/GEP_moranI.csv")

moran_dt <- moran_dt %>%
  mutate(Composition_lvl1 = if_else(grepl("U87", name), "U87", "U87_THP-1"),
         Composition_lvl2 = str_split(name, "_", Inf, T)[,2],
         Perturbation = if_else(grepl("uG", name), "uG", "KSC"))


fit <- lm(I ~ X + Composition_lvl2 + Perturbation, data = moran_dt)
summary(fit)
# PerturbationuG        0.08953    0.03536   2.532   0.0183 *
# overall, microgravity increases spatial heterogeneity

fit_int <- lm(I ~ X * Perturbation + Composition_lvl2, data = moran_dt)
summary(fit_int)

emm <- emmeans(fit_int, ~ Perturbation | X)
ug_effect_by_gep <- contrast(emm, method = "revpairwise") # revpairwise yields inverse effect (i.e. uG - control)
summary(ug_effect_by_gep, infer = TRUE, adjust = "holm")

stat.test <- summary(
  ug_effect_by_gep,
  infer = TRUE,
  adjust = "holm"
) |>
  as.data.frame() |>
  transmute(
    X = X,                     # grouping variable
    group1 = "KSC",           # control
    group2 = "uG",             # treatment
    p = p.value
  ) %>%
  left_join(
    moran_dt |>
      group_by(X) |>
      summarise(y.position = max(I, na.rm = TRUE) * 1.15),
    by = "X"
  ) %>%
  mutate(p = signif(p, digits = 3))

stat.test$p <- p.adjust(stat.test$p)
moran_dt$X <- gsub("_score", "", moran_dt$X)

p <- ggerrorplot(
  moran_dt,
  x = "Perturbation",
  y = "I",
  color = "Perturbation",
  desc_stat = "mean_se",
  add = "dotplot"
) +
  stat_pvalue_manual(stat.test, label = "p", hide.ns = T) +
  facet_wrap(. ~ X, nrow=1) +
  theme_cowplot(14) +
  scale_color_manual(values = COLPAL_UG_KSC) +
  labs(x = NULL, y = "Moran's I") +
  theme(
    strip.background = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(p, path="02_st-xenium/figures", filename= "gep_morans_I.pdf", width=11*1.1, height=6*1.1, dpi = 700, units = "cm")

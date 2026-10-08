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


## ===== Cells per layer =====

color_palette <- c(
  core      = "#8E0152",  # plum
  middle    = "#C51B7D",  # magenta-rose
  periphery = "#FDE0EF",  # pale blush
  outside   = "#BEBADA"   # desaturated lavender
)

plot_df <- dt %>%
  #filter(radial_layer_plot != "outside") %>%
  group_by(region, radial_layer_plot) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  mutate(Freq = Count / sum(Count)) %>%
  ungroup() %>%
  mutate(region = str_split(region, "-", Inf, T)[,2]) %>%
  mutate(Composition_lvl1 = if_else(grepl("U87", region), "U87", "U87_THP-1"),
         Composition_lvl2 = str_split(region, "_", Inf, T)[,2],
         Perturbation = if_else(grepl("uG", region), "uG", "KSC"))

plot_df <- plot_df %>%
  left_join(radial_layer_areas, by = c("region" = "sample", "radial_layer_plot" = "layer")) %>%
  filter(radial_layer_plot!="outside") %>%
  group_by(region) %>%
  mutate(total_cells = sum(Count, na.rm = TRUE),
            total_area  = sum(exact_area, na.rm = TRUE),
            expected    = total_cells * (exact_area / total_area),
            # expected to observed ratio
            obs_exp     = Count / expected)

library(emmeans)
fit <- lm(obs_exp ~ radial_layer_plot + region, data = plot_df)

emm <- emmeans(fit, ~ radial_layer_plot)          # returns model-adjusted means
emm_test <- test(emm, null = 1)                               # H0: mean = 1 for each layer
emm_test
# core               0.798 0.0458 14    1  -4.415  0.0006

pval <- emm_test$p.value[1]

y_pos <- plot_df %>%
  filter(radial_layer_plot == "core") %>%
  summarise(y = max(obs_exp, na.rm = TRUE)) %>%
  pull(y) %>% max()

p <- ggerrorplot(plot_df, x = "radial_layer_plot", y = "obs_exp",
                 color = "Perturbation",
                 desc_stat = "mean_se",
                 add = "dotplot") +
  labs(x = NULL, y = "Fold enrichment") +
  annotate(
    "text",
    x = 2,
    y = y_pos * 1.3,
    label = paste0("p = ", signif(pval, 2)),
    size = 5
  ) +
  theme_cowplot(14) +
  scale_color_manual(values = COLPAL_UG_KSC) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

ggsave(p, path="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures", filename= "cell_counts_per_layer.pdf", width=9*1.1, height=10*1.1, dpi = 700, units = "cm")


# compare ug vs ksc

fit_int <- lm(obs_exp ~ radial_layer_plot*Perturbation + Composition_lvl2, data = plot_df)
summary(fit_int)

emm <- emmeans(fit_int, ~ Perturbation | radial_layer_plot)
ug_effect_by_gep <- contrast(emm, method = "revpairwise") # revpairwise yields inverse effect (i.e. uG - control)
summary(ug_effect_by_gep, infer = TRUE, adjust = "holm")

stat.test <- summary(
  ug_effect_by_gep,
  infer = TRUE,
  adjust = "holm"
) |>
  as.data.frame() |>
  transmute(
    radial_layer_plot = radial_layer_plot,                     # grouping variable
    group1 = "KSC",           # control
    group2 = "uG",             # treatment
    p = p.value
  ) %>%
  left_join(
    plot_df |>
      group_by(radial_layer_plot) |>
      summarise(y.position = max(obs_exp, na.rm = TRUE) * 1.15),
    by = "radial_layer_plot"
  ) %>%
  mutate(p = signif(p, digits = 3))


p <- ggerrorplot(
  plot_df,
  x = "Perturbation",
  y = "obs_exp",
  color = "Perturbation",
  facet.by = "radial_layer_plot",
  desc_stat = "mean_se",
  add = "dotplot"
) +
  labs(x = NULL, y = "Fold enrichment") +
  stat_pvalue_manual(
    stat.test,
    label = "p",
    tip.length = 0.01,
    size = 4
  ) +
  theme_bw(14)

ggsave(p, path="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures", filename= "cell_counts_per_layer_ug_vs_ksc.pdf", width=15, height=12, dpi = 700, units = "cm")






## ===== GEPs per layer =====

plot_df <- dt %>%
  filter(radial_layer_plot != "outside") %>%
  group_by(region, Perturbation, Composition_lvl2, radial_layer_plot) %>%
  summarise(GEP1_mean = mean(GEP1_score),
            GEP2_mean = mean(GEP2_score),
            GEP3_mean = mean(GEP3_score),
            GEP4_mean = mean(GEP4_score)) %>%
  pivot_longer(starts_with("GEP"), names_to = "GEP", values_to = "mean_score")

p1 <- ggplot(plot_df, aes(
  x = Perturbation,
  y = mean_score,
  group = radial_layer_plot,
  color = radial_layer_plot
)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(Composition_lvl2 ~ GEP, scales = "free_y") +
  theme_bw(14) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = NULL,
    y = "GEP signature score",
    color = "Group"
  ) + 
  scale_color_brewer(palette = "Set2")
p1


## ====Binary GEP analysis======

library(introdataviz)

p <- ggplot(dt %>% pivot_longer(ends_with("_score"), names_to = "GEP", values_to = "GEP_score") %>% filter(radial_layer_plot != "outside"), aes(x = GEP, y = GEP_score, fill = Perturbation)) +
  introdataviz::geom_split_violin(alpha = .4, trim = FALSE) +
  geom_boxplot(width = .1, alpha = .6, show.legend = FALSE, outliers = FALSE) +
  facet_grid(radial_layer_plot ~ Composition_lvl2) +
  scale_fill_brewer(palette = "Set1", name = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = p, path="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures", filename= "GEP_binned_zonation.pdf", width=30, height=25, dpi = 700, units = "cm")



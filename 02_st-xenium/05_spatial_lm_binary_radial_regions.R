library(tidyverse)
library(ggpubr)
library(rstatix)
library(emmeans)
library(ggrepel)

dt <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/spatial_stats.csv")
radial_layer_areas <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/radial_layer_areas.csv")

dt <- dt %>%
  mutate(Composition_lvl1 = if_else(grepl("U87", region), "U87", "U87_THP-1"),
         Composition_lvl2 = str_split(region, "_", Inf, T)[,3],
         Perturbation = if_else(grepl("uG", region), "uG", "KSC"))

## ====Linear models: region-specific effects=======

lm_data <- dt %>% select(ends_with("_score"), radial_layer_plot, Composition_lvl2, Perturbation) %>% pivot_longer(ends_with("_score"), names_to = "GEP", values_to = "score") %>% filter(radial_layer_plot!="outside")

lm_data_list <- lm_data %>% group_by(Composition_lvl2, GEP) %>% group_split()
names(lm_data_list) <- lapply(lm_data_list, function(x) {paste0(unique(x[["Composition_lvl2"]]), "_", unique(x[["GEP"]]))})

lm_emmeans <- lapply(lm_data_list, function(x) {
  int_fit <- lm(score ~ radial_layer_plot * Perturbation, data = x)
  emm <- emmeans(int_fit, ~ Perturbation | radial_layer_plot)
  ug_effect_by_layer <- contrast(emm, method = "revpairwise")
  
  summary(ug_effect_by_layer, infer = TRUE, adjust = "holm") %>%
    as.data.frame() %>%
    transmute(
      radial_layer_plot = radial_layer_plot,
      group1 = "KSC",
      group2 = "uG",
      estimate = estimate,
      SE = SE,
      p = p.value)
})

lm_emmeans <- lm_emmeans %>%
  bind_rows(.id = "Group") %>%
  mutate(Composition_lvl2 = str_split(Group, "_", 2, T)[,1], GEP = str_split(Group, "_", 2, T)[,2], p_adj = p.adjust(p))

plot_df <- lm_emmeans %>%
  mutate(
    radial_layer_plot = factor(radial_layer_plot,
                               levels = c("periphery", "middle", "core"),
                               ordered = TRUE),
    sig = case_when(
      p_adj <= 1e-4 ~ "****",
      p_adj <= 1e-3 ~ "***",
      p_adj <= 1e-2 ~ "**",
      p_adj <= 5e-2 ~ "*",
      TRUE          ~ "ns"
    ),
    # only show labels when significant (per your earlier rule p_adj < 0.05)
    sig_lab = ifelse(p_adj < 0.05, sig, NA_character_)
  )

pd <- position_dodge(width = 0)

p <- ggplot(plot_df, aes(x = radial_layer_plot, y = estimate, color = GEP, group = GEP)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(position = pd, linewidth = 0.6, alpha = 0.5) +
  geom_point(position = pd, size = 2) +
  geom_errorbar(aes(ymin = estimate - SE, ymax = estimate + SE),
                position = pd, width = 0.5, linewidth = 0.5) +
  facet_wrap(~ Composition_lvl2) +
  labs(x = NULL, y = "Estimate (uG-KSC)", color = "GEP") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  geom_text_repel(
    data = subset(plot_df, p_adj < 0.05),
    aes(label = sig_lab),
    position = pd,
    size = 5,
    show.legend = FALSE,
    max.overlaps = Inf
  )

ggsave(plot = p, path="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures", filename= "GEP_region_effects_lm.pdf", width=20*.9, height=13*.9, dpi = 700, units = "cm")


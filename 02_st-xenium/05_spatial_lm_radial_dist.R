library(tidyverse)
library(ggpubr)
library(rstatix)
library(emmeans)
source("utils/utils.r")

# ------ Read data ------

dt <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/spatial_stats.csv")

dt <- dt %>%
  mutate(Composition_lvl1 = if_else(grepl("U87", region), "U87", "U87_THP-1"),
         Composition_lvl2 = str_split(region, "_", Inf, T)[,3],
         Perturbation = if_else(grepl("uG", region), "uG", "KSC"))

dt$radial_norm_boundary <- 1 - dt$radial_norm_boundary

# ----- Fit linear models ----

lm_data <- dt %>%
  filter(radial_layer_plot!="outside") %>%
  select(ends_with("_score"), radial_norm_boundary, Composition_lvl2, Perturbation) %>%
  pivot_longer(ends_with("_score"), names_to = "GEP", values_to = "score")

lm_data_list <- lm_data %>% group_by(Composition_lvl2, GEP) %>% group_split()
names(lm_data_list) <- lapply(lm_data_list, function(x) {paste0(unique(x[["Composition_lvl2"]]), "_", unique(x[["GEP"]]))})

lm_emmeans <- lapply(lm_data_list, function(x) {
  int_fit <- lm(score ~ radial_norm_boundary * Perturbation, data = x)
  slopes <- emtrends(int_fit, ~ Perturbation, var = "radial_norm_boundary")
  
  summary(slopes, infer = TRUE) %>%
    as.data.frame() %>%
    transmute(
      Perturbation = Perturbation,
      Effect = radial_norm_boundary.trend,
      SE = SE,
      p = p.value,
      p_diff = as.data.frame(pairs(slopes))$p.value)
})

lm_emmeans <- lm_emmeans %>%
  bind_rows(.id = "Group") %>%
  mutate(Composition_lvl2 = str_split(Group, "_", 2, T)[,1], GEP = str_split(Group, "_", 2, T)[,2], p_diff_adj = p.adjust(p_diff))

plot_df <- lm_emmeans %>%
  mutate(
    Perturbation = factor(Perturbation,
                          levels = c("KSC", "uG"),
                          ordered = TRUE),
    GEP = gsub("_score", "", GEP),
    sig = case_when(
      p_diff_adj <= 1e-4 ~ "****",
      p_diff_adj <= 1e-3 ~ "***",
      p_diff_adj <= 1e-2 ~ "**",
      p_diff_adj <= 5e-2 ~ "*",
      TRUE          ~ "ns"
    ),
    sig_lab = ifelse(p_diff_adj < 0.05, sig, NA_character_)
  )

ann_stat <- plot_df %>%
  group_by(GEP, Composition_lvl2) %>%   # keep facet var + x var
  summarise(
    y_top    = max(Effect + SE, na.rm = TRUE),
    y_bottom = min(Effect - SE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    group1 = "KSC",
    group2 = "uG",
    # if both bars are negative -> place bracket below
    y.position = ifelse(y_top < 0, y_bottom - 0.03, y_top + 0.03)
  ) %>%
  left_join(
    plot_df %>%
      select(GEP, Composition_lvl2, label = sig_lab) %>%
      distinct(),
    by = c("GEP", "Composition_lvl2")
  ) %>%
  add_x_position(x = "Composition_lvl2", group = "Perturbation", dodge = 0.8)

plot_df$Composition_lvl2 <- factor(plot_df$Composition_lvl2)
levels(plot_df$Composition_lvl2) <- c("GBM", "+M2", "+MΦ", "+Mono")
pd <- position_dodge(width = 0.8)
p <- ggplot(plot_df, aes(x = Composition_lvl2, y = Effect, fill = Perturbation)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_col(position = pd, width = 0.7, color = "black") +
  geom_errorbar(
    aes(ymin = Effect - SE, ymax = Effect + SE),
    position = pd, width = 0.2, linewidth = 0.5
  ) +
  scale_y_continuous(limits = c(-.55, .55)) +
  facet_wrap(~ GEP, nrow = 1) +
  labs(x = NULL, y = "Effect size", fill = NULL) +
  scale_fill_manual(values = COLPAL_UG_KSC) +
  custom_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_pvalue_manual(ann_stat, tip.length = 0)

#ggsave(plot = p, path="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures", device = cairo_pdf, filename= "GEP_radial_distance_lm.pdf", width=20*.9, height=6*1.2*.9, dpi = 700, units = "cm")

ggsave(
  plot = p,
  filename = "GEP_radial_distance_lm.svg",
  path = "/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures",
  device = svglite::svglite,
  width = 20*.9, height = 6*1.2*.9, units = "cm"
)

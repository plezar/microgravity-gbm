library(tidyverse)

dt <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/spatial_stats.csv")

dt <- dt %>%
  mutate(Composition_lvl1 = if_else(grepl("U87", region), "U87", "U87_THP-1"),
         Composition_lvl2 = str_split(region, "_", Inf, T)[,3],
         Perturbation = if_else(grepl("uG", region), "uG", "KSC"))


## ===== Cells per layer =====

plot_df <- dt %>%
  filter(radial_layer_plot != "outside") %>%
  group_by(region, radial_layer_plot) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  mutate(Freq = Count / sum(Count)) %>%
  ungroup()

p <- ggplot(plot_df, aes(x = region, y = Freq, fill = radial_layer_plot)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Region",
       y = "Fraction",
       fill = "Radial Layer") +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(p, path="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures", filename= "cell_counts_per_layer.pdf", width=30, height=15, dpi = 700, units = "cm")


## ===== GEPs Moran I =====

moran_dt <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/GEP_moranI.csv")

moran_dt <- moran_dt %>%
  mutate(Composition_lvl1 = if_else(grepl("U87", name), "U87", "U87_THP-1"),
         Composition_lvl2 = str_split(name, "_", Inf, T)[,2],
         Perturbation = if_else(grepl("uG", name), "uG", "KSC"))


moran_dt %>% ggboxplot(x = "X", y = "I", color = "Perturbation")
moran_dt %>% ggboxplot(x = "X", y = "I", color = "Composition_lvl2", add = "jitter")

stat.test <- moran_dt %>%
  group_by(X) %>%
  rstatix::wilcox_test(I ~ Perturbation) %>%
  add_significance("p") %>%
  add_xy_position(x = "X", fun = "mean_se")

p <- ggerrorplot(moran_dt, x = "X", y = "I",
                 color = "Perturbation",
                 desc_stat = "mean_se",
                 add = "jitter") +
  labs(x = NULL, y = "Moran's I") +
  stat_pvalue_manual(stat.test, label = "p") +
  theme_bw(14)


res <- lm(I ~ X + Perturbation + Composition_lvl2, data = moran_dt)
summary(res)
# PerturbationuG        0.08953    0.03536   2.532   0.0183 *
# overall, microgravity increases spatial heterogeneity

ggsave(p, path="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures", filename= "gep_morans_I.pdf", width=30, height=15, dpi = 700, units = "cm")



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


## ====GAM======

library(mgcv)
library(gratia)
library(ggplot2)
library(mgcv)
library(ggplot2)

dt$group <- interaction(dt$Perturbation, dt$Composition_lvl2, drop = TRUE)
dt$group <- as.factor(dt$group)
dt$Perturbation <- factor(dt$Perturbation)

fit_list <- list()
for (gep in paste0("GEP", 1:4, "_score")) {
  
  f <- as.formula(paste0(gep, " ~ s(radial_norm_boundary, by = group) + group"))
  fit <- gam(f,
             data = dt,
             family = gaussian())
  fit_list[[gep]] <- fit
}

## ====Partial effect smooths======

smooth_estimates_res <- lapply(fit_list, smooth_estimates)
smooth_estimates_res_df <- bind_rows(smooth_estimates_res, .id = "GEP")

p <- smooth_estimates_res_df %>%
  add_confint() %>%
  mutate(Composition_lvl1 = if_else(grepl("U87", group), "U87", "U87_THP-1"),
         Composition_lvl2 = str_split(group, "\\.", Inf, T)[,2],
         Perturbation = if_else(grepl("uG", group), "uG", "KSC")) %>%
  ggplot(aes(y = .estimate, x = radial_norm_boundary, color = Perturbation)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci),
              alpha = 0.2,
  ) +
  facet_grid(GEP ~ Composition_lvl2) +
  geom_line() +
  labs(
    y = "Partial effect",
    x = "Radial norm distance"
  ) +
  theme_bw(14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p, path="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures", filename= "smooth_curves.pdf", width=20, height=15, dpi = 700, units = "cm")


## ====Absolute effect smooths======

dt_nona <- dt %>% drop_na(radial_norm_boundary)

# limiting domain
min_plotted <- quantile(dt_nona$radial_norm_boundary, 0.01)
max_plotted <- quantile(dt_nona$radial_norm_boundary, 0.99)

# make a grid of radial positions for each condition
newdat <- expand.grid(
  radial_norm_boundary = seq(min_plotted, max_plotted, length.out = 200),
  group = unique(dt$group)
)

p_list <- lapply(fit_list, predict, newdata = newdat, type = "response", se.fit = TRUE)
p <- bind_rows(p_list, .id = "GEP")

newdat <- do.call(rbind, replicate(4, newdat, simplify = FALSE))
newdat$fit <- p$fit
newdat$GEP <- p$GEP
newdat$se  <- p$se.fit
newdat$lower <- newdat$fit - 1.96 * newdat$se
newdat$upper <- newdat$fit + 1.96 * newdat$se

newdat <- newdat %>%
  mutate(Composition_lvl1 = if_else(grepl("U87", group), "U87", "U87_THP-1"),
         Composition_lvl2 = str_split(group, "\\.", Inf, T)[,2],
         Perturbation = if_else(grepl("uG", group), "uG", "KSC"))

p <- ggplot(newdat, aes(x = radial_norm_boundary, y = fit, colour = Perturbation, fill = Perturbation)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  geom_line(size = 1) +
  geom_point(data = dt %>% pivot_longer(ends_with("_score"), names_to = "GEP", values_to = "GEP_score"), aes(x = radial_norm_boundary, y = GEP_score, color = Perturbation), alpha = 0.05, size = 0.5) +
  facet_grid(GEP ~ Composition_lvl2) +
  labs(
    x = "Normalized distance to boundary",
    y = "Fitted GEP score",
    title = "Radial GEP smooths for KSC vs uG"
  ) +
  theme_bw(base_size = 14)

ggsave(p, path="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures", filename= "absolute_smooth_curves.pdf", width=20, height=15, dpi = 700, units = "cm")

# UTD: lower GEP1 (Prolif/DDR), higher GEP3 (Inflam., in spheroid periphery)
# UTUD: slightly lower GEP3 (Inflam.) in the periphery
# UM2: higher GEP3 (Inflam., in spheroid periphery), low GEP2 (MES) in periphery but higher in the middle
# U87: lower GEP2 (MES) in the periphery of the spheroid, higher GEP3 (Inflam) in the periphery


## ====Overall effect of uG on spatial heterogeneity======

# parametric coefficients: intercepts relative to the reference level (U87_KSC)
# these are average GEP levels disregarding spatial gradients
# smooth terms: zonation analysis
# edf: effective degrees of freedom. measures spatial complexity and heterogeneity, not direction
# F/pval: test whether edf!=1

# test for the overall trend
fit2 <- gam(
  GEP1_score ~ 
    s(radial_norm_boundary, by = Perturbation) +
    Perturbation,
  data = dt,
  family = gaussian()
)

summary(fit2)
# edf Ref.df      F  p-value    
# s(radial_norm_boundary):PerturbationKSC 5.291  6.420  4.799 5.07e-05 ***
# s(radial_norm_boundary):PerturbationuG  6.225  7.394 38.510  < 2e-16 ***
# both conditions show significant radial zonation, but uG has a more complex smooth
# massive difference in F (which measures how strongly the radial smooth deviates from flatness)



## ====Does a model with separate smooths fit better than a model with a single smooth?======

f_stat_df <- data.frame()
for (gep in paste0("GEP", 1:4, "_score")) {
  fit_shared <- gam(
    as.formula(paste0(gep, " ~ s(radial_norm_boundary) + Perturbation")),
    data = dt,
    family = gaussian()
  )
  
  fit_sep <- gam(
    as.formula(paste0(gep, " ~ s(radial_norm_boundary, by = Perturbation) + Perturbation")),
    data = dt,
    family = gaussian()
  )
  
  anova_res <- anova(fit_shared, fit_sep, test = "F")
  
  f_stat_sub <- summary(fit_sep)$s.table %>% as.data.frame() %>% mutate(GEP = gep, anova_P = anova_res$`Pr(>F)`[2])
  f_stat_sub$Perturbation <- str_split(rownames(f_stat_sub), "Perturbation", Inf, T)[,2]
  f_stat_df <- rbind(f_stat_df, f_stat_sub)
}

mystat.test <- tibble(X = paste0("GEP", 1:4, "_score"),
           .y. = "F",
           group1 = "KSC",
           group2 = "uG",
           p = signif(unique(f_stat_df$anova_P), 3),
           y.position = f_stat_df %>% group_by(GEP) %>% summarise(max_val = max(F)+5) %>% pull(max_val),
           x = stat.test$x,
           xmin = stat.test$xmin,
           xmax = stat.test$xmax,
           groups = stat.test$groups)

p <- f_stat_df %>% ggbarplot("GEP", "F",
          fill = "Perturbation", color = "Perturbation",
          position = position_dodge(0.9)) +
  stat_pvalue_manual(mystat.test, label = "p")
ggsave(p, path="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures", filename= "F_stat_comparison.pdf", width=20, height=15, dpi = 700, units = "cm")

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

library(tidyverse)
library(ggpubr)
library(rstatix)
library(emmeans)
library(mgcv)
library(gratia)
library(ggplot2)
library(mgcv)
library(ggplot2)
source("utils/utils.R")

# ======= Prepare data ========

dt <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/spatial_stats.csv")
radial_layer_areas <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/radial_layer_areas.csv")

dt <- dt %>%
  mutate(Composition_lvl1 = if_else(grepl("U87", region), "U87", "U87_THP-1"),
         Composition_lvl2 = str_split(region, "_", Inf, T)[,3],
         Perturbation = if_else(grepl("uG", region), "uG", "KSC"))

## ====GAM======

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
         Perturbation = if_else(grepl("uG", group), "uG", "KSC"),
         GEP = gsub("_score", "", GEP)) %>%
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
  custom_theme +
  scale_color_manual(values = COLPAL_UG_KSC) +
  scale_x_continuous(breaks = c(0, 1), limits = c(0, 1))

ggsave(p, path="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures", filename= "smooth_curves.pdf", width=15, height=15, dpi = 700, units = "cm")


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
         Perturbation = if_else(grepl("uG", group), "uG", "KSC"),
         GEP = gsub("_score", "", GEP))

p <- ggplot(newdat, aes(x = radial_norm_boundary, y = fit, colour = Perturbation, fill = Perturbation)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  geom_line(size = 1) +
  geom_point(data = dt %>% pivot_longer(ends_with("_score"), names_to = "GEP", values_to = "GEP_score") %>% mutate(GEP = gsub("_score", "", GEP)), aes(x = radial_norm_boundary, y = GEP_score, color = Perturbation), alpha = 0.05, size = 0.01) +
  facet_grid(Composition_lvl2 ~ GEP, scales = "free_y") +
  labs(
    x = "Norm. radial distance",
    y = "Score (fitted)",
  ) +
  custom_theme +
  scale_color_manual(values = COLPAL_UG_KSC) +
  scale_fill_manual(values = COLPAL_UG_KSC) +
  scale_x_continuous(breaks = c(min_plotted, max_plotted), labels = c(1, 0), limits = c(min_plotted, max_plotted)) +
  scale_y_continuous(breaks = seq(-.25, 1, by = 0.5), limits = c(-.25, 1))

p
ggsave(p,
       path="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures",filename= "absolute_smooth_curves.pdf", width=22*.9, height=20.5*.9, dpi = 700, units = "cm")

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
# i.e. do we have evidence that uG and KSC need different shapes, or can they share the same slope?

f_stat_df <- data.frame()
edf_df <- data.frame()
for (group in unique(dt$Composition_lvl2)) {
  for (gep in paste0("GEP", 1:4, "_score")) {
    fit_shared <- gam(
      as.formula(paste0(gep, " ~ s(radial_norm_boundary) + Perturbation")),
      data = dt[dt$Composition_lvl2==group,],
      family = gaussian()
    )
    
    fit_sep <- gam(
      as.formula(paste0(gep, " ~ s(radial_norm_boundary, by = Perturbation) + Perturbation")),
      data = dt[dt$Composition_lvl2==group,],
      family = gaussian()
    )
    
    anova_res <- anova(fit_shared, fit_sep, test = "F")
    
    f_stat_sub <- data.frame(F = anova_res$F[2], P = anova_res$`Pr(>F)`[2], GEP = gep, Group = group)
    f_stat_df <- rbind(f_stat_df, f_stat_sub)
    
    edf_df <- 
  }
}

f_stat_df_plot <- f_stat_df %>%
  mutate(
    neglog10P = -log10(P),
    neglog10P_cap = pmax(neglog10P, -log10(0.05)),
    sig = P < 0.05
  )

ggplot(f_stat_df_plot,
       aes(x = GEP, y = Group,
           size = F,
           color = neglog10P_cap)) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(
    range = c(1.5, 8),
    name = "F statistic"
  ) +
  scale_color_gradient(
    low = "grey80",
    high = "#b2182b",    # deep red
    name = expression(-log[10](p)),
    breaks = c(-log10(0.05)),
    labels = c("0.05")
  ) +
  labs(
    x = "GEP",
    y = "Composition level 2"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey95"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



## ======= Does uG change the overall level of the GEP, averaged over radial position? ======

f_stat_df <- data.frame()
intercept_df <- data.frame()
for (group in unique(dt$Composition_lvl2)) {
  for (gep in paste0("GEP", 1:4, "_score")) {
    fit_no_base <- gam(as.formula(paste0(gep, " ~ s(radial_norm_boundary, by = Perturbation)")),
                       data = dt[dt$Composition_lvl2==group,],
                       family = gaussian()
    )
    
    fit_with_base <- gam(as.formula(paste0(gep, " ~ s(radial_norm_boundary, by = Perturbation) + Perturbation")),
                         data = dt[dt$Composition_lvl2==group,],
                         family = gaussian()
    )
    
    anova_res <- anova(fit_no_base, fit_with_base, test = "F")
    
    f_stat_sub <- data.frame(F = anova_res$F[2], P = anova_res$`Pr(>F)`[2], GEP = gep, Group = group)
    f_stat_df <- rbind(f_stat_df, f_stat_sub)
    
    ksc_intercept <- summary(fit_with_base)$p.table[,"Estimate"][1]
    effect_size <- summary(fit_with_base)$p.table[,"Estimate"][2]
    
    
    
    intercept_df <- rbind(intercept_df, data.frame(Condition = c("KSC", "UG"),
                                                   Intercept = c(ksc_intercept, ksc_intercept+effect_size),
                                                   SE = c(summary(fit_with_base)$p.table[,"Std. Error"][1], summary(fit_with_base)$p.table[,"Std. Error"][2]),
                                                   GEP = gep,
                                                   group = group))
  }
}

# 1. Ensure group is a factor
intercept_df$group <- as.factor(intercept_df$group)
group_levels <- levels(intercept_df$group)

# 2. Calculate max height using Bar + SE (Crucial Step!)
max_heights <- intercept_df %>%
  group_by(GEP, group) %>%
  # We take the top of the error bar (Intercept + SE) as the max height
  summarise(max_val = max(Intercept + SE, na.rm = TRUE), .groups = "drop")

# 3. Create the stats dataframe
stat.test <- f_stat_df %>%
  rename(group = Group, p = P) %>% 
  left_join(max_heights, by = c("GEP", "group")) %>%
  mutate(
    group1 = "KSC",
    group2 = "UG",
    p.signif = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ "ns"
    ),
    y.position = max_val * 1.1,  # Place bracket 10% above the error bar
    group_id = as.numeric(factor(group, levels = group_levels)),
    xmin = group_id - 0.2, 
    xmax = group_id + 0.2
  )

ggbarplot(intercept_df, 
          x = "group", 
          y = "Intercept", 
          facet.by = "GEP", 
          color = "Condition", 
          fill = "Condition",
          position = position_dodge(0.9)) +
  # --- NEW: Add Error Bars ---
  geom_errorbar(
    aes(ymin = Intercept - SE, ymax = Intercept + SE, group = Condition),
    position = position_dodge(0.9), 
    width = 0.2,       # Width of the whiskers
    color = "black"    # Keep error bars black for contrast
  ) +
  # ---------------------------
stat_pvalue_manual(
  stat.test, 
  label = "p.signif", 
  xmin = "xmin", 
  xmax = "xmax", 
  y.position = "y.position", 
  tip.length = 0.01,
  hide.ns = FALSE
)



## ====== Effect of spheroid conditions =======

dt$comp <- factor(dt$Composition_lvl2)

f_stat_df <- data.frame()
for (gep in paste0("GEP", 1:4, "_score")) {
  fit_H0 <- gam(
    as.formula(paste0(gep, " ~ comp + Perturbation + s(radial_norm_boundary) + s(radial_norm_boundary, by = Perturbation)")),
    data = dt,
    family = gaussian()
  )
  # null model where all spheroids share one baseline curve f(x)
  # uG has one extra smooth g(x) that's the same across all spheroid conditions
  # spheroids differ in intercepts only; the shapes of the smooths are the same
  
  fit_H1 <- gam(
    as.formula(paste0(gep, " ~ group + s(radial_norm_boundary, by = group)")),
    data = dt,
    family = gaussian()
  )
  # full interaction model where uG changes shape differently in each spheroid
  
  anova_res <- anova(fit_H0, fit_H1, test = "F")
  
  f_stat_sub <- data.frame(GEP = gep, F = anova_res$F[2], P = anova_res$`Pr(>F)`[2])
  f_stat_df <- rbind(f_stat_df, f_stat_sub)
}

f_stat_df_plot <- f_stat_df %>%
  mutate(
    neglog10P = -log10(P),
    neglog10P_cap = pmax(neglog10P, -log10(0.05)),
    sig = P < 0.05
  )

# ggplot displays Inf as grey. I am capping this value to min(P) in order to display it with a proper color
f_stat_df_plot$neglog10P_cap[f_stat_df_plot$neglog10P_cap==Inf] <- max(f_stat_df_plot$neglog10P_cap[f_stat_df_plot$neglog10P_cap!=Inf], na.rm = T)

ggplot(f_stat_df_plot,
       aes(x = GEP, y = "",
           size = F,
           color = neglog10P_cap)) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(
    range = c(1.5, 8),
    name = "F statistic"
  ) +
  scale_color_gradient(
    low = "grey80",
    high = "#b2182b",    # deep red
    name = expression(-log[10](p))
  ) +
  labs(
    x = "GEP",
    y = "Composition level 2"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey95"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

## ======Comparing wiggliness=======

### =====Using curvature energy (integrated squared curvature)====

gam_list <- lapply(paste0("GEP", 1:4, "_score"), function (x) {
  gam(
    as.formula(paste0(x, " ~ group + s(radial_norm_boundary, by = group)")),
    data = dt,
    family = gaussian()
  )
})

names(gam_list) <- paste0("GEP", 1:4, "_score")

d2_list <- lapply(gam_list, derivatives, order = 2)


d2 <- bind_rows(d2_list, .id = "GEP")

ce_data <- d2 %>%
  group_by(GEP, group) %>%
  summarise(CE = mean(.derivative^2)) %>%
  left_join(dt %>% distinct(Perturbation, group, comp), "group")


p1 <- ggplot(ce_data,
             aes(x = comp,
                 y = CE,
                 fill = Perturbation)) +
  geom_col(
    position = position_dodge(width = 0.8),
    color = "black",
    width = 0.7
  ) +
  facet_wrap(~ GEP, scales = "free_y") +
  scale_fill_brewer(
    palette = "Paired",
    name = "Perturbation"
  ) +
  labs(
    x = NULL,
    y = "Curvature Energy"
  ) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

### =====Using EDFs====

# verify that the chosen basis is the same
# lapply(gam_list, gam.check)

edf_data <- bind_rows(lapply(gam_list, function(x) {y <- as.data.frame(summary(x)$s.table); y$smooth <- rownames(y); return(y)}), .id = "GEP")

edf_data$group <- str_split(edf_data$smooth, "group", Inf, T)[,2]

edf_data <- edf_data %>% left_join(dt %>% distinct(Perturbation, group, comp), "group")

p2 <- ggplot(edf_data,
             aes(x = comp,
                 y = edf,
                 fill = Perturbation)) +
  geom_col(
    position = position_dodge(width = 0.8),
    color = "black",
    width = 0.7
  ) +
  facet_wrap(~ GEP, scales = "free_y") +
  scale_fill_brewer(
    palette = "Paired",
    name = "Perturbation"
  ) +
  labs(
    x = NULL,
    y = "EDF"
  ) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

p1 + p2


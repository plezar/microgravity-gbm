library(mgcv)
library(gratia)
library(tidyverse)
library(mgcv)
library(emmeans)
source("utils/utils.R")

# Xenium sc expression
norm_d <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/data/filtered_norm_xenium_counts.csv", row.names = 1)
norm_d <- t(norm_d)
d_meta <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/data/filtered_norm_xenium_counts_obs_metadata.csv", row.names = 1)

# GEP assignments
gep_assignments <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/xenium_cNMF_k4_gep_assignments.csv", row.names = 1)

# GEP gene loadings
gep_scores <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/xenium_cNMF_k4_gep_scores.csv", row.names = 1)

# DA proteins
de_proteins <- readLines("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/03_olink/results/olink_da_proteins.txt")

# select proteins that are part of the panel
de_proteins <- de_proteins[de_proteins %in% rownames(norm_d)]

# Subset expression data
gex_d <- t(norm_d[de_proteins,])
colnames(gex_d) <- paste0(colnames(gex_d), "_GEX")

# Add metadata
d <- cbind(gex_d, gep_assignments)
d <- d %>% mutate(GEP1 = as.numeric(grepl("1", GEP_assignment)),
             GEP2 = as.numeric(grepl("2", GEP_assignment)),
             GEP3 = as.numeric(grepl("3", GEP_assignment)),
             GEP4 = as.numeric(grepl("4", GEP_assignment)))

d_long <- d %>% pivot_longer(ends_with("_GEX"), names_to = "gene_name", values_to = "E")

## ===== Heatmap of GEP scores ========
gep_scores_sub <- gep_scores[de_proteins,]
colnames(gep_scores_sub) <- paste0("GEP", 1:4)

pdf("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures/GEP_weights_olink_hits.pdf", width = 3, height = 3)
Heatmap(gep_scores_sub,
        show_row_dend = F,
        show_column_dend = F,
        border = T,
        name = "W")
dev.off()

## ======= Average Expression in GEP-assigned cells =======
gep_E <- lapply(paste0("GEP", 1:4), function (x) {
  d_long %>%
    group_by(gene_name, batch, .data[[x]]) %>%
    summarise(Mean_E = mean(E)) %>%
    filter(.data[[x]] == 1) %>%
    dplyr::rename(GEP = .data[[x]]) %>%
    mutate(GEP = x)
})

gep_E <- bind_rows(gep_E)

gep_E$Perturbation <- str_split(str_split(gep_E$batch, "_", Inf, T)[,1], "-", Inf, T)[,2]

p <- ggerrorplot(gep_E, x = "GEP", y = "Mean_E",
            color = "GEP",
            facet.by = "gene_name",
            desc_stat = "mean_se",
            add = "dotplot") +
  ylab("Norm. counts")
ggsave(p, path="02_st-xenium/figures", filename= "olink_hits_exp_by_GEP.pdf", width=10, height=6, units = "in", dpi = 700)

## ======= Enrichment Analysis =======

d <- d %>%
  mutate(
    across(
      ends_with("_GEX"),
      ~ if_else(.x > 0, 1L, 0L),
      .names = "{.col}_bin"
    )
  )

fisher.test(table(d$IL33_GEX_bin, d$GEP3))
fisher.test(table(d$CSF2_GEX_bin, d$GEP2))
#fisher.test(table(d$IL33_GEX_bin, d$GEP4))

table(d$IL33_GEX_bin, d$GEP3, d$batch)


## ====== Spatial Analysis =======

spatial_dt <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/spatial_stats.csv")
d <- cbind(gex_d, spatial_dt)
d_long <- d %>% pivot_longer(ends_with("_GEX"), names_to = "gene_name", values_to = "E") %>% filter(radial_layer_plot != "outside")
d_long$radial_layer_plot <- factor(d_long$radial_layer_plot, levels = c("periphery", "middle", "core"))
d_long$Perturbation <- str_split(str_split(d_long$region, "-", Inf, T)[,2], "_", Inf, T)[,1]
d_long$Composition_lvl2 <- str_split(str_split(d_long$region, "-", Inf, T)[,2], "_", Inf, T)[,2]

d_long_summarised <- d_long %>% group_by(region, gene_name, radial_layer_plot, Perturbation) %>%
  summarise(Mean_E = mean(E))

p <- ggerrorplot(d_long_summarised, x = "radial_layer_plot", y = "Mean_E",
            color = "Perturbation",
            facet.by = "gene_name",
            desc_stat = "mean_se",
            add = "dotplot") +
  ylab("Norm. counts")
ggsave(p, path="02_st-xenium/figures", filename= "olink_hits_exp_by_radial_layer.pdf", width=10, height=6, units = "in", dpi = 700)

## ======= IL33 ========
il33_fit <- lm(E ~ Perturbation * radial_layer_plot,
               data = d_long[d_long$gene_name == "IL33_GEX", ])

emm <- emmeans(il33_fit, ~ Perturbation | radial_layer_plot)
ug_effect_by_radial_layer <- contrast(emm, method = "revpairwise")

# tidy for ggplot
plot_df <- as.data.frame(summary(ug_effect_by_radial_layer, infer = TRUE, adjust = "holm")) %>%
  mutate(
    radial_layer_plot = factor(radial_layer_plot, levels = unique(radial_layer_plot)),
    contrast = gsub(" - ", " vs ", contrast)
  )

# y position for p-value labels
pd <- position_dodge(width = 0.6)

plot_df$radial_layer_plot <- factor(plot_df$radial_layer_plot, levels = c("core", "middle", "periphery"))
p <- ggplot(plot_df, aes(x = radial_layer_plot, y = estimate, color = contrast)) +
  geom_hline(yintercept = 0, linewidth = 0.5, lty=2) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0, linewidth = 0.8,
                position = pd) +
  geom_point(size = 4, position = pd) +
  geom_text(
    aes(y = y_pos * 1.05, label = if_else(p.value < 0.1, paste0("p=", signif(p.value, 2)), "")),
    position = position_dodge2(width = 0.6, preserve = "single"),
    size = 5, show.legend = FALSE
  ) +
  labs(x = NULL, y = "Effect size") +
  theme_cowplot(14) +
  scale_color_manual(values = c("black")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "none"
  )

p
ggsave(p, path="02_st-xenium/figures", filename= "il33_lm_emmeans.pdf", width=2.75, height=4, units = "in", dpi = 700)

## ======= CSF2 ========
csf2_fit <- lm(E ~ Perturbation * radial_layer_plot,
               data = d_long[d_long$gene_name == "CSF2_GEX", ])

emm <- emmeans(csf2_fit, ~ Perturbation | radial_layer_plot)
ug_effect_by_radial_layer <- contrast(emm, method = "revpairwise")

# tidy for ggplot
plot_df <- as.data.frame(summary(ug_effect_by_radial_layer, infer = TRUE, adjust = "holm")) %>%
  mutate(
    radial_layer_plot = factor(radial_layer_plot, levels = unique(radial_layer_plot)),
    contrast = gsub(" - ", " vs ", contrast)
  )

# y position for p-value labels
pd <- position_dodge(width = 0.6)

plot_df$radial_layer_plot <- factor(plot_df$radial_layer_plot, levels = c("core", "middle", "periphery"))
p <- ggplot(plot_df, aes(x = radial_layer_plot, y = estimate, color = contrast)) +
  geom_hline(yintercept = 0, linewidth = 0.5, lty=2) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0, linewidth = 0.8,
                position = pd) +
  geom_point(size = 4, position = pd) +
  geom_text(
    aes(y = if_else(radial_layer_plot=="periphery", y_pos * 1.05, y_pos * 0.9),
        label = case_when(
          p.value < 0.001 ~ "p<0.001",
          p.value < 0.1   ~ paste0("p=", signif(p.value, 2)),
          TRUE            ~ ""
        )),
    position = position_dodge2(width = 0.6, preserve = "single"),
    size = 5, show.legend = FALSE
  ) +
  labs(x = NULL, y = "Effect size") +
  theme_cowplot(14) +
  scale_color_manual(values = c("black")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "none"
  )

p
ggsave(p, path="02_st-xenium/figures", filename= "csf2_lm_emmeans.pdf", width=2.75, height=4, units = "in", dpi = 700)



## ====== GAMs =========

ftrs <- c("CSF2_GEX", "IL33_GEX")
d_long$group <- interaction(d_long$Perturbation, d_long$Composition_lvl2, drop = TRUE)
d_long$group <- as.factor(d_long$group)
d_long$Perturbation <- factor(d_long$Perturbation)

fit_list <- list()
for (gene in ftrs) {
  
  f <- as.formula(paste0("E ~ s(radial_norm_boundary, by = group) + group"))
  fit <- gam(f,
             data = d_long[d_long$gene_name==gene,],
             family = gaussian())
  fit_list[[gene]] <- fit
}

## ====Absolute effect smooths======

d_long_nona <- d_long %>% drop_na(radial_norm_boundary)

# limiting domain
min_plotted <- quantile(d_long_nona$radial_norm_boundary, 0.01)
max_plotted <- quantile(d_long_nona$radial_norm_boundary, 0.99)

# make a grid of radial positions for each condition
newdat <- expand.grid(
  radial_norm_boundary = seq(min_plotted, max_plotted, length.out = 200),
  group = unique(d_long_nona$group)
)

p_list <- lapply(fit_list, predict, newdata = newdat, type = "response", se.fit = TRUE)
p <- bind_rows(p_list, .id = "gene_name")

newdat <- do.call(rbind, replicate(4, newdat, simplify = FALSE))
newdat$fit <- p$fit
newdat$gene_name <- p$gene_name
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
  geom_point(data = d_long[d_long$gene_name%in%ftrs,], aes(x = radial_norm_boundary, y = E, color = Perturbation), alpha = 0.5, size = 0.5) +
  facet_grid(gene_name ~ Composition_lvl2) +
  labs(
    x = "Normalized distance to boundary",
    y = "Fitted GEP score",
    title = "Radial GEP smooths for KSC vs uG"
  ) +
  theme_bw(base_size = 14)

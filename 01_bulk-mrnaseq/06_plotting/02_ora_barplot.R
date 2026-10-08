library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)
library(scales)
library(patchwork)
library(cowplot)
source("01_bulk-mrnaseq/util.R")


# === U87 ====

up_d <- readRDS("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/results/enrichment/deseq2_uG_vs_KSC_U87/ORA/GO_ORA_UP_ALL.rds")
down_d <- readRDS("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/results/enrichment/deseq2_uG_vs_KSC_U87/ORA/GO_ORA_DOWN_ALL.rds")

down_d_s <- simplify(down_d)
up_d_s <- simplify(up_d)

#down_d_s@result %>% View()
down_go_terms <- c(
  "GO:0040036", "GO:2000648", "GO:0043568", "GO:0001837",
  "GO:0001570", "GO:0010906", "GO:0019318", "GO:0002053",
  "GO:0099170", "GO:0043410", "GO:0043122", "GO:0072089"
)

# data
down_df <- prep_enrich(down_d_s@result %>% filter(ID %in% down_go_terms), n = 10, ontology = "BP")
up_df   <- prep_enrich(up_d_s@result,   n = 10, ontology = c("BP", "CC", "MF"))

# Low/High around base_blue
blue_low  <- "#A0C4D0"  # lighter
blue_high <- "#1D5162"  # darker

# Low/High for black
black_low  <- "#BFBFBF" # light gray
black_high <- "#000000" # black

# Handy vectors for ggplot2 gradients
blue_grad  <- c(low = blue_low,  high = blue_high)
black_grad <- c(low = black_low, high = black_high)

# two separate plots (each legend on the right)
p_down <- make_panel(
  down_df, direction = "down",
  fill_low = "#d6e6f2", fill_high = "#0b3c6f",
  #fill_low = black_low, fill_high = black_high,
  legend_pos = "bottom",
  axis_shrink = 3.5,
  show_y_right = FALSE
) + xlab(NULL)

p_up <- make_panel(
  up_df, direction = "up",
  fill_low = "#f6d8cb", fill_high = "#7b0f0f",
  #fill_low = blue_low, fill_high = blue_high,
  legend_pos = "bottom",
  axis_shrink = 2.9,
  show_y_right = FALSE
) + xlab(NULL)

# stack in 1 column
p <- p_down + p_up + plot_layout(widths = c(1, 1))
p

ggsave(p, path="01_bulk-mrnaseq/results/figures", filename= "U87_ORA_barplot.pdf", width=14*1.1, height=7*1.1, dpi = 700, units = "cm")


# === U87 + THP ====

up_d <- readRDS("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/results/enrichment/deseq2_uG_vs_KSC_U87THP1/ORA/GO_ORA_UP_ALL.rds")
down_d <- readRDS("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/results/enrichment/deseq2_uG_vs_KSC_U87THP1/ORA/GO_ORA_DOWN_ALL.rds")

down_d_s <- simplify(down_d)
up_d_s <- simplify(up_d)

#up_d_s@result %>% View()
up_go_terms <- c("GO:0006935", "GO:1990266", "GO:0002224", "GO:0002544", "GO:0042104", "GO:0035296", "GO:2001242", "GO:0050798", "GO:0048013", "GO:0032411")


# data
down_df <- prep_enrich(down_d_s@result, n = 10, ontology = c("CC", "MF"))
up_df   <- prep_enrich(up_d_s@result %>% filter(ID %in% up_go_terms),   n = 10, ontology = "BP")

# Low/High around base_blue
#blue_low  <- "#A0C4D0"  # lighter
#blue_high <- "#1D5162"  # darker

# Low/High for black
#black_low  <- "#BFBFBF" # light gray
#black_high <- "#000000" # black

# Handy vectors for ggplot2 gradients
blue_grad  <- c(low = blue_low,  high = blue_high)
black_grad <- c(low = black_low, high = black_high)

# two separate plots (each legend on the right)
p_down <- make_panel(
  down_df, direction = "down",
  fill_low = "#d6e6f2", fill_high = "#0b3c6f",
  #fill_low = black_low, fill_high = black_high,
  legend_pos = "bottom",
  axis_shrink = 3.5,
  show_y_right = FALSE
) +
  xlab(NULL)

p_up <- make_panel(
  up_df, direction = "up",
  fill_low = "#f6d8cb", fill_high = "#7b0f0f",
  #fill_low = blue_low, fill_high = blue_high,
  legend_pos = "bottom",
  axis_shrink = 3.5,
  show_y_right = FALSE
) + xlab(NULL)

# stack in 1 column
p <- p_down + p_up + plot_layout(widths = c(1, 1))
p

ggsave(p, path="01_bulk-mrnaseq/results/figures", filename= "U87THP1_ORA_barplot.pdf", width=15.4, height=7.7, dpi = 700, units = "cm")


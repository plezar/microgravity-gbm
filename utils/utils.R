library(cowplot)
library(ggplot2)

custom_theme <- theme_cowplot(14) +
  theme(
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank()
  )

# ======== color palettes =========

COLPAL_UG_KSC <- c("KSC" = "lightgrey", "uG" = "black")
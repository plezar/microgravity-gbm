library(spacexr)
library(tidyverse)
library(cowplot)

RCTD <- readRDS("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/data/RCTD_obj.rds")
weights <- as.data.frame(as.matrix(RCTD@results$weights))

p <- ggplot(weights, aes(x = THP1, y = U87)) +
  geom_point(alpha = 0.2, size = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "THP-1 weight",
    y = "U87 weight"
  ) +
  xlim(0, 1.45) +
  ylim(0, 1.45) +
  theme_cowplot(14) +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

ggsave(p, path="02_st-xenium/figures", filename= "RCTD_weights.pdf", width=10, height=8.4, units="cm", dpi = 700)
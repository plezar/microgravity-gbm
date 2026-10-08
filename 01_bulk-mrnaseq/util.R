prep_enrich <- function(res, n = 10, ontology = "BP", clip_width = 38) {
  # ontology can be a vector
  res %>%
    as.data.frame() %>%
    dplyr::filter(ONTOLOGY %in% c(ontology)) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = n) %>%
    dplyr::mutate(
      Description = stringr::str_trunc(Description, width = clip_width, side = "right", ellipsis = "..."),
      FDR = p.adjust,
      neglogFDR = -log10(FDR)
    )
}

make_panel <- function(df,
                       direction = c("down", "up"),
                       fill_low, fill_high,
                       legend_pos = "right",
                       x_breaks = NULL,
                       axis_shrink = 3,
                       show_y_right = FALSE) {
  direction <- match.arg(direction)
  
  df <- df %>%
    mutate(
      FE_signed   = if (direction == "down") -FoldEnrichment else FoldEnrichment,
      Description = fct_reorder(Description, FE_signed)
    )
  
  # small offset so the label sits just past the bar tip
  max_fe <- max(abs(df$FE_signed), na.rm = TRUE)
  x_off <- 0.02 * max_fe
  
  p <- ggplot(df, aes(x = FE_signed, y = Description, fill = neglogFDR)) +
    geom_col(width = 0.8, color = "black") +
    # put the "y-axis labels" as text at the top edge of each bar (slightly above)
    geom_text(
      data = df,
      aes(
        x = FE_signed + ifelse(direction == "down", -x_off, x_off),
        y = Description,
        label = Description
      ),
      inherit.aes = FALSE,
      size = 3,
      nudge_y = -.1,                         # "top of the bar"
      hjust = if (direction == "down") 1 else 0,
      vjust = 0
    ) +
    scale_fill_gradient(
      low = fill_low, high = fill_high,
      labels = label_number(accuracy = 0.1),
      name = expression(-log[10](P))
    ) +
    theme_cowplot(14) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line = element_blank(),
      legend.position = legend_pos,
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9),
      legend.key.height = unit(7.5, "pt"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
      #plot.margin = margin(5.5, 30, 5.5, 5.5)  # room for labels near panel edge
    ) +
    #guides(fill = guide_colorbar(barwidth = unit(10, "pt"))) +
    labs(x = "Fold Enrichment") +
    coord_cartesian(clip = "off")
  
  if (is.null(x_breaks)) {
    max_fe <- ceiling(max(df$FoldEnrichment, na.rm = TRUE))
    
    if (max_fe < 100) {
      x_breaks <- seq(0, max_fe, by = max(1, round(max_fe / 3)))
    } else {
      x_breaks <- seq(0, max_fe, by = max(1, round(max_fe / 2)))
    }
  }
  
  # "Down" extends left, but labels are positive
  if (direction == "down") {
    p <- p +
      scale_x_continuous(
        breaks = -x_breaks,
        labels = x_breaks,
        expand = expansion(mult = c(0.02, 0.02)),
        limits = c(-axis_shrink*max_fe, 0)
      )
  } else {
    p <- p +
      scale_x_continuous(
        breaks = x_breaks,
        labels = x_breaks,
        expand = expansion(mult = c(0.02, 0.02)),
        limits = c(0, axis_shrink*max_fe)
      )
  }
  
  # keep this behavior if you still want; axis labels are hidden anyway
  if (show_y_right) {
    p <- p + scale_y_discrete(position = "right")
  }
  
  p
}
library(tidyverse)
library(ggrepel)
library(cowplot)
library(patchwork)

# ==== Parse de novo hits ====

extract_nonredundant_headers <- function(motif_file) {
  x <- readLines(motif_file, warn = FALSE)
  hdr <- x[grepl("^>", x)]
  
  parse_one <- function(line) {
    line2 <- sub("^>", "", line)
    fields <- strsplit(line2, "\t", fixed = TRUE)[[1]]
    
    # fields[2] contains "...BestGuess:XXX/YYY..."
    bestguess_full <- sub(".*BestGuess:([^\\t]+).*", "\\1", fields[2])
    
    # TF is between BestGuess: and first "/"
    tf <- sub(".*BestGuess:([^/]+)/.*", "\\1", fields[2])
    
    # HOMER logP column (in your example it's field 4)
    logP <- suppressWarnings(as.numeric(fields[4]))
    
    # P-value lives in the last field, like "...P:1e-10"
    p_str <- sub(".*\\bP:([^\\s]+).*", "\\1", fields[length(fields)])
    p_val <- suppressWarnings(as.numeric(p_str))
    
    data.frame(
      TF = tf,
      BestGuess_full = bestguess_full,
      logP = logP,
      P = p_val,
      header = line,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, lapply(hdr, parse_one))
}

df <- extract_nonredundant_headers("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/results/homer_out/U87_DE_DOWN/nonRedundant.motifs")
df <- df %>% arrange(logP)
df$Rank <- 1:nrow(df)

ggplot(df, aes(x=Rank, y=logP)) +
  geom_point(size=0.5, shape=21) +
  geom_text_repel(data = df %>% filter(TF %in% c("ELF3")),
                  aes(x = Rank, y = logP, label = TF), color="black",
                  #min.segment.length = 1000,
                  seed = 24,
                  box.padding = .5,
                  #nudge_x = 100,
                  #nudge_y = 2,
                  #force = 300,
                  min.segment.length = 0,
                  show.legend = FALSE,
                  max.overlaps =15,
                  size=5
  ) +
  theme_cowplot(10) +
  labs(x="Rank",
       y="logP") + 
  theme(legend.position = "none"
  )  + scale_y_reverse() + scale_x_continuous(breaks=c(0, 17))


# ==== Parse known hits ====

## ==== U87 =====

known_hits <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/results/homer_out/U87_DE_DOWN/knownResults.txt", sep="\t")
known_hits$Rank <- 1:nrow(known_hits)
known_hits$Motif.Name <- str_split(known_hits$Motif.Name, "/", Inf, T)[,1]
known_hits$Color_Group <- if_else(known_hits$Log.P.value < log(0.01), "Yes", "No")
label_targets <- c("TEAD1(TEAD)", "TEAD2(TEA)")

p1 <- ggplot(known_hits, aes(x=Rank, y=Log.P.value, color = Color_Group)) +
  annotate(
    "rect",
    xmin =  -Inf, xmax =  Inf,
    ymin = log(0.01), ymax = -Inf,
    fill  = "#F6E27F",
    alpha = 0.30,
    colour = NA
  ) +
  geom_point(size=0.5, shape=21) +
  geom_text_repel(data = known_hits %>% filter(Motif.Name %in% label_targets),
                  aes(x = Rank, y = Log.P.value, label = Motif.Name), color="black",
                  #min.segment.length = 1000,
                  seed = 24,
                  box.padding = .5,
                  #nudge_x = 0.1,
                  #nudge_y = -1,
                  #force = 500,
                  min.segment.length = 0,
                  show.legend = FALSE,
                  max.overlaps =15,
                  size=5
  ) +
  scale_color_manual(values = c("Yes" = "black", "No" = "grey85")) +
  theme_cowplot(14) +
  labs(x="Rank",
       y="logP") + 
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.title = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
        #panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        ) +
  scale_y_reverse(limits = c(0, -10)) +
  scale_x_continuous(breaks=c(0, 400),  limits = c(-50, 500)) +
  ggtitle("GBM")

p1

## ==== U87 THP ====

known_hits <- read.csv("01_bulk-mrnaseq/results/homer_out/U87THP1_DE_UP/knownResults.txt", sep="\t")
known_hits$Rank <- 1:nrow(known_hits)
known_hits$Motif.Name <- str_split(known_hits$Motif.Name, "/", Inf, T)[,1]
known_hits$Color_Group <- if_else(known_hits$Log.P.value < log(0.01), "Yes", "No")

label_targets <- c("Smad3(MAD)")
p2 <- ggplot(known_hits, aes(x=Rank, y=Log.P.value, color = Color_Group)) +
  annotate(
    "rect",
    xmin =  -Inf, xmax =  Inf,
    ymin = log(0.01), ymax = -Inf,
    fill  = "#F6E27F",
    alpha = 0.30,
    colour = NA
  ) +
  geom_point(size=0.5, shape=21) +
  geom_text_repel(data = known_hits %>% filter(Motif.Name %in% label_targets),
                  aes(x = Rank, y = Log.P.value, label = Motif.Name), color="black",
                  #min.segment.length = 1000,
                  seed = 24,
                  box.padding = .5,
                  #nudge_x = 0.1,
                  #nudge_y = 2,
                  #force = 500,
                  min.segment.length = 0,
                  show.legend = FALSE,
                  max.overlaps =15,
                  size=5
  ) +
  scale_color_manual(values = c("Yes" = "black", "No" = "grey85")) +
  theme_cowplot(14) +
  labs(x="Rank",
       y=NULL) + 
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
        #panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        ) +
  scale_y_reverse(limits = c(0, -10)) +
  scale_x_continuous(breaks=c(0, 400), limits = c(-50, 500)) +
  ggtitle("GBM+Mono")

p <- p1 + plot_spacer() + p2 + plot_layout(widths = c(1, 0, 1))
ggsave(plot = p, path="01_bulk-mrnaseq/results/figures", filename= "HOMER_top_known_hits.pdf", width=16, height=7, dpi = 700, units = "cm")


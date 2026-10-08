library(patchwork)
source("utils/volcano_plot.R")

# === GBM + Mono ====

d <- read.csv("01_bulk-mrnaseq/results/deseq2_uG_vs_KSC_U87THP1.csv", row.names = 1)
d$GeneSymbol <- rownames(d)

summary(-log10(d$padj))
summary(d$log2FoldChange)

p1 <- volcano_plot(d,
             pval_colname = "pvalue",
             padj_colname = "padj",
             logfc_colname = "log2FoldChange",
             genesymbol_colname = "GeneSymbol",
             title = "GBM + Mono",
             label_top_n = 20,
             fc_threshold = 0,
             ylim_min = 0,
             ylim_max = 5,
             xlim_min = -6.5,
             xlim_max = 6.5,
             ggrepel_force = 50,
             ggrepel_max_overlaps = 5,
             col_vec = setNames(c("#0b3c6f", "grey85", "#7b0f0f"), 
                                c(-1, 0, 1)),
             rect_col_vec = setNames(c("#d6e6f2", "#f6d8cb"), 
                                     c(-1, 1)),
             padj_threshold = 0.05) 

ggsave(p1, path="01_bulk-mrnaseq/results/figures", filename= "U87THP1_volcanic_plot.pdf", width=10, height=10, dpi = 700, units = "cm")


# === GBM alone ====

d <- read.csv("01_bulk-mrnaseq/results/deseq2_uG_vs_KSC_U87.csv", row.names = 1)
d$GeneSymbol <- rownames(d)

summary(-log10(d$padj))
summary(d$log2FoldChange)

theme_genes <- c(
  # Stress / hypoxia / ISR
  "ATF4",      # integrated stress response (ISR)
  "P4HA1",     # collagen proline hydroxylation; ECM remodeling; often hypoxia/HIF-linked
  "ANGPTL4",   # secreted stress/hypoxia-responsive factor; vascular/lipid-handling contexts
  "CYBRD1",    # redox/iron-handling-associated stress contexts
  
  # Adhesion / cytoskeleton / migration
  "ITGAV",     # integrin alpha V; ECM adhesion / migration
  "FLNB",      # actin-binding cytoskeletal scaffold; motility/mechanics
  "TPM4",      # actin filament regulation; contractility/motility
  "EPB41L2",   # cytoskeletal organization / membrane-cortex linkage
  "SHTN1",     # cytoskeletal dynamics; neurite/migration-associated in some contexts
  "GJA1",      # Connexin43; gap junctions; cell-cell communication; migration-associated
  
  # RNA processing / translation / surveillance
  "UPF3A",     # nonsense-mediated decay (NMD)
  "UPF3B",     # nonsense-mediated decay (NMD)
  "HNRNPU",    # RNA binding; splicing/processing
  "ZRSR2",     # spliceosome component; U12-type intron splicing
  "EIF5B",     # translation initiation factor
  "EIF3J",      # translation initiation factor
  
  # Inflammatory output
  "CXCL3"    # chemokine; inflammatory/activated programs (uncomment if you want it included)
)

p2 <- volcano_plot(d,
                   pval_colname = "pvalue",
                   padj_colname = "padj",
                   logfc_colname = "log2FoldChange",
                   genesymbol_colname = "GeneSymbol",
                   title = "GBM",
                   #ylab = NULL,
                   label_top_n = 20,
                   fc_threshold = 0,
                   ylim_min = 0,
                   ylim_max = 5,
                   xlim_min = -6.5,
                   xlim_max = 6.5,
                   ggrepel_force = 50,
                   ggrepel_max_overlaps = 5,
                   genes_to_label = theme_genes,
                   col_vec = setNames(c("#0b3c6f", "grey85", "#7b0f0f"), 
                                      c(-1, 0, 1)),
                   rect_col_vec = setNames(c("#d6e6f2", "#f6d8cb"), 
                                           c(-1, 1)),
                   padj_threshold = 0.05) 

#p2 <- p2 + theme(axis.text.y = element_blank())
ggsave(p2, path="01_bulk-mrnaseq/results/figures", filename= "U87_volcanic_plot.pdf", width=10, height=10, dpi = 700, units = "cm")


# === Interaction ====

d <- read.csv("01_bulk-mrnaseq/results/deseq2_interaction.csv", row.names = 1)
d$GeneSymbol <- rownames(d)

summary(-log10(d$padj))
summary(d$log2FoldChange)

theme_genes <- c(
  # Stress / hypoxia / ISR
  "ATF4",      # integrated stress response (ISR)
  "P4HA1",     # collagen proline hydroxylation; ECM remodeling; often hypoxia/HIF-linked
  "ANGPTL4",   # secreted stress/hypoxia-responsive factor; vascular/lipid-handling contexts
  "CYBRD1",    # redox/iron-handling-associated stress contexts
  
  # Adhesion / cytoskeleton / migration
  "ITGAV",     # integrin alpha V; ECM adhesion / migration
  "FLNB",      # actin-binding cytoskeletal scaffold; motility/mechanics
  "TPM4",      # actin filament regulation; contractility/motility
  "EPB41L2",   # cytoskeletal organization / membrane-cortex linkage
  "SHTN1",     # cytoskeletal dynamics; neurite/migration-associated in some contexts
  "GJA1",      # Connexin43; gap junctions; cell-cell communication; migration-associated
  
  # RNA processing / translation / surveillance
  "UPF3A",     # nonsense-mediated decay (NMD)
  "UPF3B",     # nonsense-mediated decay (NMD)
  "HNRNPU",    # RNA binding; splicing/processing
  "ZRSR2",     # spliceosome component; U12-type intron splicing
  "EIF5B",     # translation initiation factor
  "EIF3J",      # translation initiation factor
  
  # Inflammatory output
  "CXCL3"    # chemokine; inflammatory/activated programs (uncomment if you want it included)
)

p3 <- volcano_plot(d,
                   pval_colname = "pvalue",
                   padj_colname = "padj",
                   logfc_colname = "log2FoldChange",
                   genesymbol_colname = "GeneSymbol",
                   title = "Interaction",
                   #ylab = NULL,
                   label_top_n = 20,
                   fc_threshold = 0,
                   ylim_min = 0,
                   ylim_max = 5,
                   xlim_min = -6.5,
                   xlim_max = 6.5,
                   ggrepel_force = 50,
                   ggrepel_max_overlaps = 5,
                   #genes_to_label = theme_genes,
                   col_vec = setNames(c("#0b3c6f", "grey85", "#7b0f0f"), 
                                      c(-1, 0, 1)),
                   rect_col_vec = setNames(c("#d6e6f2", "#f6d8cb"), 
                                           c(-1, 1)),
                   padj_threshold = 0.05)

#p3 <- p3 + theme(axis.text.y = element_blank())
ggsave(p3, path="01_bulk-mrnaseq/results/figures", filename= "Interaction_volcanic_plot.pdf", width=10, height=7, dpi = 700, units = "cm")

# === Combined plot ===

p <- p2 + plot_spacer() + p1 + plot_layout(widths = c(1, 0.025, 1))
ggsave(p, path="01_bulk-mrnaseq/results/figures", filename= "combined_volcanic_plot.pdf", width=32, height=9, dpi = 700, units = "cm")

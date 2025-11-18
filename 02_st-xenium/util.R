run_go_ora <- function(
    genes,
    universe,
    species = c("mouse", "human"),
    ont = "BP",
    p_adj = "BH",
    p_cutoff = 0.05,
    q_cutoff = 0.1
) {
  species <- match.arg(species)
  
  # Select annotation DB -------------------------------
  OrgDb <- switch(
    species,
    mouse = org.Mm.eg.db,
    human = org.Hs.eg.db
  )
  
  # Convert to Entrez ----------------------------------
  genes_entrez <- bitr(
    genes,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = OrgDb
  )$ENTREZID
  
  universe_entrez <- bitr(
    universe,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = OrgDb
  )$ENTREZID
  
  # Run ORA --------------------------------------------
  ego <- enrichGO(
    gene          = genes_entrez,
    universe      = universe_entrez,
    OrgDb         = OrgDb,
    ont           = ont,
    pAdjustMethod = p_adj,
    pvalueCutoff  = p_cutoff,
    qvalueCutoff  = q_cutoff,
    readable      = TRUE
  )
  
  return(ego)
}

run_go_gsea <- function(
    stats_vec,
    species = c("mouse", "human"),
    ont = "BP",
    n_perm = 10000,
    p_adj = "BH",
    p_cutoff = 0.05
) {
  species <- match.arg(species)
  
  # choose OrgDb -------------------------------------------------
  OrgDb <- switch(
    species,
    mouse = org.Mm.eg.db,
    human = org.Hs.eg.db
  )
  
  # Ensure vector is named --------------------------------------
  if (is.null(names(stats_vec))) {
    stop("stats_vec must be a named numeric vector of gene scores.")
  }
  
  # Convert gene symbols → Entrez and remap names ---------------
  df <- bitr(
    names(stats_vec),
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = OrgDb
  )
  
  # keep only genes that successfully mapped
  stats_vec <- stats_vec[df$SYMBOL]
  names(stats_vec) <- df$ENTREZID
  
  # Sort decreasing (required by GSEA) --------------------------
  stats_vec <- sort(stats_vec, decreasing = TRUE)
  
  # Run GSEA -----------------------------------------------------
  gsea_res <- gseGO(
    geneList     = stats_vec,
    OrgDb        = OrgDb,
    ont          = ont,
    minGSSize    = 10,
    maxGSSize    = 500,
    pvalueCutoff = p_cutoff,
    pAdjustMethod = p_adj,
    eps          = 0,          # recommended for reproducibility
    nPerm        = n_perm
  )
  
  return(gsea_res)
}

library(clusterProfiler)
library(msigdbr)

run_hallmark_gsea <- function(
    stats_vec,
    species = c("mouse", "human"),
    n_perm = 10000,
    p_adj = "BH",
    p_cutoff = 0.05
) {
  species <- match.arg(species)
  
  # Map species to msigdbr format ------------------------------------
  msig_species <- switch(
    species,
    mouse = "Mus musculus",
    human = "Homo sapiens"
  )
  
  # Choose OrgDb ------------------------------------------------------
  OrgDb <- switch(
    species,
    mouse = org.Mm.eg.db,
    human = org.Hs.eg.db
  )
  
  # Ensure named vector ----------------------------------------------
  if (is.null(names(stats_vec))) {
    stop("stats_vec must be a named numeric vector (names = gene symbols).")
  }
  
  # Get HALLMARK sets -------------------------------------------------
  hallmark_df <- msigdbr(species = msig_species, category = "H")
  
  # Convert symbols → Entrez ------------------------------------------
  conv <- bitr(
    names(stats_vec),
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = OrgDb
  )
  
  # Remap stats vector to ENTREZ --------------------------------------
  stats_vec <- stats_vec[conv$SYMBOL]
  names(stats_vec) <- conv$ENTREZID
  stats_vec <- sort(stats_vec, decreasing = TRUE)
  
  # Prepare TERM2GENE mapping -----------------------------------------
  hallmark_term2gene <- hallmark_df[, c("gs_name", "entrez_gene")]
  
  # Run GSEA -----------------------------------------------------------
  gsea_res <- GSEA(
    geneList     = stats_vec,
    TERM2GENE    = hallmark_term2gene,
    minGSSize    = 10,
    maxGSSize    = 500,
    pvalueCutoff = p_cutoff,
    pAdjustMethod = p_adj,
    verbose      = FALSE,
    eps          = 0,
    nPerm        = n_perm
  )
  
  return(gsea_res)
}

run_c8_gsea <- function(
    stats_vec,
    species = c("mouse", "human"),
    n_perm = 10000,
    p_adj = "BH",
    p_cutoff = 0.05
) {
  species <- match.arg(species)
  
  # map species to msigdbr nomenclature
  msig_species <- switch(
    species,
    mouse = "Mus musculus",
    human = "Homo sapiens"
  )
  
  # choose OrgDb for symbol → Entrez mapping
  OrgDb <- switch(
    species,
    mouse = org.Mm.eg.db,
    human = org.Hs.eg.db
  )
  
  # ensure named numeric vector
  if (is.null(names(stats_vec))) {
    stop("stats_vec must be a named vector of gene scores (names = symbols).")
  }
  
  # get C8 gene sets
  c8_df <- msigdbr(species = msig_species, category = "C8")
  
  # gene ID conversion
  conv <- bitr(
    names(stats_vec),
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = OrgDb
  )
  
  # remap stats vector to ENTREZ
  stats_vec <- stats_vec[conv$SYMBOL]
  names(stats_vec) <- conv$ENTREZID
  stats_vec <- sort(stats_vec, decreasing = TRUE)
  
  # TERM2GENE mapping
  c8_term2gene <- c8_df[, c("gs_name", "entrez_gene")]
  
  # GSEA
  gsea_res <- GSEA(
    geneList      = stats_vec,
    TERM2GENE     = c8_term2gene,
    minGSSize     = 10,
    maxGSSize     = 500,
    pvalueCutoff  = p_cutoff,
    pAdjustMethod = p_adj,
    nPerm         = n_perm,
    eps           = 0,
    verbose       = FALSE
  )
  
  return(gsea_res)
}
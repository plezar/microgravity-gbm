run_gsea_cp <- function(sig,
                        de_d,
                        stat_col = "stat",
                        pvalueCutoff = 1,
                        minGSSize = 5,
                        maxGSSize = 500,
                        seed = TRUE,
                        verbose = FALSE) {
  stopifnot(!missing(sig), !missing(de_d))
  
  # deps
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)
  requireNamespace("clusterProfiler", quietly = TRUE)
  
  # 1) TERM2GENE
  term2gene <- dplyr::bind_rows(lapply(names(sig), function(term) {
    genes <- as.character(sig[[term]])
    genes <- genes[!is.na(genes) & genes != ""]
    tibble::tibble(term = term, gene = unique(genes))
  }))
  
  # 2) ranked geneList
  if (!(stat_col %in% colnames(de_d))) {
    stop("stat_col not found in de_d: ", stat_col, "\nAvailable: ", paste(colnames(de_d), collapse = ", "))
  }
  
  geneList <- de_d %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(!is.na(.data[[stat_col]])) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    dplyr::arrange(dplyr::desc(.data[[stat_col]])) %>%
    { stats <- .[[stat_col]]; names(stats) <- .$gene; stats }
  
  geneList <- sort(geneList, decreasing = TRUE)
  
  # 3) run GSEA
  gsea_cp <- clusterProfiler::GSEA(
    geneList     = geneList,
    TERM2GENE    = term2gene,
    pvalueCutoff = pvalueCutoff,
    minGSSize    = minGSSize,
    maxGSSize    = maxGSSize,
    verbose      = verbose,
    seed         = seed
  )
  
  list(
    gsea_cp = gsea_cp,
    term2gene = term2gene,
    geneList = geneList,
    sig = sig
  )
}
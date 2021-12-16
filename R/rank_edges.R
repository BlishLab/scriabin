

#' Identify genes variable across axis of interest
#'
#' @param seu A seurat object
#' @param assay Which assay to use
#' @param slot Which slot to use
#' @param n.gene Number of variable genes to return
#' @param group.by Name of meta.data column corresponding to how dataset should be split
#' @param filter_quality Remove quality-associated genes like mitochondrial, ribosomal, etc.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' var_genes <- IDVariantGenes(seu)
#' }
IDVariantGenes <- function(seu, assay = "SCT", slot = "data", n.gene = 2000,
                           group.by = "time.orig", filter_quality = F) {
  var_genes <- AverageExpression(seu, assays = assay, slot = slot)[[assay]]
  var_genes_sd <- data.frame(x = rownames(var_genes), y = rowSds(var_genes))
  if(filter_quality) {
    bad.genes <- c(grep("^RPS", rownames(var_genes), value = T),
                   grep("^RPL", rownames(var_genes), value = T),
                   grep("^MT-", rownames(var_genes), value = T),
                   grep("^MTR", rownames(var_genes), value = T),
                   grep("MALAT1", rownames(var_genes), value = T),
                   grep("^RNA18S5", rownames(var_genes), value = T),
                   grep("^RNA28S5", rownames(var_genes), value = T),
                   grep("^ENSMNEG", rownames(var_genes), value = T))
    var_genes_sd %<>% dplyr::filter(x %notin% bad.genes)
  }
  var_genes_sd %<>% top_n(n = n.gene, wt = y) %>% pull(x)
  return(var_genes_sd)
}


#' Create single-cell gene signature
#'
#' @param seu A seurat object
#' @param variant_features Character vector of variable genes returned by IDVariantGenes
#' @param quantile_threshold Genes a distance less than this threshold will be considered part of the gene signature
#'
#' @return A list n cells long with nearest gene-cell distances for each cell
#' @export
#'
#' @examples
#' \dontrun{
#' ranked_genes <- crGeneSig(seu, variant_features = var_genes)
#' }
crGeneSig <- function(seu, variant_features = NULL, quantile_threshold = 0.1) {
  seu <- RunMCA(seu, features = variant_features)
  ds2 <- GetCellGeneRanking(seu, reduction = "mca")
  ranked_genes <- lapply(ds2, FUN = function(x) {x[x<quantile(x,quantile_threshold)]})
  return(ranked_genes)
}

PrioritizeLigands <- function(seu, gene_rankings = NULL, min.pct = 0.025, assay = "RNA", slot = "counts") {
  exprs <- GetAssayData(seu, assay = assay, slot = slot)
  expressed_genes <- rownames(exprs)[(Matrix::rowSums(exprs !=0)/ncol(exprs))>min.pct]
  background_expressed_genes <- expressed_genes %>% .[. %in% rownames(ligand_target_matrix)]
  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()
  expressed_ligands = intersect(ligands,expressed_genes)
  expressed_receptors = intersect(receptors,expressed_genes)
  potential_ligands = lr_network %>% dplyr::filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

  message("Beginning NicheNet")

  nichenet.results <- pblapply(gene_rankings, function(x) {
    activities <- predict_ligand_activities(geneset = names(x),
                                            ligand_target_matrix = ligand_target_matrix,
                                            potential_ligands = potential_ligands,
                                            background_expressed_genes = background_expressed_genes)
    best_upstream_ligands = activities %>%
      top_n(20, pearson) %>%
      arrange(-pearson) %>%
      pull(test_ligand) %>%
      unique()
    active_ligand_target_links_df = best_upstream_ligands %>%
      lapply(get_weighted_ligand_target_links,
             geneset = names(x),
             ligand_target_matrix = ligand_target_matrix, n = 200) %>%
      bind_rows %>% drop_na()
    return(list(activities,active_ligand_target_links_df))
  })
  return(nichenet.results)
}

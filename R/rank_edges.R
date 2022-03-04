

#' Identify genes variable across axis of interest
#'
#' @param seu A seurat object
#' @param assay Which assay to use
#' @param slot Which slot to use
#' @param n.gene Number of variable genes to return
#' @param group.by Name of meta.data column corresponding to how dataset should be split.
#' This corresponds to the axis of biologically interesting variation.
#' @param filter_quality Remove quality-associated genes like mitochondrial, ribosomal, etc.
#'
#' @return
#' @import genefilter Seurat dplyr
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



#' Predict ligand activity using cell-resolved gene signatures
#'
#' @param seu A seurat object
#' @param variant_genes character vector of genes to use for gene-signature calculation. This list of genes should comprise genes that are variable across the biological axis of interest.
#' Eg. for single datasets, this may simply be the dataset's HVGs. For a longitudinal dataset, these may be genes that change over time. Use `IDVariantGenes` to identify genes that vary across time, between groups, etc.
#' @param dq numeric. Distance quantile. Genes that are further away than this quantile threshold will not be considered part of a cell's gene signature and not used for ligand activity prediction
#' @param species character. Name of species from which to load ligand-receptor databases. One of: "human", "mouse", "rat". Default: "human"
#' @param database Name of ligand-receptor database to use. Default: "OmniPath"
#' When species is "human", one of: OmniPath, CellChatDB, CellPhoneDB, Ramilowski2015, Baccin2019, LRdb, Kirouac2010, ICELLNET, iTALK, EMBRACE, HPMR, Guide2Pharma, connectomeDB2020, talklr, CellTalkDB
#' When species is "mouse" or "rat", only "OmniPath" is supported.
#' To pass a custom ligand-receptor database to this function, set database = "custom"
#' @param ligands Character vector of custom ligands to use for interaction graph generation. Ignored unless database = "custom"
#' When ligands is supplied, recepts must also be supplied and equidimensional.
#' @param recepts Character vector of custom receptors to use for interaction graph generation. Ignored unless database = "custom"
#' When recepts is supplied, ligands must also be supplied and equidimensional.
#'
#' @return Returns a matrix where columns are cells, rows are potential ligands, and values are pearson coefficients corresponding to each ligand's predicted activity in that cell.
#' @references Browaeys, et al. Nat Methods (2019); Cortal, et al. Nat Biotech (2021)
#' @import nichenetr dplyr scales
#' @export
#'
#' @examples
RankActiveLigands <- function(seu, variant_genes, dq = 0.5,
                                   species = "human", database = "OmniPath",
                                   ligands = NULL, recepts = NULL, ...) {
  if(species %notin% c("human","mouse","rat") & database != "custom") {
    stop("Only human, mouse, and rat are currently supported as species\nTo use a custom ligand-receptor pair list please set database = 'custom'")
  }

  seu <- RunMCA(seu, features = variant_genes)
  ds2 <- do.call(rbind,GetCellGeneRanking(seu, reduction = "mca"))
  ds2s <- scales::rescale(ds2, from = c(min(ds2),quantile(ds2,dq)), to = c(0,1))
  ds2s[ds2s>1] <- 1
  dsp <- 1-ds2s
  if(species != "human") {
    warning("Warning: NicheNet's ligand-target matrix is built only on human observations. Use caution when extrapolating the data in this database to non-human datasets")
    colnames(dsp) <- nichenetr::convert_mouse_to_human_symbols(colnames(dsp))
    potential_ligands <- IDPotentialLigands(seu, species = species, database = database, ligands = ligands, recepts = recepts, ...)
    potential_ligands <- nichenetr::convert_mouse_to_human_symbols(potential_ligands)
  }
  else {
    potential_ligands <- IDPotentialLigands(seu, species = species, database = database, ligands = ligands, recepts = recepts, ...)
  }
  shared_targets <- intersect(rownames(ligand_target_matrix),colnames(dsp))
  shared_targets <- shared_targets[shared_targets %in% potential_ligands[[2]]]
  dsp <- t(dsp)[shared_targets,]
  ltm <- ligand_target_matrix[shared_targets,]
  message("Calculating active ligands")
  ligands_map <- potential_ligands[[1]][potential_ligands[[1]] %in% ligands_for_optim]
  preds <- cor(dsp,ltm[,colnames(ltm) %in% ligands_map])
  preds[is.na(preds)] <- 0
  if(species != "human") {
    colnames(preds) <- nichenetr::convert_human_to_mouse_symbols(rownames(dsp))
  }
  return(t(preds))
}










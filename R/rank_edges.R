

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
#' @import Seurat dplyr
#' @importFrom genefilter rowSds
#' @export
#'
#' @examples
#' \dontrun{
#' var_genes <- IDVariantGenes(seu)
#' }
IDVariantGenes <- function(seu, assay = "SCT", slot = "data", n.gene = 2000,
                           group.by = "orig.ident", filter_quality = F) {
  var_genes <- AverageExpression(seu, assays = assay, slot = slot)[[assay]]
  var_genes_sd <- data.frame(x = rownames(var_genes), y = genefilter::rowSds(var_genes))
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


#' Generate single-cell gene signature for ligand activity ranking
#'
#' @param seu A seurat object
#' @param variant_genes character vector of genes to use for gene-signature calculation. This list of genes should comprise genes that are variable across the biological axis of interest.
#' Eg. for single datasets, this may simply be the dataset's HVGs. For a longitudinal dataset, these may be genes that change over time. Use `IDVariantGenes` to identify genes that vary across time, between groups, etc.
#' @param dq numeric. Distance quantile. Genes that are further away than this quantile threshold will not be considered part of a cell's gene signature and not used for ligand activity prediction
#'
#' @return Returns a matrix where columns are cells, rows are genes, and values are binary corresponding to if that gene appears in that cell's gene signature
#' @references Cortal, et al. Nat Biotech (2021)
#' @import dplyr CelliD
#' @importFrom matrixStats rowQuantiles
#' @export
#'
#' @examples
GenerateCellSignature <- function(seu, variant_genes, dq = 0.05) {
  seu <- RunMCA(seu, features = variant_genes)
  ds2 <- t(CelliD:::GetCellGeneDistance(seu, reduction = "mca", dims = 1:30))
  rq <- matrixStats::rowQuantiles(ds2, probs = dq)
  dsp <- ds2
  dsp[rq[row(dsp)]<dsp] <- 0
  dsp[dsp>0] <- 1
  return(dsp)
}

#' Predict ligand activity using cell-resolved gene signatures
#'
#' @param seu A seurat object
#' @param signature_matrix The output of `GenerateCellSignature`
#' @param species character. Name of species from which to load ligand-receptor databases. One of: "human", "mouse", "rat". Default: "human"
#' @param database Name of ligand-receptor database to use. Default: "OmniPath"
#' When species is "human", one of: OmniPath, CellChatDB, CellPhoneDB, Ramilowski2015, Baccin2019, LRdb, Kirouac2010, ICELLNET, iTALK, EMBRACE, HPMR, Guide2Pharma, connectomeDB2020, talklr, CellTalkDB
#' When species is "mouse" or "rat", only "OmniPath" is supported.
#' To pass a custom ligand-receptor database to this function, set database = "custom"
#' @param ligands Character vector of custom ligands to use for interaction graph generation. Ignored unless database = "custom"
#' When ligands is supplied, recepts must also be supplied and equidimensional.
#' @param recepts Character vector of custom receptors to use for interaction graph generation. Ignored unless database = "custom"
#' When recepts is supplied, ligands must also be supplied and equidimensional.
#' @param potential_ligands character vector of ligands to include in active ligand ranking.
#' @param ... Additional arguments passed to IDPotentialLigands
#'
#' @return Returns a matrix where columns are cells, rows are potential ligands, and values are pearson coefficients corresponding to each ligand's predicted activity in that cell.
#' @references Browaeys, et al. Nat Methods (2019)
#' @import nichenetr dplyr scales
#' @export
#'
#' @examples
RankActiveLigands <- function(seu, signature_matrix, potential_ligands = NULL,
                              species = "human", database = "OmniPath",
                              ligands = NULL, recepts = NULL, ...) {
  if(species %notin% c("human","mouse","rat") & database != "custom") {
    stop("Only human, mouse, and rat are currently supported as species\nTo use a custom ligand-receptor pair list please set database = 'custom'")
  }
  library(nichenetr)
  dsp <- signature_matrix

  if(species != "human") {
    warning("Warning: NicheNet's ligand-target matrix is built only on human observations. Use caution when extrapolating the data in this database to non-human datasets\n")
    colnames(dsp) <- nichenetr::convert_mouse_to_human_symbols(colnames(dsp))
    potl_ligs <- IDPotentialLigands(seu, species = species, database = database, ligands = ligands, recepts = recepts, ...)

  }
  else {
    potl_ligs <- IDPotentialLigands(seu, species = species, database = database, ligands = ligands, recepts = recepts, ...)
  }
  beg <- potl_ligs[[2]]
  if(is.null(potential_ligands)) {
    potential_ligands <- potl_ligs[[1]]
  }
  shared_targets <- intersect(rownames(ligand_target_matrix),colnames(dsp))
  shared_targets <- shared_targets[shared_targets %in% beg]
  dsp <- t(dsp)[shared_targets,]
  ltm <- ligand_target_matrix[shared_targets,]
  message("Calculating active ligands")
  preds <- cor(dsp,ltm[,colnames(ltm) %in% potential_ligands])
  preds[is.na(preds)] <- 0
  if(species != "human") {
    colnames(preds) <- nichenetr::convert_human_to_mouse_symbols(colnames(preds))
  }
  return(t(preds))
}

#' Predict ligand activity using cell-resolved gene signatures
#'
#' @param seu A seurat object
#' @param signature_matrix The output of `GenerateCellSignature`
#' @param species character. Name of species from which to load ligand-receptor databases. One of: "human", "mouse", "rat". Default: "human"
#' @param database Name of ligand-receptor database to use. Default: "OmniPath"
#' When species is "human", one of: OmniPath, CellChatDB, CellPhoneDB, Ramilowski2015, Baccin2019, LRdb, Kirouac2010, ICELLNET, iTALK, EMBRACE, HPMR, Guide2Pharma, connectomeDB2020, talklr, CellTalkDB
#' When species is "mouse" or "rat", only "OmniPath" is supported.
#' To pass a custom ligand-receptor database to this function, set database = "custom"
#' @param ligands Character vector of custom ligands to use for interaction graph generation. Ignored unless database = "custom"
#' When ligands is supplied, recepts must also be supplied and equidimensional.
#' @param recepts Character vector of custom receptors to use for interaction graph generation. Ignored unless database = "custom"
#' When recepts is supplied, ligands must also be supplied and equidimensional.
#' @param potential_ligands character vector of ligands to include in active ligand ranking.
#' @param ntargets The number of top targets from the ligand-target matrix to consider as active targets of that ligand
#' @param ... Additional arguments passed to IDPotentialLigands
#'
#' @return Returns a matrix where columns are cells, rows are potential ligands, and values are pearson coefficients corresponding to each ligand's predicted activity in that cell.
#' @references Browaeys, et al. Nat Methods (2019)
#' @import nichenetr dplyr scales pbapply
#' @export
#'
#' @examples
RankLigandTargets <- function(seu, signature_matrix, potential_ligands = NULL,
                              species = "human", database = "OmniPath",
                              ligands = NULL, recepts = NULL,
                              ntargets = 250, ...) {
  if(species %notin% c("human","mouse","rat") & database != "custom") {
    stop("Only human, mouse, and rat are currently supported as species\nTo use a custom ligand-receptor pair list please set database = 'custom'")
  }
  library(nichenetr)
  dsp <- signature_matrix
  if(species != "human") {
    warning("Warning: NicheNet's ligand-target matrix is built only on human observations. Use caution when extrapolating the data in this database to non-human datasets\n")
    colnames(dsp) <- nichenetr::convert_mouse_to_human_symbols(colnames(dsp))
    potl_ligs <- IDPotentialLigands(seu, species = species, database = database, ligands = ligands, recepts = recepts, ...)
  }
  else {
    potl_ligs <- IDPotentialLigands(seu, species = species, database = database, ligands = ligands, recepts = recepts, ...)
  }
  beg <- potl_ligs[[2]]
  if(is.null(potential_ligands)) {
    potential_ligands <- potl_ligs[[1]]
  }
  shared_targets <- intersect(rownames(ligand_target_matrix),colnames(dsp))
  shared_targets <- shared_targets[shared_targets %in% beg]
  dsp <- dsp[,shared_targets]
  ltm <- ligand_target_matrix[shared_targets,]
  message("Calculating ligand-target links")
  ltl <- pblapply(seq_along(1:length(potential_ligands)), function(x) {
    ligand = potential_ligands[x]
    ltm_vec <- ligand_target_matrix[,ligand]
    top_n_score = ltm_vec %>% sort(decreasing = T) %>%
      head(ntargets) %>% min()
    ltm_vec[ltm_vec<top_n_score] <- 0
    ltm_vec <- ltm_vec[shared_targets]
    tm <- matrix(rep(ltm_vec,nrow(dsp)), ncol = ncol(dsp))
    tm[dsp==0] <- 0
    colnames(tm) <- colnames(dsp)
    rownames(tm) <- rownames(dsp)
    return(tm)
  })
  if(species != "human") {
    rownames(tm) <- nichenetr::convert_human_to_mouse_symbols(rownames(tm))
  }
  return(tm)
}

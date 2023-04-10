

#' Identify genes variable across axis of interest
#'
#' @param seu A seurat object
#' @param assay Which assay to use
#' @param slot Which slot to use
#' @param n.gene Number of variable genes to return (default: 500)
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
IDVariantGenes <- function(seu, assay = "SCT", slot = "data", n.gene = 500,
                           group.by = "orig.ident", filter_quality = F) {
  Idents(seu) <- group.by
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
                              receiver_cells = NULL,
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
  ltm <- ligand_target_matrix[shared_targets,]
  message("Calculating ligand-target links")
  if(!is.null(receiver_cells)) {
    dsp_sub <- t(dsp)[shared_targets,receiver_cells]
  } else {
    dsp_sub <- t(dsp)[shared_targets,]
  }

  ltl <- pblapply(seq_along(1:length(potential_ligands)), function(x) {
    ligand = potential_ligands[x]
    ltm_vec <- ligand_target_matrix[,ligand]
    top_n_score = ltm_vec %>% sort(decreasing = T) %>%
      head(ntargets) %>% min()
    ltm_vec[ltm_vec<top_n_score] <- 0
    ltm_vec <- ltm_vec[shared_targets]
    tm <- matrix(rep(ltm_vec,ncol(dsp_sub)), ncol = ncol(dsp_sub))
    tm[dsp_sub==0] <- 0
    colnames(tm) <- colnames(dsp_sub)
    rownames(tm) <- rownames(dsp_sub)
    if(species != "human") {
      rownames(tm) <- nichenetr::convert_human_to_mouse_symbols(rownames(tm))
    }
    return(tm)
  })
  names(ltl) <- potential_ligands
  return(ltl)
}


#' Plot most active ligands in receiver cell type
#'
#' @param seu A seurat object
#' @param active_ligands The output of `RankActiveLigands`
#' @param sender Character vector specifying the name of cells to consider as senders in the `group.by` meta.data column
#' @param receiver Character vector specifying the name of cells to consider as receivers in the `group.by` meta.data column
#' @param group.by Name of meta.data column to search for `sender` and `receiver`
#' @param pearson.threshold Pearson coefficient for calculating the percentage of receiver cells with ligand activity (default: 0.1)
#' @param ligands.display Number of top ligands to display, ranked by mean ligand activity (default: 25)
#' @param assay Assay of Seurat object to calculate average expression
#'
#' @return A ggplot object
#' @import dplyr tibble cowplot ggpubr Seurat
#' @export
#'
#' @examples
TopLigandsByIdent <- function(seu, active_ligands = NULL,
                              sender = NULL, receiver = NULL, group.by = "orig.ident",
                              pearson.threshold = 0.1, ligands.display = 25,
                              assay = "SCT") {

  receiver_cells <- colnames(seu)[seu@meta.data[,group.by]==receiver]
  ral_sub <- active_ligands[,receiver_cells]
  ral_means <- rowMeans(ral_sub)
  ral_pct <- ral_sub
  ral_pct[ral_pct>pearson.threshold] <- 1
  ral_pct[ral_pct<1] <- 0
  ral_pct <- 100*rowSums(ral_pct)/ncol(ral_sub)
  ral_res <- data.frame(ligand = names(ral_means), mean = ral_means, pct = ral_pct) %>%
    as_tibble() %>% top_n(wt = mean, n = ligands.display)


  avg_exprs <- AverageExpression(seu, features = unique(ral_res$ligand),
                                 group.by = group.by, assay = assay)[[1]][,sender]
  avg_exprs <- data.frame(ligand = names(avg_exprs), exprs = avg_exprs)
  ral_res <- merge(ral_res, avg_exprs, by = "ligand") %>% arrange(-mean)
  ral_res$ligand <- factor(ral_res$ligand, levels = ral_res$ligand)

  ggplot(ral_res, aes(x = ligand, y = mean, size = pct, color = exprs)) +
    geom_point() + scale_radius() + theme_cowplot() +
    scale_color_gradientn(colours = c("grey90","yellow","red3")) +
    ggpubr::rotate_x_text() +
    labs(x = "Ligand", y = paste0("Mean ligand activity in\n",receiver," receiver cells"),
         color = paste0("Expression by\n",sender,"\nsender cells"),
         size = paste0("Percent receiver cells\nwith activity > ",pearson.threshold))
}



#' Plot connections between ligands and predicted target genes
#' This function plots an alluvium displaying ligand-target linkages. If you consider a subset of `receiver_cells`, and a set of `ligands_of_interest`, perhaps those expressed by a particular sender cell type of interest, this function will find the top `n_targets` target genes predicted to be upregulated by those ligands. The only target genes shown are those that also appear in the individual cell's gene signature.
#'
#' @param seu A seurat object
#' @param signature_matrix The output of `GenerateCellSignature`
#' @param active_ligands The output of `RankActiveLigands`
#' @param receiver_cells Character vector of receiver cells to query
#' @param ligands_of_interest Ligands to plot
#' @param n_targets Top ligand-target gene linkages for each ligand to plot
#'
#' @return A ggplot object
#' @export
#' @import dplyr ggalluvial ggfittext
#'
#' @examples
PlotLigandTargetAlluvium <- function(seu, signature_matrix = NULL,
                                     active_ligands = active_ligands,
                                     receiver_cells = NULL, ligands_of_interest = NULL,
                                     n_targets = 250) {
  dsp <- signature_matrix
  loi <- ligands_of_interest

  lts <- reshape2::melt(ligand_target_matrix[,loi]) %>%
    group_by(Var2) %>% top_n(n = n_targets, wt = value)
  colnames(lts) <- c("target","ligand","target_weight")

  liscs <- reshape2::melt(active_ligands[loi,receiver_cells])
  colnames(liscs) <- c("ligand","cell","ligand_weight")
  liscs <- liscs %>% dplyr::mutate(cell_ligand = paste(cell,ligand)) %>%
    select(cell_ligand,ligand_weight)

  tiscs <- reshape2::melt(t(dsp)[,receiver_cells]) %>% dplyr::filter(value==1) %>%
    select(Var1,Var2)
  colnames(tiscs) <- c("target","cell")

  lts_merge <- merge(lts,tiscs,by = "target") %>%
    dplyr::mutate(cell_ligand = paste(cell,ligand)) %>%
    select(-ligand) %>% select(-cell)

  final_merge <- merge(lts_merge,liscs, by = "cell_ligand") %>%
    separate(cell_ligand, into = c("cell","ligand"), sep = " ") %>%
    dplyr::mutate(total_weight = target_weight*ligand_weight) %>%
    group_by(ligand,target) %>%
    summarise(total = sum(total_weight))

  colnames(final_merge) <- c("Ligand","Target","total")
  connectome_lodes <- to_lodes_form(final_merge, axes = 1:2)

  ggplot(connectome_lodes, aes(x = x, stratum = stratum, alluvium = alluvium,
                               y = total, label = stratum)) +
    geom_flow(color = "darkgray") +
    stat_stratum(aes(fill = stratum), width = 1/3) +
    ggfittext::geom_fit_text(aes(label = stratum),
                             stat = ggalluvial::StatStratum, width = 1/3, min.size = 3) +
    theme_cowplot() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    labs(x = NULL, y = "Regulatory weight") + NoLegend()
}


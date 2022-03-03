### This file contains functions that are currently unsupported or undergoing further testing ###


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


#' Title
#'
#' @param seu
#' @param gene_rankings
#' @param min.pct
#' @param assay
#' @param slot
#'
#' @return
#' @export
#'
#' @examples
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



#' Title
#'
#' @param seu
#' @param by
#' @param name
#' @param split.by
#'
#' @return
#' @import reshape2 pbapply
#' @export
#'
#' @examples
AssembleInteractionGraphs <- function(seu, by = "weighted", split.by = "time.orig") {
  seu_split <- SplitObject(seu, split.by = split.by)
  if(by=="weighted") {
    seu_split <- pblapply(seu_split, BuildWeightedInteraction)
    interaction_graphs <- pblapply(seu_split, function(x) {
      p <- as.matrix(x@graphs$weighted_interaction)
      merge_ids <- seu$bins
      ids_name <- merge_ids[names(merge_ids) %in% colnames(p)]
      identical(names(ids_name),colnames(p))
      p_summ <- apply(p,1,FUN = function(x) {aggregate(x,by = list(ids_name),FUN = median)})
      p_summ <- bind_rows(p_summ, .id = "column_label")
      p_summ <- reshape2::dcast(p_summ,column_label~Group.1,value.var="x") %>% column_to_rownames("column_label")
      p_summ <- apply(p_summ,2,FUN = function(x) {aggregate(x,by = list(ids_name),FUN = median)})
      p_summ <- bind_rows(p_summ, .id = "column_label")
      p_summ <- reshape2::dcast(p_summ,column_label~Group.1,value.var="x") %>% column_to_rownames("column_label")
      p_summ <- t(p_summ[str_sort(rownames(p_summ),numeric=T),str_sort(colnames(p_summ),numeric=T)])
      return(p_summ)
    })
  }
  if(by=="prior") {
    seu_split <- pblapply(seu_split, BuildPriorInteraction)
    interaction_graphs <- pblapply(seu_split, function(x) {
      p <- as.matrix(x@graphs$prior_interaction)
      merge_ids <- seu$bins
      ids_name <- merge_ids[names(merge_ids) %in% colnames(p)]
      identical(names(ids_name),colnames(p))
      p_summ <- apply(p,1,FUN = function(x) {aggregate(x,by = list(ids_name),FUN = median)})
      p_summ <- bind_rows(p_summ, .id = "column_label")
      p_summ <- reshape2::dcast(p_summ,column_label~Group.1,value.var="x") %>% column_to_rownames("column_label")
      p_summ <- apply(p_summ,2,FUN = function(x) {aggregate(x,by = list(ids_name),FUN = median)})
      p_summ <- bind_rows(p_summ, .id = "column_label")
      p_summ <- reshape2::dcast(p_summ,column_label~Group.1,value.var="x") %>% column_to_rownames("column_label")
      p_summ <- t(p_summ[str_sort(rownames(p_summ),numeric=T),str_sort(colnames(p_summ),numeric=T)])
      return(p_summ)
    })
  }

  return(interaction_graphs)
}



#' Title
#'
#' @param object
#' @param assay
#' @param slot
#' @param database
#' @param ligands
#' @param recepts
#' @param specific
#' @param ranked_genes
#' @param cell.type.calls
#' @param dims_umap
#' @param dims_neighbors
#' @param cluster_resolution
#'
#' @return
#' @import dplyr
#' @export
#'
#' @examples
scInteraction <- function(object, assay = "SCT", slot = "data",
                          database = "OmniPath", ligands = NULL, recepts = NULL,
                          specific = F, ranked_genes = NULL,
                          cell.type.calls = "celltype.l2", dims_umap = 1:50,
                          dims_neighbors = 1:50, cluster_resolution = 0.1) {
  if(cell.type.calls %notin% colnames(object@meta.data)) {
    stop("Argument cell.type.calls not found in provided Seurat object")
  }
  if(database=="custom") {
    message("Using custom database")
    ligands <- ligands
    recepts <- recepts
    lit.put <- data.frame(pair.name = paste(ligands,recepts,sep="_"), ligands = ligands, recepts = recepts)
  }
  else {
    all <- readRDS(system.file(package = "scriabin", "lr_resources.rds"))
    if(database %notin% names(all)) {
      stop("Database must be one of: OmniPath, CellChatDB, CellPhoneDB, Ramilowski2015, Baccin2019, LRdb, Kirouac2010, ICELLNET, iTALK, EMBRACE, HPMR, Guide2Pharma, connectomeDB2020, talklr, CellTalkDB")
    }
    message(paste("Using database",database))
    pairs <- as.data.frame(all[[database]][,c("source_genesymbol","target_genesymbol")] %>% mutate_all(as.character))
    lit.put <- pairs %>% dplyr::mutate(pair = paste(source_genesymbol,target_genesymbol, sep = "_"))
    lit.put <- as.data.frame(lit.put[,c("pair","source_genesymbol","target_genesymbol")])
    ligands <- as.character(lit.put[, "source_genesymbol"])
    recepts <- as.character(lit.put[, "target_genesymbol"])
  }
  ligands.use <- intersect(ligands, rownames(object@assays[[assay]]))
  recepts.use <- intersect(recepts, rownames(object@assays[[assay]]))
  genes.use = union(ligands.use, recepts.use)

  if(specific) {
    message("Only considering genes in per-cell gene signature")
    ranked_names <- lapply(ranked_genes, function(x) {
      names(x)
    })
    ranked_mat <- as.matrix(reshape2::dcast(reshape2::melt(t(bind_rows(ranked_names))), formula = value~Var1) %>% column_to_rownames("value"))
    genes.use <- intersect(genes.use,rownames(ranked_mat))
    ranked_mat <- ranked_mat[genes.use,]
    cell.exprs <- GetAssayData(object, assay = assay, slot = slot)[genes.use,]
    cell.exprs[is.na(ranked_mat)] <- 0
    cell.exprs <- as.data.frame(cell.exprs) %>% rownames_to_column(var = "gene")
  }
  else {
    cell.exprs <- as.data.frame(GetAssayData(object, assay = assay, slot = slot)[genes.use,]) %>% rownames_to_column(var = "gene")
  }

  ligands.df <- data.frame(ligands)
  ligands.df$id <- 1:nrow(ligands.df)
  recepts.df <- data.frame(recepts)
  recepts.df$id <- 1:nrow(recepts.df)
  cell.exprs.rec <- merge(recepts.df, cell.exprs,
                          by.x = "recepts", by.y = "gene", all.x = T)
  cell.exprs.rec <- cell.exprs.rec[order(cell.exprs.rec$id),
                                   ]
  cell.exprs.lig <- merge(ligands.df, cell.exprs,
                          by.x = "ligands", by.y = "gene", all.x = T)
  cell.exprs.lig <- cell.exprs.lig[order(cell.exprs.lig$id),
                                   ]
  sources <- colnames(object)
  targets <- colnames(object)

  message(paste("\nGenerating Interaction Matrix..."))

  a <- as.matrix(cell.exprs.lig[,3:ncol(cell.exprs.lig)])
  a[is.na(a)] <- 0
  b <- as.matrix(cell.exprs.rec[,3:ncol(cell.exprs.rec)])
  b[is.na(b)] <- 0

  m <- as.sparse(sqrt(t(pbsapply(1:nrow(a), function(i) tcrossprod(a[i, ], b[i, ])))))
  rownames(m) <- paste(cell.exprs.lig$ligands, cell.exprs.rec$recepts, sep = "=")
  cna <- rep(colnames(object),ncol(object))
  cnb <- rep(colnames(object),each=ncol(object))
  colnames(m) <- paste(cna,cnb,sep = "=")

  message("Creating single-cell interaction object")

  edge_seu <- CreateSeuratObject(m)
  edge_assay <- CreateAssayObject(data = m)
  edge_seu[["edge"]] <- edge_assay
  DefaultAssay(edge_seu) <- "edge"
  edge_seu <- ScaleData(edge_seu)
  edge_seu@assays$RNA <- NULL

  message("Running PCA")
  edge_seu <- RunPCA(edge_seu, assay = "edge", features = rownames(edge_seu))
  message("Running UMAP")
  edge_seu <- RunUMAP(edge_seu, assay = "edge", dims = dims_umap)
  message("Constructing neighbor graphs")
  edge_seu <- FindNeighbors(edge_seu, assay = "edge", dims = dims_neighbors)
  message("Clustering")
  edge_seu <- FindClusters(edge_seu, assay = "edge", resolution = cluster_resolution)

  message("Mapping metadata")
  edge_seu$sender <- sub("\\=.*", "", colnames(edge_seu))
  edge_seu$receiver <- sub('.*=', '', colnames(edge_seu))

  edge_seu$sender_ct <- mapvalues(edge_seu$sender, from = colnames(object),
                                  to = as.character(object@meta.data[,cell.type.calls]),
                                  warn_missing = F)
  edge_seu$receiver_ct <- mapvalues(edge_seu$receiver, from = colnames(object),
                                    to = as.character(object@meta.data[,cell.type.calls]),
                                    warn_missing = F)

  return(edge_seu)
}





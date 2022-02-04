
#' Title
#'
#' @param object
#' @param assay
#' @param slot
#' @param database
#' @param ligands
#' @param recepts
#' @param senders
#' @param receivers
#'
#' @return
#' @import dplyr pbapply
#' @export
#'
#' @examples
GenerateCCIM <- function(object, assay = "SCT", slot = "data",
                         database = "OmniPath", ligands = NULL,
                         recepts = NULL, senders = NULL, receivers = NULL) {
  #code to connect with l-r databases adapted from Connectome (Raredon, et al.)
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

  if(is.null(senders) & is.null(receivers)){
    senders <- receivers <- colnames(object)
    cell.exprs <- as.data.frame(GetAssayData(object, assay = assay, slot = slot)[genes.use,]) %>% rownames_to_column(var = "gene")
    ligands.df <- data.frame(ligands)
    ligands.df$id <- 1:nrow(ligands.df)
    recepts.df <- data.frame(recepts)
    recepts.df$id <- 1:nrow(recepts.df)
    cell.exprs.lig <- merge(ligands.df, cell.exprs,
                            by.x = "ligands", by.y = "gene", all.x = T)
    cell.exprs.lig <- cell.exprs.lig[order(cell.exprs.lig$id),
                                     ]
    cell.exprs.rec <- merge(recepts.df, cell.exprs,
                            by.x = "recepts", by.y = "gene", all.x = T)
    cell.exprs.rec <- cell.exprs.rec[order(cell.exprs.rec$id),
                                     ]

  }
  else {
    if(is.null(senders)) {
      senders <- colnames(object)
    }
    if(is.null(receivers)) {
      receivers <- colnames(object)
    }
    cell.exprsl <- as.data.frame(GetAssayData(object, assay = assay, slot = slot)[genes.use,senders]) %>% rownames_to_column(var = "gene")
    cell.exprsr <- as.data.frame(GetAssayData(object, assay = assay, slot = slot)[genes.use,receivers]) %>% rownames_to_column(var = "gene")
    ligands.df <- data.frame(ligands)
    ligands.df$id <- 1:nrow(ligands.df)
    recepts.df <- data.frame(recepts)
    recepts.df$id <- 1:nrow(recepts.df)
    cell.exprs.lig <- merge(ligands.df, cell.exprsl,
                            by.x = "ligands", by.y = "gene", all.x = T)
    cell.exprs.lig <- cell.exprs.lig[order(cell.exprs.lig$id),
                                     ]
    cell.exprs.rec <- merge(recepts.df, cell.exprsr,
                            by.x = "recepts", by.y = "gene", all.x = T)
    cell.exprs.rec <- cell.exprs.rec[order(cell.exprs.rec$id),
                                     ]
  }

  a <- as.matrix(cell.exprs.lig[,3:ncol(cell.exprs.lig)])
  a[is.na(a)] <- 0
  b <- as.matrix(cell.exprs.rec[,3:ncol(cell.exprs.rec)])
  b[is.na(b)] <- 0

  message(paste("Calculating CCIM between",length(senders),"senders and",length(receivers),"receivers"))
  message(paste("\nGenerating Interaction Matrix..."))
  m <- sqrt(as.sparse((pbsapply(1:nrow(a), function(i) tcrossprod(a[i, ], b[i, ])))))
  colnames(m) <- paste(cell.exprs.lig$ligands, cell.exprs.rec$recepts, sep = "=")
  cna <- rep(senders,length(receivers))
  cnb <- rep(receivers,each=length(senders))
  rownames(m) <- paste(cna,cnb,sep = "=")
  seu <- CreateSeuratObject(counts = t(m), assay = "CCIM")
  return(seu)
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
#' @param correct.depth
#' @param graph_name
#'
#' @return
#' @import dplyr
#' @export
#'
#' @examples
BuildPriorInteraction <- function (object, assay = "SCT", slot = "data",
                                   database = "OmniPath", ligands = NULL, recepts = NULL,
                                   specific = F, ranked_genes = NULL,
                                   correct.depth = T, graph_name = "prior_interaction") {
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
  m <- crossprod(sqrt(a),sqrt(b))

  # m <- outer (
  #   cell.exprs.lig[,3:ncol(cell.exprs.lig)],     # First dimension:  the rows     (x)
  #   cell.exprs.rec[,3:ncol(cell.exprs.rec)],     # Second dimension: the columns  (y)
  #   Vectorize(function (x, y)   sum(x*y,na.rm=T))
  #   ## Rows represent "senders" and columns represent "receivers". Eg. m[1,2] is the interaction potential between cell #1 as a sender and cell #2 as a receiver.
  # )

  if(correct.depth) {
    message("Correcting for sequencing depth . . . ")
    depth <- outer(object$nCount_RNA,object$nCount_RNA,FUN = "+")
    lm_data <- data.frame(score=as.vector(m),depth=as.vector(depth))
    model <- lm(score~depth, data = lm_data)

    results <- matrix(model$residuals, nrow = ncol(m))
    dimnames(results) <- dimnames(m)

    results <- as.Graph(results)
    DefaultAssay(results) <- assay
  }
  else {
    results <- as.Graph(m)
    DefaultAssay(results) <- assay
  }

  object[[graph_name]] <- results
  return(object)
}





#' Title
#'
#' @param object
#' @param nichenet_results
#' @param assay
#' @param slot
#' @param pearson.cutoff
#' @param scale.factors
#' @param database
#' @param ligands
#' @param recepts
#' @param correct.depth
#' @param graph_name
#'
#' @return
#' @import dplyr tidyft tibble
#' @export
#'
#' @examples
BuildWeightedInteraction <- function (object, nichenet_results = late1.nnr, assay = "SCT", slot = "data",
                                      pearson.cutoff = 0.1, scale.factors = c(1.5,3),
                                      database = "OmniPath", ligands = NULL, recepts = NULL,
                                      correct.depth = T, graph_name = "weighted_interaction") {
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
  cell.exprs <- as.data.frame(GetAssayData(object, assay = assay, slot = slot)[genes.use,]) %>% rownames_to_column(var = "gene")
  ligands.df <- data.frame(lit.put[, c(1,2)]) %>% mutate_all(as.character)
  colnames(ligands.df) <- c("pair","ligands")
  ligands.df$id <- 1:nrow(ligands.df)
  recepts.df <- data.frame(lit.put[, c(1,3)]) %>% mutate_all(as.character)
  colnames(recepts.df) <- c("pair","recepts")
  recepts.df$id <- 1:nrow(recepts.df)
  cell.exprs.rec <- merge(recepts.df, cell.exprs,
                          by.x = "recepts", by.y = "gene", all.x = T)
  cell.exprs.rec <- cell.exprs.rec[order(cell.exprs.rec$id),
                                   ]
  cell.exprs.lig <- merge(ligands.df, cell.exprs,
                          by.x = "ligands", by.y = "gene", all.x = T)
  cell.exprs.lig <- cell.exprs.lig[order(cell.exprs.lig$id),
                                   ]
  message("Weighting interaction matrix")
  nichenet_results <- nichenet_results[names(nichenet_results) %in% colnames(object)]
  nichenet_results <- lapply(nichenet_results, FUN = function(x) {
    results <- x[[1]] %>% pull(pearson)
    names(results) <- x[[1]] %>% pull(test_ligand)
    return(results)
  })
  test <- bind_rows(nichenet_results)
  rownames(test) <- names(nichenet_results)
  test <- reshape2::melt(test %>% rownames_to_column(var = "cell")) %>% dplyr::filter(value>pearson.cutoff)
  colnames(test) <- c("cell","ligand","pearson")
  test$weight_factor <- scales::rescale(test$pearson, scale.factors)
  int.to.merge <- lit.put[,c(1,2,3)]
  colnames(int.to.merge) <- c("pair","ligand","recepts")
  test <- merge(test,int.to.merge,by = "ligand", all.x=T)
  test <- test[!is.na(test$recepts),]
  test2 <- reshape2::dcast(test, recepts~cell, value.var = "weight_factor", fun.aggregate = sum)
  test3 <- merge(cell.exprs.rec[,1:3],test2,all.x = T) %>% arrange(id)

  cell.exprs.rec.m <- cell.exprs.rec[,4:ncol(cell.exprs.rec)]
  weight.m <- as.matrix(test3[,4:ncol(test3)])
  weight.m[weight.m==0] <- 1
  weight.m[is.na(weight.m)] <- 1
  cells.bind <- colnames(cell.exprs.rec.m)[colnames(cell.exprs.rec.m) %notin% colnames(weight.m)]
  tobind <- matrix(1,nrow = nrow(weight.m),ncol = length(cells.bind))
  colnames(tobind) <- cells.bind
  weight.m <- cbind(weight.m,tobind)
  weight.m <- weight.m[,match(colnames(cell.exprs.rec.m),colnames(weight.m))]

  final <- weight.m * as.matrix(cell.exprs.rec.m)
  weighted.rec <- cbind(cell.exprs.rec.m[,1:3],final)

  message(paste("\nGenerating Interaction Matrix..."))
  a <- as.matrix(cell.exprs.lig[,4:ncol(cell.exprs.lig)])
  a[is.na(a)] <- 0
  b <- as.matrix(weighted.rec[,4:ncol(weighted.rec)])
  b[is.na(b)] <- 0
  m <- crossprod(sqrt(a),sqrt(b))

  #
  # m <- outer (
  #   cell.exprs.lig[,4:ncol(cell.exprs.lig)],     # First dimension:  the rows     (x)
  #   weighted.rec[,4:ncol(weighted.rec)],     # Second dimension: the columns  (y)
  #   Vectorize(function (x, y)   sum(x*y,na.rm=T))
  #   ## Rows represent "senders" and columns represent "receivers". Eg. m[1,2] is the interaction potential between cell #1 as a sender and cell #2 as a receiver.
  # )

  if(correct.depth) {
    message("Correcting for sequencing depth . . . ")
    depth <- outer(object$nCount_RNA,object$nCount_RNA,FUN = "+")
    lm_data <- data.frame(score=as.vector(m),depth=as.vector(depth))
    model <- lm(score~depth, data = lm_data)

    results <- matrix(model$residuals, nrow = ncol(m))
    dimnames(results) <- dimnames(m)

    results <- as.Graph(results)
    DefaultAssay(results) <- assay
  }
  else {
    results <- as.Graph(m)
    DefaultAssay(results) <- assay
  }

  object[[graph_name]] <- results
  return(object)
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
AssembleInteractionGraphs <- function(seu, by = "weighted", name = "rv", split.by = "time.orig") {
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


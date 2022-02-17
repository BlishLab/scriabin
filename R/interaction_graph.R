
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
  m <- m[,Matrix::colSums(m)>0]
  return(CreateSeuratObject(counts = Matrix::t(m), assay = "CCIM"))
}

#' Title
#'
#' @param ccim_seu
#' @param seu
#' @param columns_map
#'
#' @return
#' @import stringr
#' @export
#'
#' @examples
MapMetaData <- function(ccim_seu, seu, columns_map = NULL) {
  if(is.null(columns_map)) {
    message("Mapping all metadata columns")
    columns_map <- colnames(seu@meta.data)
  }
  ccim_seu$sender <- stringr::word(colnames(ccim_seu),1,sep = "=")
  ccim_seu$receiver <- stringr::word(colnames(ccim_seu),2,sep = "=")
  #convert any factors to character
  classes <- unlist(lapply(seu@meta.data,class))[columns_map]
  col_convert <- names(classes)[classes=="factor"]
  for(i in 1:length(col_convert)) {
    seu@meta.data[,col_convert[i]] <- as.character(seu@meta.data[,col_convert[i]])
  }
  for (i in 1:length(columns_map)) {
    ccim_seu@meta.data[,paste("sender",columns_map[i],sep = "_")] <- scriabin::mapvalues(ccim_seu$sender, from = colnames(seu), to = seu@meta.data[,columns_map[i]], warn_missing = F)
    ccim_seu@meta.data[,paste("receiver",columns_map[i],sep = "_")] <- scriabin::mapvalues(ccim_seu$receiver, from = colnames(seu), to = seu@meta.data[,columns_map[i]], warn_missing = F)
  }
  return(ccim_seu)
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



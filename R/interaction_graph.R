
#' Calculate cell-cell interaction matrix
#'
#' @param object A seurat object
#' @param assay Assay in Seurat object from which to pull expression values
#' @param slot Slot within assay from which to pull expression values
#' @param species character. Name of species from which to load ligand-receptor databases. One of: "human", "mouse", "rat". Default: "human"
#' @param database Name of ligand-receptor database to use. Default: "OmniPath"
#' When species is "human", one of: OmniPath, CellChatDB, CellPhoneDB, Ramilowski2015, Baccin2019, LRdb, Kirouac2010, ICELLNET, iTALK, EMBRACE, HPMR, Guide2Pharma, connectomeDB2020, talklr, CellTalkDB
#' When species is "mouse" or "rat", only "OmniPath" is supported.
#' To pass a custom ligand-receptor database to this function, set database = "custom"
#' @param ligands Character vector of custom ligands to use for interaction graph generation. Ignored unless database = "custom"
#' When ligands is supplied, recepts must also be supplied and equidimensional.
#' @param recepts Character vector of custom receptors to use for interaction graph generation. Ignored unless database = "custom"
#' When recepts is supplied, ligands must also be supplied and equidimensional.
#' @param senders Character vector of cell names to be used as sender cells for interaction graph generation. By default, all cells in object are used as senders.
#' @param receivers Character vector of cell naems to be used as receiver cells for interaction graph generation. By default, all cells in object are used as receivers.
#' @param weighted logical. Weight the cell-cell interaction graph by ligands predicted to be active by NicheNet? Default: FALSE
#' @param nichenet_results List or matrix. Predicted ligand activities using Scriabin's implementation of NicheNet in RankActiveLigands
#' @param pearson.cutoff numeric. Threshold for determining which ligand activities are "active". Ligands below this threshold will be considered inactive and not used for weighting. Default: 0.075
#' @param scale.factors numeric. Determines the magnitude of ligand and receptor expression values weighting by their predicted activities. Given scale factors of c(x,y), ligand-receptor pairs where the ligand's pearson = pearson.cutoff will be weighted by a factor of x, and the ligand with the highest pearson will be weighted by a factor of y. Default: c(1.5,3).
#' @param weight.method One of "sum" or "product". Method to use for weighting CCIM by ligand activities. "Sum" sums receptor expression values with ligand activities passing the `pearson.cutoff` threshold. Thus, a zero value for receptor expression may still result in a positive value for the ligand-receptor pair mechanism. "Product" takes the product of the receptor expression value and the predicted ligand activities. A zero value in the receptor expression will always result in a zero value for that ligand-receptor pair.
#'
#' @return Returns a Seurat object with assay "CCIM" where columns are cell-cell pairs and rows are ligand-receptor pairs. "Expression" values are the geometric mean expression value between each pair of sender and receiver cells. By default, names ligand-receptor pairs and cell-cell pairs are separated by "="
#' @import dplyr pbapply Matrix Seurat
#' @importFrom tibble rownames_to_column
#' @importFrom tibble column_to_rownames
#' @references Browaeys, et al. Nature Methods (2019)
#' @export
#'
#' @examples
GenerateCCIM <- function(object, assay = "SCT", slot = "data",
                         species = "human", database = "OmniPath",
                         ligands = NULL, recepts = NULL,
                         senders = NULL, receivers = NULL,
                         weighted = F, nichenet_results = NULL,
                         pearson.cutoff = 0.075, scale.factors = c(1.5,3),
                         weight.method = "sum") {
  #code to connect with l-r databases adapted from Connectome (Raredon, et al.)
  if(database=="custom") {
    if(is.null(ligands) | is.null(recepts)) {
      stop("To use custom database, please supply equidimensional character vectors of ligands and recepts")
    }
    message("Using custom database")
    ligands <- ligands
    recepts <- recepts
    lit.put <- data.frame(pair.name = paste(ligands,recepts,sep="_"), source_genesymbol = ligands, target_genesymbol = recepts)
  } else if((!is.null(ligands) | !is.null(recepts)) & database != "custom") {
    stop("To use custom ligand or receptor lists, set database = 'custom'")
  } else {
    lit.put <- scriabin::LoadLR(species = species, database = database)
    ligands <- as.character(lit.put[, "source_genesymbol"])
    recepts <- as.character(lit.put[, "target_genesymbol"])
  }

  ligands.use <- intersect(ligands, rownames(object@assays[[assay]]))
  recepts.use <- intersect(recepts, rownames(object@assays[[assay]]))
  genes.use = union(ligands.use, recepts.use)

  lit.put <- lit.put %>%
    dplyr::filter(source_genesymbol %in% ligands.use) %>%
    dplyr::filter(target_genesymbol %in% recepts.use)
  ligands <- as.character(lit.put[, "source_genesymbol"])
  recepts <- as.character(lit.put[, "target_genesymbol"])

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

  if(weighted) {
    if(is.null(nichenet_results)) {
      stop("NicheNet results must be included to perform CCIM weighting")
    }
    message("Weighting interaction graph")
    if(class(nichenet_results) %notin% c("list","matrix")) {
      stop("NicheNet results must be supplied as one of list or matrix")
    }
    if(class(nichenet_results)=="list") {
      nichenet_results <- nichenet_results[names(nichenet_results) %in% receivers]
      nichenet_results <- bind_rows(lapply(nichenet_results, FUN = function(x) {
        results <- x[[1]] %>% pull(pearson)
        names(results) <- x[[1]] %>% pull(test_ligand)
        return(results)
      }), .id = "cell")
      nnm <- reshape2::melt(nichenet_results) %>% dplyr::filter(value>pearson.cutoff)
      colnames(nnm) <- c("cell","ligand","pearson")
    }
    if(class(nichenet_results)[[1]]=="matrix") {
      nichenet_results <- nichenet_results[,receivers]
      nnm <- reshape2::melt(nichenet_results) %>% dplyr::filter(value>pearson.cutoff)
      colnames(nnm) <- c("ligand","cell","pearson")
    }
    if(nrow(nnm)==0) {
      warning("No active ligands to weight. Please check specified pearson cutoff")
      message("Using unweighted ligand-receptor matrices")
      a <- as.matrix(cell.exprs.lig[,3:ncol(cell.exprs.lig)])
      a[is.na(a)] <- 0
      b <- as.matrix(cell.exprs.rec[,3:ncol(cell.exprs.rec)])
      b[is.na(b)] <- 0
    }
    else {
      message(paste("Found active ligands to weight. Weighting with pearson cutoff ",pearson.cutoff))
      if(weight.method %notin% c("product","sum")) {
        stop("weight.method must be one of product or sum")
      }
      if(weight.method=="product") {
        message("Using product method to weight CCIM")
        nnm$weight_factor <- scales::rescale(nnm$pearson, scale.factors)
        int.to.merge <- lit.put[,c(1,2,3)]
        colnames(int.to.merge) <- c("pair","ligand","recepts")
        nnm <- merge(nnm,int.to.merge,by = "ligand", all.x=T)
        nnm <- nnm[!is.na(nnm$recepts),]
        test2 <- reshape2::dcast(nnm, recepts~cell, value.var = "weight_factor", fun.aggregate = sum)
        test3 <- merge(cell.exprs.rec[,1:2],test2,all.x = T) %>% arrange(id)

        cell.exprs.rec.m <- cell.exprs.rec[,3:ncol(cell.exprs.rec)]
        weight.m <- as.matrix(test3[,3:ncol(test3)])
        weight.m[weight.m==0] <- 1
        weight.m[is.na(weight.m)] <- 1
        cells.bind <- colnames(cell.exprs.rec.m)[colnames(cell.exprs.rec.m) %notin% colnames(weight.m)]
        tobind <- matrix(1,nrow = nrow(weight.m),ncol = length(cells.bind))
        colnames(tobind) <- cells.bind
        weight.m <- cbind(weight.m,tobind)
        weight.m <- weight.m[,match(colnames(cell.exprs.rec.m),colnames(weight.m))]

        final <- weight.m * as.matrix(cell.exprs.rec.m)
        # weighted.rec <- cbind(cell.exprs.rec.m[,1:3],final)

        a <- as.matrix(cell.exprs.lig[,3:ncol(cell.exprs.lig)])
        a[is.na(a)] <- 0
        b <- as.matrix(final)
        b[is.na(b)] <- 0
      }
      if(weight.method=="sum") {
        message("Using sum method to weight CCIM")
        nnm$weight_factor <- scales::rescale(nnm$pearson, scale.factors)
        int.to.merge <- lit.put[,c(1,2,3)]
        colnames(int.to.merge) <- c("pair","ligand","recepts")
        nnm <- merge(nnm,int.to.merge,by = "ligand", all.x=T)
        nnm <- nnm[!is.na(nnm$recepts),]
        test2 <- reshape2::dcast(nnm, recepts~cell, value.var = "weight_factor", fun.aggregate = sum)
        test3 <- merge(cell.exprs.rec[,1:2],test2,all.x = T) %>% arrange(id)

        cell.exprs.rec.m <- cell.exprs.rec[,3:ncol(cell.exprs.rec)]
        weight.m <- as.matrix(test3[,3:ncol(test3)])
        # weight.m[weight.m==0] <- 1
        weight.m[is.na(weight.m)] <- 0
        cells.bind <- colnames(cell.exprs.rec.m)[colnames(cell.exprs.rec.m) %notin% colnames(weight.m)]
        tobind <- matrix(0,nrow = nrow(weight.m),ncol = length(cells.bind))
        colnames(tobind) <- cells.bind
        weight.m <- cbind(weight.m,tobind)
        weight.m <- weight.m[,match(colnames(cell.exprs.rec.m),colnames(weight.m))]

        final <- weight.m + as.matrix(cell.exprs.rec.m)
        # weighted.rec <- cbind(cell.exprs.rec.m[,1:3],final)

        a <- as.matrix(cell.exprs.lig[,3:ncol(cell.exprs.lig)])
        a[is.na(a)] <- 0
        b <- as.matrix(final)
        b[is.na(b)] <- 0
      }

    }
  }


  else {
    message("Using unweighted ligand-receptor matrices")
    a <- as.matrix(cell.exprs.lig[,3:ncol(cell.exprs.lig)])
    a[is.na(a)] <- 0
    b <- as.matrix(cell.exprs.rec[,3:ncol(cell.exprs.rec)])
    b[is.na(b)] <- 0
  }

  message(paste("Calculating CCIM between",length(senders),"senders and",length(receivers),"receivers"))
  message(paste("\nGenerating Interaction Matrix..."))
  m <- sqrt(as.sparse((pbsapply(1:nrow(a), function(i) tcrossprod(a[i, ], b[i, ])))))

  colnames(m) <- paste(cell.exprs.lig$ligands, cell.exprs.rec$recepts, sep = "=")
  cna <- rep(senders,length(receivers))
  cnb <- rep(receivers,each=length(senders))

  rownames(m) <- paste(cna,cnb,sep = "=")
  m <- m[,Matrix::colSums(m)>0]
  seu <- CreateSeuratObject(counts = Matrix::t(m), assay = "CCIM")
  seu <- MapMetaData(ccim_seu = seu, seu = object)
  return(seu)
}


#' Map metadata from Seurat object to CCIM object
#'
#' @param ccim_seu Seurat object generated by GenerateCCIM that contains a cell-cell interaction matrix
#' @param seu Parent Seurat object from which the cell-cell interaction matrix was generated
#' @param columns_map Character vector of names of metadata columns to map from parent Seurat object
#'
#' @return Returns a Seurat object with metadata from each sender and receiver cell mapped from the parent Seurat object. Metadata about the sender cell in each cell-cell pair is prepended with "sender_" and metadata about the receiver cell prepended with "receiver_"
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
  if(length(col_convert)>0) {
    for(i in 1:length(col_convert)) {
      seu@meta.data[,col_convert[i]] <- as.character(seu@meta.data[,col_convert[i]])
    }
  }

  for (i in 1:length(columns_map)) {
    ccim_seu@meta.data[,paste("sender",columns_map[i],sep = "_")] <- scriabin::mapvalues(ccim_seu$sender, from = colnames(seu), to = seu@meta.data[,columns_map[i]], warn_missing = F)
    ccim_seu@meta.data[,paste("receiver",columns_map[i],sep = "_")] <- scriabin::mapvalues(ccim_seu$receiver, from = colnames(seu), to = seu@meta.data[,columns_map[i]], warn_missing = F)
  }
  return(ccim_seu)
}


#' Build unweighted summarized interaction graph
#'
#' @param object A Seurat object
#' @param assay Assay in Seurat object from which to pull expression values
#' @param slot Slot within assay from which to pull expression values
#' @param species Name of species for which to pull ligand-receptor interactions. One of "human", "mouse", or "rat"
#' @param database Name of ligand-receptor database to use. Default: "OmniPath"
#' When species is "human", one of: OmniPath, CellChatDB, CellPhoneDB, Ramilowski2015, Baccin2019, LRdb, Kirouac2010, ICELLNET, iTALK, EMBRACE, HPMR, Guide2Pharma, connectomeDB2020, talklr, CellTalkDB
#' When species is "mouse" or "rat", only "OmniPath" is supported.
#' To pass a custom ligand-receptor database to this function, set database = "custom"
#' @param ligands Character vector of custom ligands to use for interaction graph generation. Ignored unless database = "custom"
#' When ligands is supplied, recepts must also be supplied and equidimensional.
#' @param recepts Character vector of custom receptors to use for interaction graph generation. Ignored unless database = "custom"
#' When recepts is supplied, ligands must also be supplied and equidimensional.
#' @param specific logical. When TRUE, consider only the genes in each cell's predefined gene signature (see crGeneSig) as expressed. Default FALSE
#' @param ranked_genes Cell-resolved gene signatures, used only when specific = T
#' @param correct.depth Correct summarized interaction graph for sequencing depth by linear regression. The sequencing depth of a cell-cell pair is the sum of UMI counts for each cell
#' @param graph_name Name of summarized interaction graph to place into output. Default "prior_interaction"
#'
#' @return Returns a Seurat object with an unweighted summarized interaction graph in the Graphs slot
#' @import dplyr Seurat
#' @importFrom tibble rownames_to_column
#' @importFrom tibble column_to_rownames
#' @export
#'
#' @examples
BuildPriorInteraction <- function (object, assay = "SCT", slot = "data",
                                   species = "human", database = "OmniPath",
                                   ligands = NULL, recepts = NULL,
                                   specific = F, ranked_genes = NULL,
                                   correct.depth = T, graph_name = "prior_interaction") {
  if(database=="custom") {
    if(is.null(ligands) | is.null(recepts)) {
      stop("To use custom database, please supply equidimensional character vectors of ligands and recepts")
    }
    message("Using custom database")
    ligands <- ligands
    recepts <- recepts
    lit.put <- data.frame(pair.name = paste(ligands,recepts,sep="_"), source_genesymbol = ligands, target_genesymbol = recepts)
  } else if((!is.null(ligands) | !is.null(recepts)) & database != "custom") {
    stop("To use custom ligand or receptor lists, set database = 'custom'")
  } else {
    lit.put <- scriabin::LoadLR(species = species, database = database)
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



#' Build summarized interaction graph weighted by predicted ligand activities
#'
#' @param object A Seurat object
#' @param nichenet_results List or matrix. Predicted ligand activities using Scriabin's implementation of NicheNet in RankActiveLigands
#' @param pearson.cutoff numeric. Threshold for determining which ligand activities are "active". Ligands below this threshold will be considered inactive and not used for weighting. Default: 0.075
#' @param scale.factors numeric. Determines the magnitude of ligand and receptor expression values weighting by their predicted activities. Given scale factors of c(x,y), ligand-receptor pairs where the ligand's pearson = pearson.cutoff will be weighted by a factor of x, and the ligand with the highest pearson will be weighted by a factor of y. Default: c(1.5,3).
#' @param assay Assay in Seurat object from which to pull expression values
#' @param slot Slot within assay from which to pull expression values
#' @param species Name of species for which to pull ligand-receptor interactions. One of "human", "mouse", or "rat"
#' @param database Name of ligand-receptor database to use. Default: "OmniPath"
#' When species is "human", one of: OmniPath, CellChatDB, CellPhoneDB, Ramilowski2015, Baccin2019, LRdb, Kirouac2010, ICELLNET, iTALK, EMBRACE, HPMR, Guide2Pharma, connectomeDB2020, talklr, CellTalkDB
#' When species is "mouse" or "rat", only "OmniPath" is supported.
#' To pass a custom ligand-receptor database to this function, set database = "custom"
#' @param ligands Character vector of custom ligands to use for interaction graph generation. Ignored unless database = "custom"
#' When ligands is supplied, recepts must also be supplied and equidimensional.
#' @param recepts Character vector of custom receptors to use for interaction graph generation. Ignored unless database = "custom"
#' When recepts is supplied, ligands must also be supplied and equidimensional.
#' @param ranked_genes Cell-resolved gene signatures, used only when specific = T
#' @param correct.depth Correct summarized interaction graph for sequencing depth by linear regression. The sequencing depth of a cell-cell pair is the sum of UMI counts for each cell
#' @param graph_name Name of summarized interaction graph to place into output. Default "weighted_interaction"
#'
#' @return Returns a Seurat object with a weighted summarized interaction graph in the Graphs slot
#' @import dplyr Seurat
#' @importFrom tibble rownames_to_column
#' @importFrom tibble column_to_rownames
#' @references Browaeys, et al. Nature Methods (2019)

#' @export
#'
#' @examples
BuildWeightedInteraction <- function (object, nichenet_results, assay = "SCT", slot = "data",
                                      pearson.cutoff = 0.075, scale.factors = c(1.5,3),
                                      species = "human", database = "OmniPath",
                                      ligands = NULL, recepts = NULL,
                                      correct.depth = T, graph_name = "weighted_interaction") {
  if(database=="custom") {
    if(is.null(ligands) | is.null(recepts)) {
      stop("To use custom database, please supply equidimensional character vectors of ligands and recepts")
    }
    message("Using custom database")
    ligands <- ligands
    recepts <- recepts
    lit.put <- data.frame(pair.name = paste(ligands,recepts,sep="_"), source_genesymbol = ligands, target_genesymbol = recepts)
  } else if((!is.null(ligands) | !is.null(recepts)) & database != "custom") {
    stop("To use custom ligand or receptor lists, set database = 'custom'")
  } else {
    lit.put <- scriabin::LoadLR(species = species, database = database)
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
  cell.exprs.rec <- merge(recepts.df[,2:3], cell.exprs,
                          by.x = "recepts", by.y = "gene", all.x = T)
  cell.exprs.rec <- cell.exprs.rec[order(cell.exprs.rec$id),
                                   ]
  cell.exprs.lig <- merge(ligands.df[,2:3], cell.exprs,
                          by.x = "ligands", by.y = "gene", all.x = T)
  cell.exprs.lig <- cell.exprs.lig[order(cell.exprs.lig$id),
                                   ]
  message("Weighting interaction matrix")

  if(is.null(nichenet_results)) {
    stop("NicheNet results must be included to perform SIG weighting")
  }
  message("Weighting interaction graph")
  if(class(nichenet_results) %notin% c("list","matrix")) {
    stop("NicheNet results must be supplied as one of list or matrix")
  }
  if(class(nichenet_results)=="list") {
    nichenet_results <- nichenet_results[names(nichenet_results) %in% colnames(object)]
    nichenet_results <- bind_rows(lapply(nichenet_results, FUN = function(x) {
      results <- x[[1]] %>% pull(pearson)
      names(results) <- x[[1]] %>% pull(test_ligand)
      return(results)
    }), .id = "cell")
    nnm <- reshape2::melt(nichenet_results) %>% dplyr::filter(value>pearson.cutoff)
    colnames(nnm) <- c("cell","ligand","pearson")
  }
  if(class(nichenet_results)=="matrix") {
    nichenet_results <- nichenet_results[,colnames(object)]
    nnm <- reshape2::melt(nichenet_results) %>% dplyr::filter(value>pearson.cutoff)
    colnames(nnm) <- c("ligand","cell","pearson")
  }
  if(nrow(nnm)==0) {
    warning("No active ligands to weight. Please check specified pearson cutoff")
    message("Using unweighted ligand-receptor matrices")
    a <- as.matrix(cell.exprs.lig[,3:ncol(cell.exprs.lig)])
    a[is.na(a)] <- 0
    b <- as.matrix(cell.exprs.rec[,3:ncol(cell.exprs.rec)])
    b[is.na(b)] <- 0
  }
  else {
    message(paste("Found active ligands to weight. Weighting with pearson cutoff ",pearson.cutoff))
    nnm$weight_factor <- scales::rescale(nnm$pearson, scale.factors)
    int.to.merge <- lit.put[,c(1,2,3)]
    colnames(int.to.merge) <- c("pair","ligand","recepts")
    nnm <- merge(nnm,int.to.merge,by = "ligand", all.x=T)
    nnm <- nnm[!is.na(nnm$recepts),]
    test2 <- reshape2::dcast(nnm, recepts~cell, value.var = "weight_factor", fun.aggregate = sum)
    test3 <- merge(cell.exprs.rec[,1:2],test2,all.x = T) %>% arrange(id)

    cell.exprs.rec.m <- cell.exprs.rec[,3:ncol(cell.exprs.rec)]
    weight.m <- as.matrix(test3[,3:ncol(test3)])
    weight.m[weight.m==0] <- 1
    weight.m[is.na(weight.m)] <- 1
    cells.bind <- colnames(cell.exprs.rec.m)[colnames(cell.exprs.rec.m) %notin% colnames(weight.m)]
    tobind <- matrix(1,nrow = nrow(weight.m),ncol = length(cells.bind))
    colnames(tobind) <- cells.bind
    weight.m <- cbind(weight.m,tobind)
    weight.m <- weight.m[,match(colnames(cell.exprs.rec.m),colnames(weight.m))]

    final <- weight.m * as.matrix(cell.exprs.rec.m)
    # weighted.rec <- cbind(cell.exprs.rec.m[,1:3],final)

    a <- as.matrix(cell.exprs.lig[,3:ncol(cell.exprs.lig)])
    a[is.na(a)] <- 0
    b <- as.matrix(final)
    b[is.na(b)] <- 0
  }
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



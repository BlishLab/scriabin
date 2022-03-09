#' Identify circuits between two timepoints
#'
#' @param seu A binned Seurat object with binning results stored in the "bins" column of the meta.data slot
#' @param nnr Ligand and target ranking results from NicheNet. Currently only the output of `PrioritizeLigands` is supported.
#' @param ranked_genes Single-cell gene signatures, output from `crGeneSig`
#' @param tps Character vector of two timepoints between which to find circuits. Time point identities must be in the correct sequential order.
#' @param assay.use Name of assay from which to gather gene expression data
#' @param slot.use Name of slot within assay from which to gather gene expression data
#' @param pearson.cutoff numeric. Threshold for determining which ligand activities are "active". Ligands below this threshold will be considered inactive and not used for weighting. Default: 0.075
#' @param split.by Meta.data column name indicating how full dataset should be split. This corresponds to the column that contain time point information
#'
#' @return Returns a data.frame containing all single-cell level circuits between two timepoints of a dataset
#' @import Seurat Matrix dplyr
#' @importFrom tidyr separate
#' @export
#'
#' @examples
FindCircuits <- function(seu, nnr, ranked_genes, tps, split.by = "orig.ident",
                         assay.use = "SCT", slot.use = "data", pearson.cutoff = 0.075) {
  ### things to fix in here: provide target linking strategy that doesn't rely on Prioritize Ligands.
  ### provide adequate LR resource support, currently relies only on Connectome's fantom5
  nnr <- nnr[colnames(seu)]
  ranked_genes <- ranked_genes[colnames(seu)]

  gsub_function = ".*[_]([^.]+)[.].*"
  orig_cell_names <- colnames(seu)

  sample_ids <- paste("project",seu@meta.data[,split.by], sep = "_")
  new_cell_names <- paste(sample_ids, 1:ncol(seu), sep = ".")

  seu <- RenameCells(seu, new.names = new_cell_names)
  names(nnr) <- names(ranked_genes) <- colnames(seu)

  nnr_l <- lapply(nnr, FUN = function(x) {x[[1]]}) #grabs ligands
  nnr_t <- lapply(nnr, FUN = function(x) {x[[2]]}) #grabs target links

  nnr_l <- nnr_l[!(lapply(nnr_l,length)==0)]
  nnr_t <- nnr_t[!(lapply(nnr_t,length)==0)]

  cell.exprs <- GetAssayData(seu, assay = assay.use, slot = slot.use)

  #first identify ligands in the receivers of the second timepoint

  nnr_l2 <- bind_rows(nnr_l[gsub(gsub_function,"\\1",names(nnr_l))==tps[2]], .id = "cell") %>%
    dplyr::filter(pearson>pearson.cutoff) #this is a list of significant ligands seen by receivers at t2

  ### now find cells at this second timepoint who have these ligands in their gene signature and express them

  gs_nnr <- ranked_genes[gsub(gsub_function,"\\1",names(nnr_l))==tps[2]]
  gs_nnr <- lapply(gs_nnr, function(x) {
    x <- names(x)
    x <- x[x %in% unique(nnr_l2$test_ligand)]
    return(x)
  })
  gs_nnr <- gs_nnr[lapply(gs_nnr,length)>0]
  gs_exprs <- cell.exprs[unique(nnr_l2$test_ligand),names(gs_nnr)]
  gs_exprs <- sapply(seq_along(1:length(gs_nnr)), function(x) {
    y <- gs_exprs[,names(gs_nnr)[x]]
    y[names(y) %notin% gs_nnr[[x]]] <- 0
    return(y)
  })
  colnames(gs_exprs) <- names(gs_nnr)

  gs_exprs <- gs_exprs[,Matrix::colMeans(gs_exprs)>0]
  ligands_search <- reshape2::melt(t(as.matrix(gs_exprs))) %>% dplyr::filter(value>0)
  colnames(ligands_search) <- c("cell","ligand","value")
  ligands_search$bins <- scriabin::mapvalues(ligands_search$cell, from = colnames(seu), to = seu$bins, warn_missing = F) #this is a list of significant ligands sent by senders at t2

  #now identify cells in the first timepoint with these targets
  nnr_l1 <- bind_rows(nnr_l[gsub(gsub_function,"\\1",names(nnr_l))==tps[1]], .id = "cell") %>%
    dplyr::filter(pearson>pearson.cutoff)
  nnr_l1$cell_ligand <- paste(nnr_l1$cell,nnr_l1$test_ligand,sep = "=") #this is a list of significant ligands seen by receivers at t1

  nnr_t1 <- bind_rows(nnr_t[gsub(gsub_function,"\\1",names(nnr_l))==tps[1]], .id = "cell") %>%
    dplyr::filter(target %in% unique(ligands_search$ligand))
  nnr_t1$cell_ligand <- paste(nnr_t1$cell,nnr_t1$ligand,sep = "=") #this is a list of significant targets upregulated by receivers at t1

  nnr_1 <- merge(nnr_l1[,c("cell_ligand","pearson")],
                 nnr_t1[,c("cell_ligand","target","weight")],
                 by = "cell_ligand") %>%
    tidyr::separate(cell_ligand, sep = "=", into = c("cell","ligand"))

  nnr_1$bins <- mapvalues(nnr_1$cell, from = colnames(seu), to = seu$bins, warn_missing = F)
  nnr_1 %<>% dplyr::filter(bins %in% unique(ligands_search$bins)) #this is a list of significant ligand-target links in the receivers at t1

  nnr_1$geneset <- unlist(lapply(seq_along(1:nrow(nnr_1)), function(x) {
    ifelse(nnr_1[x,"target"] %in% names(ranked_genes[[nnr_1[x,"cell"]]]),T,F)
  }))
  nnr_1 %<>% dplyr::filter(geneset==T) %>% dplyr::select(-geneset)

  ligands_search$gene_bin <- paste(ligands_search$ligand,ligands_search$bins, sep = "_")
  nnr_1$gene_bin <- paste(nnr_1$target,nnr_1$bins, sep = "_")

  links <- merge(nnr_1,ligands_search,by = "gene_bin", all.x = F, all.y = F) %>%
    dplyr::select(-one_of(c("gene_bin","bins.x"))) %>% dplyr::rename(bins = bins.y) #this is a list of significant ligand-target links in the receivers at t1 that are

  #links represents the first part of a circuit. To complete it, we find which cells in the second time points could be impacted.
  rec_cells <- as.character(unique(nnr_l2$cell))
  fantom5 <- Connectome::ncomms8866_human
  lit.put <- fantom5[fantom5$Pair.Evidence %in% c("literature supported",
                                                  "putative"), ]
  pairs <- lit.put[,c(2,4)] %>%
    mutate_all(as.character) %>%
    dplyr::filter(Ligand.ApprovedSymbol %in% unique(links$ligand.y))
  colnames(pairs) <- c("ligand","receptor")

  cell.exprs <- GetAssayData(seu, assay = assay.use, slot = slot.use)[rownames(seu) %in% c(unique(pairs$ligand),unique(pairs$receptor)),c(rec_cells)]
  # cell.exprs <- cell.exprs[rowSums(cell.exprs)>0,]

  nnr_filtered <- nnr_l2 %>% dplyr::filter(test_ligand %in% pairs$ligand)
  nnr_filtered$recept <- unlist(lapply(seq_along(1:nrow(nnr_filtered)), FUN = function(x) {
    recept <- pairs %>% dplyr::filter(ligand==as.character(nnr_filtered[x,"test_ligand"])) %>% pull(receptor)
    recept <- recept[recept %in% rownames(seu)]
    return(ifelse(sum(cell.exprs[recept,nnr_filtered$cell[x]])>0,T,F))
  }))
  nnr_filtered %<>% dplyr::filter(recept==T) %>% dplyr::rename(cell.z=cell) %>% dplyr::rename(ligand.y=test_ligand)

  links <- merge(links, nnr_filtered[,c("cell.z","ligand.y")], by = "ligand.y", all.x = F, all.y = F)
  links$tp1 <- tps[1]
  links$tp2 <- tps[2]
  links$cell.x <- scriabin::mapvalues(links$cell.x, from=new_cell_names, to = orig_cell_names, warn_missing = F)
  links$cell.y <- scriabin::mapvalues(links$cell.y, from=new_cell_names, to = orig_cell_names, warn_missing = F)
  links$cell.z <- scriabin::mapvalues(links$cell.z, from=new_cell_names, to = orig_cell_names, warn_missing = F)


  return(links[,c("tp1","ligand.x","pearson","cell.x",
                  "target","bins","cell.y","cell.z","tp2")])
}


#' Find all circuits in a multi-timepoint longitudinal dataset
#'
#' @param seu A binned Seurat object with binning results stored in the "bins" column of the meta.data slot
#' @param nnr Ligand and target ranking results from NicheNet. Currently only the output of `PrioritizeLigands` is supported.
#' @param ranked_genes Single-cell gene signatures, output from `crGeneSig`
#' @param all.tps Character vector of all timepoints between which to find circuits. Time point identities must be in the correct sequential order.
#' @param subsample logical. Subsample circuits for each pair of timepoints to the threshold set in subsample.n?
#' @param subsample.n numeric. If subsample=T, subsample circuits from each timepoint pair to this level
#'
#' @return
#' @export
#'
#' @examples
FindAllCircuits <- function(seu, nnr, ranked_genes = NULL, all.tps = NULL,
                            subsample = F, subsample.n = 1e5) {
  message("Identifying circuits . . . ")
  all_circuits <- bind_rows(pblapply(seq_along(1:(length(all.tps)-1)), function(x) {
    circuits <- FindCircuits(seu, nnr,
                             ranked_genes = ranked_genes,
                             tps = c(all.tps[x],all.tps[x+1]))
  }))
  message("Finished finding circuits.")

  if(subsample) {
    message(paste0("Subsampling to ",subsample.n," circuits"))
    all_circuits <- all_circuits[sample(1:nrow(all_circuits),subsample.n),]
  }
  return(all_circuits)
}

#' Stitch multi-timepoint circuits together
#'
#' @param all_circuits
#' @param all.tps
#'
#' @return
#' @export
#'
#' @examples
FormatCircuits <- function(all_circuits = NULL, all.tps = NULL) {
  test4 <- bind_rows(pblapply(seq_along(1:(length(all.tps)-1)), function(x){
    data_tp_pre <- all_circuits %>% dplyr::filter(tp1==all.tps[x])
    data_tp_pre$ligand_receiver.x <- paste(data_tp_pre$cell.x,
                                           data_tp_pre$ligand.x, sep = "=")
    data_tp_pre$ligand_receiver.z <- paste(data_tp_pre$cell.z,
                                           data_tp_pre$target, sep = "=")
    data_tp_pre <- data_tp_pre[,c("ligand_receiver.x","cell.y","ligand_receiver.z")]

    vec <- c()
    for(i in 1:(length(all.tps)-1)) {
      vec <- c(vec,
               paste("ligand-receiver",all.tps[i],sep = "_"),
               paste("sender",all.tps[i+1],sep = "_"),
               paste("ligand-receiver",all.tps[i+1],sep = "_"))
    }
    vec <- unique(vec)
    vec <- c(vec,"dataset")
    final_df <- data.frame(matrix(ncol=length(vec), nrow=nrow(data_tp_pre)))
    colnames(final_df) <- vec

    coln <- ((x*2)-1):((x*2)+1)
    final_df[,coln] <- data_tp_pre
    final_df$dataset <- x
    return(final_df)
  }))

  message("Matching circuit pairs . . .")
  for (x in 1:(length(all.tps)-2)) {
    message("Step1")
    var = colnames(test4)[(2*x)+1]
    pre <- as.data.frame(test4 %>% dplyr::filter(dataset==x))
    post <- as.data.frame(test4 %>% dplyr::filter(dataset==x+1))

    message("Step2")
    bad_post_indices <- test4 %>%
      dplyr::filter((!!sym(var)) %notin% pre[,var][pre[,var] %in% post[,var]]) %>%
      dplyr::filter(dataset==x+1)

    message("Step3")
    pre_indices <- test4[apply(test4,1,function(x) {sum(is.na(x))})<10,]

    message("Step4")
    good_indices <- test4 %>%
      dplyr::filter((!!sym(var)) %in% pre[,var][pre[,var] %in% post[,var]])

    message("Step5")
    good_indices %<>%
      group_by((!!sym(var))) %>%
      fill(everything(), .direction = 'updown') %>%
      ungroup %>%
      distinct(.keep_all = TRUE)

    good_indices$dataset <- x+1

    message("Step6")
    test4 %<>% dplyr::filter(!dataset %in% c(x,x+1))
    test4 <- rbind(pre_indices, good_indices, bad_post_indices, test4)
    message("Finished a loop of matching")
  }

  nnas <- ncol(test4)-apply(test4,1,function(x) {sum(is.na(x))})
  test4$n_stitched <- (nnas-2)/2
  test4$dataset <- NULL
  coln <- colnames(test4)[grepl("ligand-receiver",colnames(test4))]

  for (i in 1:sum(grepl("ligand-receiver",colnames(test4)))) {
    var <- coln[i]
    new_name1 <- paste0("receiver_",sub("^[^_]*_", "", coln[i]))
    new_name2 <- paste0("ligand_",sub("^[^_]*_", "", coln[i]))
    test4 %<>% tidyr::separate((!!sym(var)), sep = "=", into = c(new_name1,new_name2))
  }


  return(test4)
}

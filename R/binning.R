#' Align Datasets through binning
#'
#' @param seuObj A seurat object
#' @param split.by Name of meta.data column containing the factor by which the dataset will be split
#' @param dims Dimensions of reduction to use as input for neighbor graph calculation
#' @param snn.reduction Name of reduction used as input to building the SNN
#' @param anchor_score_threshold Anchor pairs scoring below this threshold will be discarded.
#' @param optim_quan.threshold Percentage of poor connectivity bins to remove on each iteration of connectivity optimization
#' @param optim_k.unique Mean number of datasets represented in each bin at which connectivity optimization will be considered complete
#'
#' @return A binned Seurat object with bin assignments in the "bins" column of the meta.data slot
#' @import dplyr Seurat
#' @export
#'
#' @examples
#' seu <- AlignDatasets(seu, split.by = "time")
AlignDatasets <- function(seuObj, split.by = "time.orig",
                          dims = 1:50, snn.reduction = "pca",
                          anchor_score_threshold = 0.5,
                          optim_quan.threshold = 0.1, optim_k.unique = 6,
                          verbose = F)
{
  gsub_function = ".*[=]([^.]+)[.].*"
  orig_cell_names <- colnames(seuObj)

  sample_ids <- paste("project",seuObj@meta.data[,split.by], sep = "=")
  new_cell_names <- paste(sample_ids, 1:ncol(seuObj), sep = ".")

  seuObj <- RenameCells(seuObj, new.names = new_cell_names)

  if(is.null(snn.reduction)) {
    message("\nNo reduction for SNN construction provided. Running PCA . . . ")
    seuObj <- RunPCA(seuObj, verbose=F)
  }

  message(paste0("Constructing neighbor graphs with reduction ",snn.reduction," . . . "))
  seuObj <- FindNeighbors(seuObj, dims = dims, verbose = F, reduction = snn.reduction)

  message("Splitting object . . . ")
  object.list <- SplitObject(seuObj, split.by = split.by) #split object

  message("Preparing integration . . . ")
  integration.features = SelectIntegrationFeatures(object.list) #prepare for integration

  object.ncells <- sapply(X = object.list, FUN = function(x) dim(x = x)[2])
  adims <- min(30, min(sapply(object.list, ncol)))
  if(adims<30) {
    warning("Fewer than 30 dimensions available for integration. Recommend choosing coarser cell type labels. Integration may fail.")
  }
  k.filter <- min(200, min(sapply(object.list, ncol)))

  message("Generating anchorsets . . . ")
  full.anchors <- FindIntegrationAnchors(object.list,k.filter = k.filter,dims = 1:(adims-1)) #Identify anchors

  #map anchors to cell names
  message("Mapping anchor names . . . ")
  full_anchors <- full.anchors@anchors

  cell_list <- bind_rows(lapply(seq_along(1:length(object.list)), function(x) {
    data.frame(id = paste(x,1:ncol(object.list[[x]]),sep = "_"), cell = colnames(object.list[[x]]))
  }))

  full_anchors$cell.name1 <- scriabin::mapvalues(paste(full_anchors$dataset1,full_anchors$cell1,sep = "_"),
                                       from = cell_list$id, to = cell_list$cell, warn_missing = F)

  full_anchors$cell.name2 <- scriabin::mapvalues(paste(full_anchors$dataset2,full_anchors$cell2,sep = "_"),
                                       from = cell_list$id, to = cell_list$cell, warn_missing = F)

  anchors_map <- full_anchors[full_anchors$score>anchor_score_threshold,] #remove anchors with poor scoring
  anchor_parts <- unique(c(anchors_map$cell.name1,anchors_map$cell.name2))

  #for each cell that participates in an anchor, what are all the other cells it anchors with?
  message("Generating anchor bins . . . ")
  n.results <- list()
  for (i in 1:length(anchor_parts)) {
    n.results[[i]] <- c(anchor_parts[i],anchors_map[anchors_map$cell.name1==anchor_parts[i],"cell.name2"])
    n.results[[i]] <- unique(c(n.results[[i]],c(anchors_map[anchors_map$cell.name1==anchor_parts[i],"cell.name2"])))
  }
  names(n.results) <- 1:length(n.results)

  # #remove anchor bins that do not have anchor cells in each timepoint
  # if(filter_anchors) {
  #   testresults <- lapply(n.results, FUN = function(x) {
  #     y <- gsub(gsub_function, "\\1", x)
  #     length(unique(y))
  #   })
  #   n.results <- n.results[testresults==length(object.list)]
  # }

  message("Generating connectivity matrix . . . ")
  SNN <- as.sparse(seuObj@graphs[[grep("_snn",Graphs(seuObj),value=T)]])
  outerresults <- sapply(n.results, function(v) Matrix::rowMeans(SNN[, v]))

  message("Optimizing number of bins . . . ")
  success = F
  n=0
  while(!success) {
    #assign cells to anchorset with highest connectivity
    ids <- apply(outerresults,1,FUN = function(i){
      m <- max(i, na.rm = T)
      mi <- which(x = i == m, arr.ind = TRUE)
      closest_cluster <- sample(x = names(x = i[mi]),
                                1)
      return(closest_cluster)
    })
    names(ids) <- colnames(seuObj)
    nbs_completion <- data.frame(cell=names(ids),id=ids,ident=gsub(gsub_function, "\\1", names(ids))) %>% dplyr::group_by(id) %>% dplyr::mutate(unique_types=n_distinct(ident))
    nbs_unique <- unique(nbs_completion[,c("id","unique_types")])
    success1 <- ifelse(mean(nbs_unique$unique_types)>optim_k.unique,T,F)
    success2 <- n>50
    success <- success1|success2
    if(verbose) {
      message(paste0("Average completion score: ",mean(nbs_unique$unique_types)))
    }
    #pull out the worst bins (perhaps those with the fewest unique_types and the lowest overall connectivity)
    quan.threshold = optim_quan.threshold
    quan.level <- round(unname(quantile(nbs_unique$unique_types,quan.threshold)))
    if(sum(nbs_unique$unique_types==quan.level)/nrow(nbs_unique)>quan.threshold) {
      nnbs_remove <- round(ncol(outerresults)*quan.threshold)
      bad_nbs <- outerresults[,nbs_unique %>% dplyr::filter(unique_types==quan.level) %>% pull(id)]
      nbs_remove <- data.frame(value=colSums(bad_nbs),names=colnames(bad_nbs)) %>% top_n(nnbs_remove,wt = value) %>% pull(names)
      outerresults <- outerresults[,!colnames(outerresults) %in% nbs_remove]
    }
    else {
      outerresults <- outerresults[,!colnames(outerresults) %in% (nbs_unique %>% dplyr::filter(unique_types==quan.level) %>% pull(id))]
    }
    n=n+1
    # message("Finished an iteration")
  }

  ###Check completion

  if(check_completion(nbs_completion)) {
    message(paste0("Finished! Found ",length(unique(ids))," bins"))
    seuObj$bins <- ids
    seuObj <- RenameCells(seuObj, new.names = orig_cell_names)
    return(seuObj)
  }


  ###Assign anything in a bad neighborhood to one that's more complete
  cells_reassign <- nbs_completion %>% dplyr::filter(unique_types<optim_k.unique) %>% dplyr::pull(cell)
  bad_nbs <- nbs_unique %>% dplyr::filter(unique_types<optim_k.unique) %>% dplyr::pull(id)
  outer_reassign <- outerresults[,colnames(outerresults) %in% nbs_unique$id]
  outer_reassign <- outer_reassign[cells_reassign,colnames(outer_reassign) %notin% bad_nbs]
  reassign_ids <- apply(outer_reassign,1,FUN = function(i){
    m <- max(i, na.rm = T)
    mi <- which(x = i == m, arr.ind = TRUE)
    closest_cluster <- sample(x = names(x = i[mi]),
                              1)
    return(closest_cluster)
  })
  names(reassign_ids) <- cells_reassign
  new_ids <- ifelse(names(ids) %in% names(reassign_ids),reassign_ids,ids)
  names(new_ids) <- names(ids)
  new_completion <- data.frame(cell=names(new_ids),id=new_ids,ident=gsub(gsub_function, "\\1", names(new_ids))) %>%
    dplyr::group_by(id) %>% dplyr::mutate(unique_types=n_distinct(ident))
  new_unique <- unique(new_completion[,c("id","unique_types")])
  if(mean(new_unique$unique_types)==length(object.list)) {
    seuObj$bins <- new_ids
    return(seuObj)
  }

  ###Search for cells that can be stolen by incomplete bins
  outer_steal <- outerresults[,colnames(outerresults) %in% new_unique$id]
  success = F
  tps <- unique(gsub(gsub_function, "\\1", names(new_ids)))
  steal_ids <- new_ids
  stolen_completion <- new_completion
  stolen_unique <- new_unique

  incomplete_ids <- sample(stolen_unique %>% dplyr::filter(unique_types<length(object.list)) %>% dplyr::pull(id))
  for (j in 1:length(incomplete_ids)) {
    ##for each randomly-sampled incomplete neighborhood
    noi <- incomplete_ids[j]
    ##what timepoint(s) is/are missing?
    toi <- tps[tps %notin% unique(gsub(gsub_function, "\\1", names(steal_ids)[steal_ids==noi]))]
    for (i in 1:length(toi)) {
      ##create a ranked-list of cells from that timepoint that could be stolen
      outer_steal_working <- outer_steal[grepl(toi[i],rownames(outer_steal)),noi]
      outer_steal_working <- outer_steal_working[outer_steal_working>0]
      outer_steal_working <- rev(outer_steal_working[order(outer_steal_working)])
      cells_to_steal <- names(outer_steal_working)
      ##in order, evaluate if that cell can be stolen without making that neighborhood incomplete
      sus_cells <- unlist(lapply(cells_to_steal,FUN = function(x) {
        defending_nb <- gsub(gsub_function, "\\1", names(steal_ids[steal_ids==steal_ids[x]]))
        ifelse(length(defending_nb[defending_nb==toi[i]])>1,"Available","Unavailable")
      }))
      ##steal the first cell that can be stolen
      if(length(sus_cells[sus_cells=="Available"])>1) {
        cell_steal <- names(outer_steal_working)[sus_cells=="Available"][1]
        steal_ids[cell_steal] <- noi
        if(verbose) {message("Stole a cell!")}
      }
      else {
        if(verbose){message(paste0("Neighborhood ",noi," has no potential targets"))}
      }
    }
    ##reevaluate completion
    stolen_completion <- data.frame(cell=names(steal_ids),id=steal_ids,ident=gsub(gsub_function, "\\1", names(steal_ids))) %>%
      dplyr::group_by(id) %>% dplyr::mutate(unique_types=n_distinct(ident))
    stolen_unique <- unique(stolen_completion[,c("id","unique_types")])
  }
  incomplete_ids <- stolen_unique %>% dplyr::filter(unique_types<length(object.list)) %>% dplyr::pull(id)
  if(verbose){message(paste0("Merging partners must be found for ",length(incomplete_ids)," bins"))}


  message("Merging bins")
  merge_ids <- steal_ids
  merge_completion <- stolen_completion
  merge_unique <- stolen_unique
  ##identify incomplete bins
  incomplete_ids <- sample(incomplete_ids)
  for (i in 1:length(incomplete_ids)) {
    noi <- incomplete_ids[i]
    toi <- tps[tps %notin% unique(gsub(gsub_function, "\\1", names(merge_ids)[merge_ids==noi]))]
    outerresults_merge <- outerresults[,colnames(outerresults) %in% merge_ids]
    outerresults_cor <- cor(outerresults_merge)
    outerresults_noi <- outerresults_cor[,noi]
    outerresults_noi <- rev(outerresults_noi[order(outerresults_noi)])[2:length(outerresults_noi)]
    nb_merge <- names(outerresults_noi)
    merging_partners <- unlist(lapply(nb_merge,FUN = function(x) {
      partner_nb <- gsub(gsub_function, "\\1", names(merge_ids[merge_ids==x]))
      ifelse(sum(toi %in% partner_nb)==length(toi),"Good","Bad")
    }))
    if(length(merging_partners[merging_partners=="Good"])>1) {
      nb_to_merge <- names(outerresults_noi)[merging_partners=="Good"][1]
      merge_ids[merge_ids==noi] <- nb_to_merge
      message("Merged a neighborhood!")
      if(max(outerresults_noi)<0) {
        if(verbose){message(paste0("Warning! Neighborhood ",noi," is anti-correlated with all bins"))}
      }
    }
    else {
      message(paste0("Neighborhood ",noi," has no potential merging partners"))
    }

  }
  merge_completion <- data.frame(cell=names(merge_ids),id=merge_ids,ident=gsub(gsub_function, "\\1", names(merge_ids))) %>%
    dplyr::group_by(id) %>% dplyr::mutate(unique_types=n_distinct(ident))
  merge_unique <- unique(merge_completion[,c("id","unique_types")])

  message(paste0("Finished! Found ",length(unique(merge_ids))," bins"))

  seuObj$bins <- merge_ids

  seuObj <- RenameCells(seuObj, new.names = orig_cell_names)

  return(seuObj)
}


#' Test bin connectivity
#'
#' @param seu_oi A Seurat object with bin identities in the "bins" column of the meta.data slot
#' @param SNN The SNN to use for connectivity testing
#' @param bin Character string of the bin ID to test
#' @param split.by Name of meta.data column containing the factor by which the dataset will be split
#' @param sigtest_cell_types Name of meta.data column containing the cell type labels to use for connectivity testing
#'
#' @return p-value of bin connectivity significance
#' @export
#'
#' @examples
#' random_connectivity_test()
random_connectivity_test <- function(seu_oi = seu_oi,
                                     SNN=SNN,bin=bin,split.by=split.by,
                                     sigtest_cell_types=sigtest_cell_types) {
  true_bin <- colnames(seu_oi)[seu_oi$bins==bin]
  bin_comp <- as.data.frame(table(seu_oi@meta.data[seu_oi$bins==bin,split.by],
                                  seu_oi@meta.data[seu_oi$bins==bin,sigtest_cell_types])) %>%
    dplyr::filter(Freq>0)
  random_bin <- unlist(lapply(seq_along(1:nrow(bin_comp)), FUN = function(y) {
    sample(colnames(seu_oi)[seu_oi@meta.data[,split.by]==bin_comp[y,"Var1"] & seu_oi@meta.data[,sigtest_cell_types]==bin_comp[y,"Var2"]],bin_comp[y,"Freq"])
  }))
  a <- as.vector(SNN[true_bin,true_bin])
  b <- as.vector(SNN[random_bin,random_bin])
  return(wilcox.test(a,b,alternative = "greater",exact=F)$p.value)
}

#' Align Datasets through binning
#'
#' @param seuObj A seurat object
#' @param split.by Name of meta.data column containing the factor by which the dataset will be split
#' @param dims Dimensions of reduction to use as input for neighbor graph calculation
#' @param coarse_cell_types Name of meta.data column containing cell type calls for binning. Each set of cell type calls will be binned separately.
#' @param sigtest_cell_types Name of meta.data column containing the cell type labels to use for connectivity testing
#' @param snn.reduction Name of reduction used as input to building the SNN
#' @param anchor_score_threshold Anchor pairs scoring below this threshold will be discarded.
#' @param optim_quan.threshold Percentage of poor connectivity bins to remove on each iteration of connectivity optimization
#' @param optim_k.unique Mean number of datasets represented in each bin at which connectivity optimization will be considered complete
#'
#' @return A binned Seurat object with bin assignments in the "bins" column of the meta.data slot
#' @import dplyr Seurat pbapply
#' @export
#'
#' @examples
#' seu <- BinDatasets(seu, split.by = "time")
BinDatasets <- function(seu, split.by = "time.orig", dims = 1:50,
                        coarse_cell_types = NULL, sigtest_cell_types = NULL,
                        snn.reduction = "pca", anchor_score_threshold = 0.5,
                        optim_quan.threshold = 0.1, optim_k.unique = NULL, verbose = F)
{
  if(is.null(coarse_cell_types)) {
    warning("It is recommend to specify coarse cell types to improve bin significance testing. \nWhen not specified, Scriabin defaults to generating dataset-wide bins and testing significance based on cluster results.")
    # status <- readline(prompt = "Do you wish to continue? (yes/no): ")
    # if(status=="no"){
    #   stop("Terminated by user.")
    # }
  }
  if(!is.null(sigtest_cell_types) & sigtest_cell_types %notin% colnames(seu@meta.data)) {
    stop("sigtest_cell_types supplied but not found in meta.data slot of object")
  }
  if(is.null(optim_k.unique)) {
    n <- length(unique(seu@meta.data[,split.by]))
    optim_k.unique = n*2/3
    message("Setting optim_k.unique to ",format(round(optim_k.unique, 2), nsmall = 2))
  }

  if(!is.null(coarse_cell_types)) {
    message("Binning with coarse cell types")
    if(coarse_cell_types %notin% colnames(seu@meta.data)) {
      stop("coarse_cell_types not present in meta.data slot of object")
    }
    #check coarse cell IDs are present in all samples to be aligned
    ct_check <- as.data.frame(table(seu@meta.data[,split.by],seu@meta.data[,coarse_cell_types]))
    if(sum(ct_check$Freq==0)>0) {
      stop("Coarse cell types must be present in all samples to be binned")
    }
    seu_ct_split <- SplitObject(seu, split.by = coarse_cell_types)
    bin_ids <- lapply(seq_along(1:length(seu_ct_split)), function(x) {
      seu_oi <- AlignDatasets(seuObj = seu,
                              dims = dims,
                              anchor_score_threshold = anchor_score_threshold,
                              split.by = split.by,
                              optim_quan.threshold = optim_quan.threshold,
                              optim_k.unique = optim_k.unique,
                              snn.reduction = snn.reduction, verbose = verbose)
      message("Testing bin significance")
      SNN <- as.sparse(seu_oi@graphs[[grep("_snn",Graphs(seu_oi),value=T)]])

      if(is.null(sigtest_cell_types)) {
        seu_oi <- FindClusters(seu_oi)
        sigtest_cell_types <- "seurat_clusters"
        message("Using clusters to test bin significance")
      }

      bin_p <- unlist(pblapply(seq_along(1:length(unique(seu_oi$bins))), FUN = function(j) {
        bin = unique(seu_oi$bins)[j]
        random_distribution <- replicate(100, random_connectivity_test(seu_oi = seu_oi, SNN=SNN,
                                                                       bin=bin,
                                                                       split.by = split.by,
                                                                       sigtest_cell_types = sigtest_cell_types))
        return(sum(random_distribution>0.05)<5)
      }))
      browser()
      names(bin_p) <- unique(seu_oi$bins)

      anno.overlap <- t(reshape2::dcast(as.data.frame(table(as.character(seu_oi$bins),
                                                            seu_oi@meta.data[,sigtest_cell_types])),
                                        formula = Var1~Var2,value.var = "Freq") %>% column_to_rownames(var = "Var1"))
      anno.overlap <- t(100*anno.overlap/rowSums(anno.overlap))
      bin_max <- data.frame(max=apply(anno.overlap,1,max)) %>% rownames_to_column("bin")
      exempt_bins <- bin_max %>% dplyr::filter(max>95) %>% dplyr::pull(bin)

      nonsig_bins <- names(bin_p)[!bin_p]
      nonsig_bins <- nonsig_bins[nonsig_bins %notin% exempt_bins]

      message("Found ", length(nonsig_bins), " non-significant bin(s)")
      if(length(nonsig_bins)>0) {
        message("Merging non-significant bins")
        for(i in 1:length(nonsig_bins)) {
          bins = unique(seu_oi$bins)
          bin = nonsig_bins[i]
          bin_cells = colnames(seu_oi)[seu_oi$bins==bin]
          others = bins[bins %notin% nonsig_bins]
          con <- sapply(seq_along(1:length(others)), function(x) {
            mean(SNN[bin_cells,colnames(seu_oi)[seu_oi$bins==others[x]]])
          })
          seu_oi$bins[seu_oi$bins==bin] <- others[scriabin::resample(which(con==max(con)),1)]

        }
      }

      seu_oi$bins <- paste(unique(seu@meta.data[,coarse_cell_types])[x],seu_oi$bins,sep = "=")
      return(seu_oi)

    })

    bin_id_list <- unlist(lapply(bin_ids, function(x) {x$bins}))
    bin_rank <- as.character(rank(unique(bin_id_list)))
    names(bin_rank) <- unique(bin_id_list)
    bin_id_list <- scriabin::mapvalues(bin_id_list, from = names(bin_rank), to = bin_rank)
    seu$bins <- scriabin::mapvalues(colnames(seu), from = names(bin_id_list), to = bin_id_list)
    return(seu)
  }

  #without coarse cell types
  else {
    seu_oi <- AlignDatasets(seuObj = seu,
                            dims = dims,
                            anchor_score_threshold = anchor_score_threshold,
                            split.by = split.by,
                            optim_quan.threshold = optim_quan.threshold,
                            optim_k.unique = optim_k.unique,
                            snn.reduction = snn.reduction,
                            verbose = verbose)

    message("Testing bin significance")
    SNN <- as.sparse(seu_oi@graphs[[grep("_snn",Graphs(seu_oi),value=T)]])

    if(is.null(sigtest_cell_types)) {
      seu_oi <- FindClusters(seu_oi)
      sigtest_cell_types <- "seurat_clusters"
      message("Using clusters to test bin significance")
    }

    bin_p <- unlist(pblapply(seq_along(1:length(unique(seu_oi$bins))), FUN = function(j) {
      bin = unique(seu_oi$bins)[j]
      random_distribution <- replicate(100, random_connectivity_test(seu_oi = seu_oi, SNN=SNN,
                                                                     bin=bin,
                                                                     split.by = split.by,
                                                                     sigtest_cell_types = sigtest_cell_types))
      return(sum(random_distribution>0.05)<5)
    }))
    names(bin_p) <- unique(seu_oi$bins)

    anno.overlap <- t(reshape2::dcast(as.data.frame(table(as.character(seu_oi$bins),
                                                          seu_oi@meta.data[,sigtest_cell_types])),
                                      formula = Var1~Var2,value.var = "Freq") %>% column_to_rownames(var = "Var1"))
    anno.overlap <- t(100*anno.overlap/rowSums(anno.overlap))
    bin_max <- data.frame(max=apply(anno.overlap,1,max)) %>% rownames_to_column("bin")
    exempt_bins <- bin_max %>% dplyr::filter(max>95) %>% dplyr::pull(bin)

    nonsig_bins <- names(bin_p)[!bin_p]
    nonsig_bins <- nonsig_bins[nonsig_bins %notin% exempt_bins]

    message("Found ", length(nonsig_bins), " non-significant bin(s)")
    if(length(nonsig_bins)>0) {
      message("Merging non-significant bins")
      for(i in 1:length(nonsig_bins)) {
        bins = unique(seu_oi$bins)
        bin = nonsig_bins[i]
        bin_cells = colnames(seu_oi)[seu_oi$bins==bin]
        others = bins[bins %notin% nonsig_bins]
        con <- sapply(seq_along(1:length(others)), function(x) {
          mean(SNN[bin_cells,colnames(seu_oi)[seu_oi$bins==others[x]]])
        })

        seu_oi$bins[seu_oi$bins==bin] <- others[scriabin::resample(which(con==max(con)),1)]

      }
    }

    return(seu_oi)

  }

}

#' Helper function to check completion of bins
#'
#' @param nbs
#'
#' @return logical indicating if binning is complete
#' @export
#'
#' @examples
check_completion <- function(nbs) {
  unique(nbs$unique_types)==length(unique(nbs$ident))
}

#' Heatmap of bin-annotation overlap
#'
#' @param seu A seurat object with binning IDs in the "bins" column of meta.data
#' @param cell.type.calls Name of meta.data column containing cell type calls on which to visualize overlap
#'
#' @return A heatmap depicting cell type calls within each bin
#' @import ComplexHeatmap
#' @importFrom tibble rownames_to_column
#' @importFrom tibble column_to_rownames
#' @importFrom circlize colorRamp2
#' @export
#'
#' @examples
#' \dontrun{
#' BinAnnotationPlot(seu)
#' }
BinAnnotationPlot <- function(seu, cell.type.calls = "celltype.l2") {
  anno.overlap <- reshape2::dcast(as.data.frame(table(as.character(seu$bins),seu@meta.data[,cell.type.calls])),
                                  formula = Var1~Var2,value.var = "Freq") %>% column_to_rownames(var = "Var1")
  anno.overlap <- 100*anno.overlap/rowSums(anno.overlap)

  rowanno <- as.data.frame(table(as.character(seu$bins)))
  colnames(rowanno) <- c("bin","freq")
  bin_max <- data.frame(max=apply(anno.overlap,1,max)) %>% rownames_to_column("bin")
  rowanno <- merge(rowanno,bin_max,by = "bin") %>% column_to_rownames("bin")
  rowanno <- rowanno[match(rownames(rowanno),rownames(anno.overlap)),]

  la = rowAnnotation("Bin size" = rowanno$freq, "Max Overlap" = rowanno$max)

  Heatmap(as.matrix(anno.overlap), name = "% within single\ncell type", cluster_rows = T, cluster_columns = T,
          col = colorRamp2(c(0,50,100),c("lightsteelblue","yellow","red")), show_row_names = F,
          left_annotation = la)
}


#' Plot composition of each bin
#'
#' @param seu A seurat object with binning IDs in the "bins" column of meta.data
#' @param split.by Name of meta.data column defining how sub-datasets were binned
#' @param fill.colors Character of colors for each sub-dataset ID
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' BinCompositionAnalysis(seu)
#' }
BinCompositionAnalysis <- function(seu, split.by = "time.orig", fill.colors = NULL) {
  tp_props <- as.data.frame(table(seu$bins, seu@meta.data[,split.by]))
  tp_counts <- as.data.frame(table(seu@meta.data[,split.by]))
  tp_props$total <- as.numeric(as.character(scriabin::mapvalues(tp_props$Var2, from = tp_counts$Var1,
                                                      to = tp_counts$Freq, warn_missing = F)))
  tp_props$Freq <- 100*tp_props$Freq/tp_props$total
  colnames(tp_props) <- c("bin","tp","Freq","total")
  bin_mean <- aggregate(tp_props$Freq,by = list(tp_props$bin), FUN=mean)
  tp_props$mean <- as.numeric(as.character(scriabin::mapvalues(tp_props$bin, from = bin_mean$Group.1, to = bin_mean$x, warn_missing = F)))
  tp_props$tp_dev <- abs(tp_props$Freq-tp_props$mean)
  tp_mat <- reshape2::dcast(tp_props[,c("bin","tp","Freq")], formula = bin~tp, value.var = "Freq")

  data_plot <- data.frame(bin = tp_mat$bin)

  bin_sd <- aggregate(tp_props$Freq,by = list(tp_props$bin), FUN=sd)
  data_plot$dev <- as.numeric(as.character(scriabin::mapvalues(data_plot$bin, from = bin_sd$Group.1, to = bin_sd$x, warn_missing = F)))

  bin_size <- as.data.frame(table(seu$bins))
  data_plot$size <- as.numeric(as.character(scriabin::mapvalues(data_plot$bin, from = bin_size$Var1, to = bin_size$Freq, warn_missing = F)))

  data_plot$region <- 1:nrow(data_plot)

  data_plot <- cbind(data_plot,tp_mat[,2:ncol(tp_mat)])

  data_plot$dev <- data_plot$dev/data_plot$size
  data_plot$dev_scale <- scales::rescale(data_plot$dev, to = c(0,max(data_plot$size)))
  data_plot$dev <- data_plot$dev_scale/10

  p <- ggplot() + geom_scatterpie(aes(x = size, y = dev_scale, group = region, r = dev),
                             data = data_plot, cols = colnames(tp_mat)[-1]) +
    coord_equal() + theme_cowplot() +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
    labs(x = "Bin size", y = "Total bin deviance\nfrom expected composition", fill = "Dataset")
  if(!is.null(fill.colors)) {
    p <- p + scale_fill_manual(values = fill.colors)
  }
  else {
    p <- p + scale_fill_igv()
  }
}



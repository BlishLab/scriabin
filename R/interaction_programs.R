
#' Discover co-expressed ligand-receptor interaction programs
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
#' @param specific logical. When TRUE, consider only the genes in each cell's predefined gene signature (see crGeneSig) as expressed. Default FALSE
#' @param ranked_genes Cell-resolved gene signatures, used only when specific = T
#' @param return.mat logical. Returns ligand-receptor covariance matrix and TOM along with modules and intramodular connectivity. Required for significance testing.
#' @param softPower softPower threshold for adjacency matrix. If unspecified, will be chosen automatically based on the minimum softPower that results in a scale-free topology fitting index (R^2) of greater than r2_cutoff (by default 0.6)
#' @param min.size Minimum size of each interaction program
#' @param plot.mods Plot modules and associated dendrograms during analysis
#' @param tree.cut.quantile The dendrogram tree height quantile at which the dendrogram should be cut. Higher values lead to fewer, smaller modules.
#' @param iterate.threshold When a dataset has more cells than `iterate.threshold`, the TOM will be approximated iteratively. For each iteration, a randomly subsampled dataset of size `iterate.threshold` will be used to construct the CCIM and generate the TOM.
#' @param n.iterate For datasets larger than `iterate.threshold`, determines how many iterations to perform to approximate the TOM.
#' @param threads To enable WGCNA multi-threading, specify the number of threads to allow. This parameter may be required when running on multiple cores
#' @param r2_cutoff The softPower will be chosen as the minimum value that satisfies a scale-free topology fitting index (R^2) of r2_cutoff (by default: 0.6)
#' @param cell_types `meta.data` column name corresponding to cell types or clusters. If iteratively approximating the TOM, when specified the sequences of subsampled CCIM will be sampled proportionally to annotations present in this grouping.
#' @param min.cell When `cell_types` is specified, the minimum number of cells in a cell type that will be included when generating each iterative CCIM (default: 3)
#'
#' @return When return.mat = T, returns a list of length 4 containing ligand-receptor covariance matrix, TOM, module lists, and intramodular connectivity. Otherwise, returns a list of length 2 containing only module lists and intramodular connectivity.
#' @import qlcMatrix WGCNA flashClust dynamicTreeCut reshape2 pbapply
#' @export
#'
#' @examples
InteractionPrograms <- function(object, assay = "SCT", slot = "data",
                                species = "human", database = "OmniPath",
                                ligands = NULL, recepts = NULL,
                                iterate.threshold = 300, n.iterate = NULL,
                                specific = F, ranked_genes = NULL,
                                return.mat = T, softPower = NULL,
                                r2_cutoff = 0.6,
                                min.size = 5, plot.mods = F,
                                tree.cut.quantile = 0.4,
                                threads = NULL, cell_types = NULL,
                                min.cell = 3) {
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
    cell.exprs <- GetAssayData(object, assay = assay, slot = slot)[genes.use,]
  }

  ligands.df <- data.frame(ligands)
  ligands.df$id <- 1:nrow(ligands.df)
  recepts.df <- data.frame(recepts)
  recepts.df$id <- 1:nrow(recepts.df)


  if(ncol(object)>iterate.threshold) {
    message("\nIteratively generating interaction matrix")
    if(is.null(n.iterate)) {
      n.rep <- round(sqrt(ncol(object)/iterate.threshold))
    }
    else {
      n.rep = n.iterate
    }
    message(paste("Will perform",n.rep,"iterations to approximate TOM"))
    if (is.null(cell_types)) {
      warning("We recommend setting a cell_types parameter so that all cell types are included in each sequence of TOM generation")
    }
    mat_list <- lapply(seq_along(1:n.rep), function(z) {

      ## (1) Subsample proportional to cell type ##
      if(!is.null(cell_types)) {
        sub_prop <- (object@meta.data %>% rownames_to_column("cell"))[,c("cell",cell_types)]
        colnames(sub_prop) <- c("cell","var")
        cells <- sub_prop %>% group_by(var) %>% dplyr::mutate(prop = round(iterate.threshold*n()/nrow(.))) %>%
          dplyr::mutate(prop = ifelse(prop<min.cell,min.cell,prop)) %>% group_by(var) %>%
          dplyr::sample_n(prop) %>% pull(cell)
        cell.exprs.sub <- as.data.frame(cell.exprs[,cells]) %>% rownames_to_column(var = "gene")
      } else {
        cell.exprs.sub <- as.data.frame(cell.exprs[,sample(colnames(cell.exprs),iterate.threshold)]) %>% rownames_to_column(var = "gene")
      }

      ## (2) Define which ligands and receptors I'll accept on the other side ##
      cell.exprs.rec <- merge(recepts.df, cell.exprs.sub,
                              by.x = "recepts", by.y = "gene", all.x = T)
      cell.exprs.rec <- cell.exprs.rec[order(cell.exprs.rec$id),
                                       ]
      cell.exprs.lig <- merge(ligands.df, cell.exprs.sub,
                              by.x = "ligands", by.y = "gene", all.x = T)
      cell.exprs.lig <- cell.exprs.lig[order(cell.exprs.lig$id),
                                       ]

      a <- as.matrix(cell.exprs.lig[,3:ncol(cell.exprs.lig)])
      a[is.na(a)] <- 0
      b <- as.matrix(cell.exprs.rec[,3:ncol(cell.exprs.rec)])
      b[is.na(b)] <- 0

      m <- sqrt(as.sparse((pbsapply(1:nrow(a), function(i) tcrossprod(a[i, ], b[i, ])))))
      colnames(m) <- paste(cell.exprs.lig$ligands, cell.exprs.rec$recepts, sep = "=")
      cna <- rep(colnames(a),ncol(a))
      cnb <- rep(colnames(a),each=ncol(a))
      rownames(m) <- paste(cna,cnb,sep = "=")
      m <- m[,Matrix::colSums(m)>0]
      m_cor <- 0.5+(0.5*corSparse(m))
      rownames(m_cor) <- colnames(m)
      colnames(m_cor) <- colnames(m)
      m_cor
      # m
    })

    lr_union <- unique(unlist(lapply(mat_list,rownames)))
    mat_list_2 <- pblapply(mat_list, function(m) {
      missing <- lr_union[lr_union %notin% rownames(m)]
      if(length(missing)==0) {
        return(m[lr_union,lr_union])
      } else if(length(missing)==1) {
        m <- rbind(m,missing=NA)
        rownames(m)[nrow(m)] <- missing
        m <- cbind(m,missing=NA)
        colnames(m)[ncol(m)] <- missing
        return(m[lr_union,lr_union])
      } else if(length(missing)>1) {
        add_rows <- matrix(NA, nrow = length(missing), ncol = ncol(m))
        rownames(add_rows) <- missing
        m <- rbind(m,add_rows)

        add_cols <- matrix(NA, nrow = nrow(m), ncol = length(missing))
        colnames(add_cols) <- missing
        m <- cbind(m,add_cols)

        return(m[lr_union,lr_union])
      }
    })
    m_cor <- apply(simplify2array(mat_list_2), 1:2, function(x) {mean(x,na.rm = T)})
    # mat_list <- lapply(m_list, function(m) {
    #
    # })
    #filter to intersection
    ## (3) Union NOT intersection ##
    ## Fill missing rows with NA ##
    # cor_names <- Reduce(intersect,lapply(mat_list,colnames))
    # # m <- do.call(rbind,lapply(m_list, function(x) {x[,cor_names]}))
    # # m <- m[sample(1:nrow(m),iterate.threshold^2),]
    #
    # mat_list <- lapply(mat_list, function(x) {x[cor_names,cor_names]})
    # m_cor <- apply(simplify2array(mat_list), 1:2, median)
  }
  else {
    message(paste("\nGenerating Interaction Matrix..."))
    cell.exprs.sub <- as.data.frame(cell.exprs) %>% rownames_to_column(var = "gene")
    cell.exprs.rec <- merge(recepts.df, cell.exprs.sub,
                            by.x = "recepts", by.y = "gene", all.x = T)
    cell.exprs.rec <- cell.exprs.rec[order(cell.exprs.rec$id),
    ]
    cell.exprs.lig <- merge(ligands.df, cell.exprs.sub,
                            by.x = "ligands", by.y = "gene", all.x = T)
    cell.exprs.lig <- cell.exprs.lig[order(cell.exprs.lig$id),
    ]

    a <- as.matrix(cell.exprs.lig[,3:ncol(cell.exprs.lig)])
    a[is.na(a)] <- 0
    b <- as.matrix(cell.exprs.rec[,3:ncol(cell.exprs.rec)])
    b[is.na(b)] <- 0

    m <- sqrt(as.sparse((pbsapply(1:nrow(a), function(i) tcrossprod(a[i, ], b[i, ])))))
    colnames(m) <- paste(cell.exprs.lig$ligands, cell.exprs.rec$recepts, sep = "=")
    cna <- rep(colnames(object),ncol(object))
    cnb <- rep(colnames(object),each=ncol(object))
    rownames(m) <- paste(cna,cnb,sep = "=")
    m <- m[,Matrix::colSums(m)>0]
    m_cor <- 0.5+(0.5*corSparse(m))
    rownames(m_cor) <- colnames(m)
    colnames(m_cor) <- colnames(m)
  }
  if(!is.null(threads)) {
    enableWGCNAThreads(nThreads = threads)
  }
  if(is.null(softPower)) {
    message("Automatically selecting softPower . . .")
    sp_det <- pickSoftThreshold.fromSimilarity(m_cor, RsquaredCut = r2_cutoff)
    softPower = sp_det$powerEstimate
    if(is.na(softPower) | softPower>3) {
      warning("No appropriate softPower found to reach minimum scale free topology fit. Proceeding without soft thresholding, interpret results with caution")
      softPower = 1
    }
  }

  message("Identifying modules")
  adj <- m_cor^softPower
  tom <- TOMsimilarity(adj, TOMType = "signed")
  colnames(tom) <- rownames(tom) <- rownames(m_cor)
  geneTree = flashClust(as.dist(1-tom), method = "complete")
  dynamicMods = cutreeDynamic(dendro = geneTree,
                              method="tree", minClusterSize = min.size, cutHeight = quantile(geneTree$height, 0.5));
  dynamicColors = labels2colors(dynamicMods)
  if(plot.mods) {
    plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
  }
  module_colors= setdiff(unique(dynamicColors), "grey")
  modules = lapply(module_colors, function(x){colnames(m_cor)[which(dynamicColors==x)]})
  names(modules) = module_colors

  #Calculate module membership and intramodule connectivity

  #create color list
  module_melt <- reshape2::melt(modules)
  colors <- scriabin::mapvalues(colnames(m_cor),
                                from = module_melt$value,
                                to = module_melt$L1,
                                warn_missing = F)
  colors[colors %notin% unique(module_melt$L1)] <- "grey"

  Alldegrees1=intramodularConnectivity(adj, colors)
  names(modules) <- paste0("IP-",1:length(modules))

  if(return.mat) {
    return(list(cor_mat = m_cor, tom = tom, modules = modules, connectivity = Alldegrees1))
  }
  else {
    return(list(modules = modules, connectivity = im_results))
  }
}



#' Find all interaction programs in a multi-sample dataset
#'
#' @param seu A Seurat object
#' @param group.by Meta.data column name defining samples into which object should be split. Interaction programs will be found for each unique value in this column.
#' @param sim_threshold During module merging from different samples, the Jaccard overlap index threshold above which modules will be merged. Default: 0.15
#' @param ... Additional arguments passed to `InteractionPrograms`
#'
#' @return Returns a list of length 4 with ligand-receptor covariance matrices, TOMs, modules, and intramodular connectivity for each sample.
#' @import Seurat ade4 dplyr
#' @export
#'
#' @examples
FindAllInteractionPrograms <- function(seu, group.by = NULL, sim_threshold = 0.15, ...) {
  if(is.null(group.by)) {
    message("Grouping by current idents")
    seu$mod_grouping <- Idents(seu)
    group.by = "mod_grouping"
  }

  seu_split <- SplitObject(seu, split.by = group.by)
  q_mods <- lapply(seu_split, function(x) {
    InteractionPrograms(object = x, return.mat = T, ...)
  })

  mod_list <- unlist(lapply(q_mods,function(x){x[[3]]}), recursive = F)
  tom_list <- lapply(q_mods, function(x){x[[2]]})
  m_cor_list <- lapply(q_mods,function(x){x[[1]]})
  con_list <- lapply(q_mods,function(x){x[[4]]})
  names(con_list) <- names(tom_list) <- names(m_cor_list) <- names(q_mods)

  #Merge similar modules
  message("Merging similar modules")

  lrsum <- t(do.call(rbind,lapply(1:length(names(mod_list)), function(x) {
    unique(unlist(mod_list)) %in% mod_list[[x]]
  })))
  lrsum[lrsum==T] <- 1
  colnames(lrsum) <- names(mod_list)
  rownames(lrsum) <- unique(unlist(mod_list))

  d <- as.matrix(dist.binary(t(lrsum), method = 1))
  d <- 1-d^2

  merge_num = 1
  while(max(d[d<1])>sim_threshold) {
    tomerge <- rownames(which(d==max(d[d<1]), arr.ind = T))
    new_mod <- unique(unlist(mod_list[tomerge]))
    mod_list <- mod_list[names(mod_list) %notin% tomerge]
    mod_list[[length(mod_list)+1]] <- new_mod
    names(mod_list)[length(mod_list)] <- paste0("merge_",merge_num)
    merge_num = merge_num+1

    lrsum <- t(do.call(rbind,lapply(1:length(names(mod_list)), function(x) {
      unique(unlist(mod_list)) %in% mod_list[[x]]
    })))
    lrsum[lrsum==T] <- 1
    colnames(lrsum) <- names(mod_list)
    rownames(lrsum) <- unique(unlist(mod_list))

    d <- as.matrix(dist.binary(t(lrsum), method = 1))
    d <- 1-d^2
  }
  names(mod_list) <- paste0("IP-",1:length(mod_list))
  return(list(cor_mat = m_cor_list, tom = tom_list, modules = mod_list, connectivity = con_list))
}


#' Test Interaction program statistical significance
#'
#' @param ip_data Interaction program data, ie. the output of either `InteractionPrograms` or `FindAllInteractionPrograms`
#' @param n.replicate Number of one-sided Mann-Whitney tests to perform to assess correlation of each module. Default: 10000
#' @param min.members Minimum number of unique ligands or receptors in order to keep a module. Default: 1. Ie, if a module contains many ligand-receptor pairs that all share the same single ligand, this module will be discarded
#'
#' @return Returns a data.frame of all significant modules, their connectivity p-values for each sample, and their intramodular connectivity in each sample
#' @import pbapply dplyr
#' @importFrom plyr ldply
#' @importFrom purrr reduce
#' @export
#'
#' @examples
InteractionProgramSignificance <- function(ip_data, n.replicate = 10000, min.members = 1) {
  if(class(ip_data[[1]])[1]=="matrix") {
    m_cor_list = list(ip_data[[1]])
    tom_list = list(ip_data[[2]])
    con_list = list(ip_data[[4]])
    mod_list <- ip_data[[3]]
    sample_names = "sample"
  }
  else {
    m_cor_list <- ip_data[[1]]
    tom_list <- ip_data[[2]]
    con_list <- ip_data[[4]]
    mod_list <- ip_data[[3]]
    sample_names <- names(m_cor_list)
  }
  #Test module significance in each sample
  message("Testing module significance")
  random_connectivity_test <- function(m_cor=m_cor,mod=mod) {
    vars <- sample(colnames(m_cor),length(mod))
    a <- as.vector(m_cor[rownames(m_cor) %in% mod, colnames(m_cor) %in% mod])
    b <- as.vector(m_cor[vars,vars])
    return(wilcox.test(a,b,alternative = "greater",exact = F)$p.value)
  }
  mod_sign <- lapply(seq_along(1:length(m_cor_list)), function(x) {
    tmp_m_cor <- m_cor_list[[x]]
    tmp_mod_sign <- unlist(pblapply(seq_along(1:length(mod_list)), function(y) {
      try({
        p <- replicate(n.replicate,random_connectivity_test(m_cor = tmp_m_cor, mod = mod_list[[y]]))
        return(sum(p>0.05)/n.replicate)
      }, silent = T)
    }))
  })

  #are any modules non-significant across all samples? If so, remove.
  mod_sign_m <- t(do.call(rbind,mod_sign))
  mod_sign_m[grepl("Error",mod_sign_m)] <- 1
  mod_sign_m <- matrix(as.numeric(unlist(mod_sign_m)),nrow=nrow(mod_sign_m))

  rownames(mod_sign_m) <- names(mod_list)
  colnames(mod_sign_m) <- paste(sample_names,"pval",sep = "_")
  n_nonsig <- apply(mod_sign_m,1,function(x) {sum(x>0.05)})
  mod_list <- mod_list[n_nonsig<ncol(mod_sign_m)]

  #remove modules that contain only one ligand or receptor
  lr_counts <- data.frame(ligands=unlist(lapply(mod_list, function(x) {
    length(unique(unlist(lapply(str_split(x, pattern = "="), function(y) {y[[1]]}))))
  })),
  receptors = unlist(lapply(mod_list, function(x) {
    length(unique(unlist(lapply(str_split(x, pattern = "="), function(y) {y[[2]]}))))
  })))
  mod_list <- mod_list[lr_counts$ligands>min.members & lr_counts$receptors>min.members]

  #now return data
  #make a dataframe of module p values for each sample
  #merge this with a gene-wise modularity dataframe (remember to handle merged modules)
  range01 <- function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}
  mod_df <- reshape2::melt(plyr::ldply(mod_list, rbind), id.var = ".id") %>%
    dplyr::select(-variable) %>% dplyr::filter(!is.na(value))
  colnames(mod_df) <- c("name","lr_pair")
  # mod_df <- data.frame(lr_pair=unlist(mod_list, use.names = T)) %>%
  #   rownames_to_column("name") %>%
  #   dplyr::mutate(name = gsub('[[:digit:]]+', '', name))
  mod_df <- merge(mod_df, as.data.frame(mod_sign_m) %>%
                    rownames_to_column(var = "name"), by = "name", all.y = F)
  con_sum <- as.matrix(lapply(1:length(m_cor_list), function(x) {
    y <- as.data.frame(con_list[[x]]) %>%
      rownames_to_column(var = "lr_pair") %>% dplyr::select(lr_pair,kTotal)
    colnames(y) <- c("lr_pair",
                     paste(sample_names[x],
                           "connectivity",sep = "_"))
    return(y)
  }) %>% purrr::reduce(full_join, by = "lr_pair") %>%
    column_to_rownames(var = "lr_pair") %>%
    mutate_all(range01))
  con_sum[is.na(con_sum)] <- 0
  con_sum <- con_sum %>% as.data.frame() %>% rownames_to_column("lr_pair")
  mod_df <- merge(mod_df, con_sum, by = "lr_pair", all.y = F) %>%
    arrange(name)
  return(mod_df)

}

#' Score expression of single-cells by expression of discovered interaction programs
#'
#' @param mods Interaction program data. Either the data.frame output of `InteractionProgramSignificance`, or a list of interaction program genes, as in the output of `InteractionPrograms` or `FindAllInteractionPrograms`
#' @param object A Seurat object to be scored
#' @param return.assay
#'
#' @return A Seurat object. When `return.assay = R`, returns interaction program ligand scores as an assay `IP_ligands` and receptor scores as an assay `IP_receptors`.
#' When `return.assay = F`, Returns a Seurat object with interaction program scores as meta.data columns. Ligand and receptor expression of interaction programs are scored separately.
#' The scores for sender cells (ligand expression) are stored in columns named "ligands_[interaction program name]"
#' The scores for receiver cells (receptor expression) are stored in columns named "receptors_[interaction program name]"
#' @import dplyr stringr Seurat
#' @export
#'
#' @examples
ScoreInteractionPrograms <- function(object, mods, return.assay = T) {
  seu <- object
  if(class(mods)=="data.frame") {
    mod_names <- unique(mods$name)
    mods <- lapply(seq_along(1:length(mod_names)), function(x) {
      mods <- mods %>% dplyr::filter(name==mod_names[x]) %>% pull(lr_pair)
    })
    names(mods) <- mod_names
  }

  l_mods <- lapply(mods, function(x) {
    unique(unlist(lapply(str_split(x, pattern = "="), function(y) {y[[1]]})))
  })
  r_mods <- lapply(mods, function(x) {
    unique(unlist(lapply(str_split(x, pattern = "="), function(y) {y[[2]]})))
  })
  message("Scoring ligands")
  seu <- AddModuleScore(seu, features = l_mods, name = "ligands")
  message("Scoring receptors")
  seu <- AddModuleScore(seu, features = r_mods, name = "receptors")

  colnames(seu@meta.data)[grepl("^ligands",colnames(seu@meta.data))] <- paste("ligands",names(mods),sep = "_")
  colnames(seu@meta.data)[grepl("^receptors",colnames(seu@meta.data))] <- paste("receptors",names(mods),sep = "_")
  if(return.assay) {
    ligand_scores <- as.matrix(seu@meta.data[,grepl("^ligands",colnames(seu@meta.data))])
    colnames(ligand_scores) <- gsub("ligands_","",colnames(ligand_scores))
    ligand_assay <- CreateAssayObject(data = t(ligand_scores))

    receptor_scores <- as.matrix(seu@meta.data[,grepl("^receptors",colnames(seu@meta.data))])
    colnames(receptor_scores) <- gsub("receptors_","",colnames(receptor_scores))
    receptor_assay <- CreateAssayObject(data = t(receptor_scores))

    object[["IPligands"]] <- ligand_assay
    object[["IPreceptors"]] <- receptor_assay
    return(object)
  } else {
    return(seu)
  }

}



#' Identify most highly expressed interaction programs by cell type
#'
#' @param seu A seurat object where interaction program expression has been scored by `ScoreInteractionPrograms`
#' @param group.by `meta.data` column corresponding to the cell type or cluster annotations to use for aggregating Interaction program expression values
#'
#' @return A data.frame with average interaction program expression score for each cell type combination.
#' @import dplyr tibble
#' @export
#'
#' @examples
IPCellTypeSummary <- function(seu, group.by) {
  #matrix addition at baseline will represent autocrine signaling.
  #by shuffling the row order in one matrix, we can calculate high scores across cells
  ip_interact = list()

  ip_lig <- as.matrix(seu[["IPligands"]]@data %>% t() %>%
                                as.data.frame() %>% add_column(celltype = seu@meta.data[,group.by]) %>%
                                group_by(celltype) %>%
                                summarise_if(is.numeric, mean) %>% column_to_rownames("celltype"))

  # ip_lig <- as.matrix(seu@meta.data %>%
  #                               select(group.by,
  #                                      starts_with("ligands")) %>%
  #                               group_by((!!sym(group.by))) %>%
  #                               summarise_if(is.numeric, mean) %>% column_to_rownames("celltype"))
  ip_rec <- as.matrix(seu[["IPreceptors"]]@data %>% t() %>%
                        as.data.frame() %>% add_column(celltype = seu@meta.data[,group.by]) %>%
                        group_by(celltype) %>%
                        summarise_if(is.numeric, mean) %>% column_to_rownames("celltype"))

  row_orders = sapply(0:(nrow(ip_lig)-1), function(x) c((1 + x):nrow(ip_lig), seq_len(x))) #create matrix of row shufflings
  for (i in 1:nrow(ip_lig)){
    ip_interact[[i]] = ip_lig[row_orders[i,],] + ip_rec
  }

  #pick out indices in each where the sum is > 0.5 and export table
  interact.cutoff = 0
  all.poi.list = lapply(ip_interact, function(mat){
    hits.ind = which(mat >interact.cutoff, arr.ind = TRUE)
    hits = data.frame(sender = rownames(mat)[hits.ind[,1]],
                      receiver = rownames(ip_rec)[hits.ind[,1]],
                      program = gsub("ligands_","", x = colnames(mat)[hits.ind[,2]]),
                      additive.score = apply(hits.ind, 1, function(inds) mat[inds[1],inds[2]]))
    return(hits)
  })
  all.poi = do.call(rbind, all.poi.list)
  return(all.poi %>% as_tibble() %>% arrange(-additive.score))
}




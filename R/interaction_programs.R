
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
#' @param softPower softPower threshold for adjacency matrix
#' @param min.size Minimum size of each interaction program
#' @param plot.mods Plot modules and associated dendrograms during analysis
#' @param tree.cut.quantile The dendrogram tree height quantile at which the dendrogram should be cut. Higher values lead to fewer, smaller modules.
#'
#' @return When return.mat = T, returns a list of length 4 containing ligand-receptor covariance matrix, TOM, module lists, and intramodular connectivity. Otherwise, returns a list of length 2 containing only module lists and intramodular connectivity.
#' @import qlcMatrix WGCNA flashClust dynamicTreeCut reshape2 pbapply
#' @export
#'
#' @examples
InteractionPrograms <- function(object, assay = "SCT", slot = "data",
                                species = "human", database = "OmniPath",
                                ligands = NULL, recepts = NULL,
                                iterate.threshold = 500, n.iterate = NULL,
                                specific = F, ranked_genes = NULL,
                                return.mat = T, softPower = 1,
                                min.size = 5, plot.mods = F,
                                tree.cut.quantile = 0.4) {
  if(database=="custom") {
    if(is.null(ligands) | is.null(recepts)) {
      stop("To use custom database, please supply equidimensional character vectors of ligands and recepts")
    }
    message("Using custom database")
    ligands <- ligands
    recepts <- recepts
    lit.put <- data.frame(pair.name = paste(ligands,recepts,sep="_"), source_genesymbol = ligands, target_genesymbol = recepts)
  }
  if((!is.null(ligands) | !is.null(recepts)) & database != "custom") {
    stop("To use custom ligand or receptor lists, set database = 'custom'")
  }
  else {
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
    mat_list <- lapply(seq_along(1:n.rep), function(z) {
      cell.exprs.sub <- as.data.frame(cell.exprs[,sample(colnames(cell.exprs),iterate.threshold)]) %>% rownames_to_column(var = "gene")
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
    # mat_list <- lapply(m_list, function(m) {
    #
    # })
    #filter to intersection
    cor_names <- Reduce(intersect,lapply(mat_list,colnames))
    # m <- do.call(rbind,lapply(m_list, function(x) {x[,cor_names]}))
    # m <- m[sample(1:nrow(m),iterate.threshold^2),]

    mat_list <- lapply(mat_list, function(x) {x[cor_names,cor_names]})
    m_cor <- apply(simplify2array(mat_list), 1:2, median)
  }
  else {
    message(paste("\nGenerating Interaction Matrix..."))
    m <- sqrt(as.sparse((pbsapply(1:nrow(a), function(i) tcrossprod(a[i, ], b[i, ])))))
    colnames(m) <- paste(cell.exprs.lig$ligands, cell.exprs.rec$recepts, sep = "=")
    cna <- rep(colnames(object),ncol(object))
    cnb <- rep(colnames(object),each=ncol(object))
    rownames(m) <- paste(cna,cnb,sep = "=")
    m <- m[,Matrix::colSums(m)>0]
    m_cor <- 0.5+(0.5*corSparse(m))
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
    InteractionPrograms(object = x, return.mat = T)
  })

  mod_list <- unlist(lapply(q_mods,function(x){x[[3]]}), recursive = F)
  tom_list <- lapply(q_mods, function(x){x[[2]]})
  m_cor_list <- lapply(q_mods,function(x){x[[1]]})
  con_list <- lapply(q_mods,function(x){x[[4]]})
  names(con_list) <- names(tom_list) <- names(m_cor_list) <- names(q_mods)

  #Merge similar modules
  message("Merging similar modules")
  lrsum <- matrix(ncol = length(names(mod_list)), nrow = length(unique(unlist(mod_list))))
  colnames(lrsum) <- names(mod_list)
  rownames(lrsum) <- unique(unlist(mod_list))

  for (i in 1:nrow(lrsum)){
    for (k in colnames(lrsum)) {
      lrsum[i,k] <- ifelse(rownames(lrsum)[i] %in% mod_list[[k]],1,0)
    }
  }

  d <- as.matrix(dist.binary(t(lrsum), method = 1))
  d <- 1-d^2
  tomerge <- which(d>sim_threshold&d<1)
  partners <- unique(bind_rows(lapply(seq_along(1:length(tomerge)), function(x) {
    k <- arrayInd(tomerge[x],dim(d))
    y <- sort(c(rownames(d)[k[,1]],colnames(d)[k[,2]]))
    data.frame(a = y[1], b = y[2], ind = d[d>sim_threshold&d<1][x])
  }))) %>% arrange(-ind)
  partners <- partners[!(duplicated(partners$a)),]
  partners <- partners[!(duplicated(partners$b)),]
  for (i in 1:nrow(partners)) {
    mod_list[[length(mod_list)+1]] <- unique(c(mod_list[[partners[i,"a"]]],
                                               mod_list[[partners[i,"b"]]]))
    mod_list[partners[i,"a"]] <- NULL
    mod_list[partners[i,"b"]] <- NULL
    names(mod_list)[length(mod_list)] <- paste(partners[i,"a"],partners[i,"b"],sep = "=")
  }

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
    a <- as.vector(m_cor[mod,mod])
    b <- as.vector(m_cor[vars,vars])
    return(wilcox.test(a,b,alternative = "greater")$p.value)
  }
  mod_sign <- lapply(seq_along(1:length(m_cor_list)), function(x) {
    tmp_m_cor <- m_cor_list[[x]]
    tmp_mod_sign <- unlist(pblapply(seq_along(1:length(mod_list)), function(y) {
      # random_distribution <- replicate(n.replicate, random_connectivity(m_cor = tmp_m_cor, mod = mod_list[[y]]))
      # connectivity <- mean(tmp_m_cor[rownames(tmp_m_cor) %in% mod_list[[y]],
      #                                colnames(tmp_m_cor) %in% mod_list[[y]]])
      # sum(connectivity<random_distribution)/n.replicate
      p <- replicate(n.replicate,random_connectivity_test(m_cor = tmp_m_cor, mod = mod_list[[y]]))
      sum(p>0.05)/n.replicate
    }))
  })

  #are any modules non-significant across all samples? If so, remove.
  mod_sign_m <- t(do.call(rbind,mod_sign))
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
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  mod_df <- reshape2::melt(plyr::ldply(mod_list, rbind), id.var = ".id") %>%
    dplyr::select(-variable) %>% dplyr::filter(!is.na(value))
  colnames(mod_df) <- c("name","lr_pair")
  # mod_df <- data.frame(lr_pair=unlist(mod_list, use.names = T)) %>%
  #   rownames_to_column("name") %>%
  #   dplyr::mutate(name = gsub('[[:digit:]]+', '', name))
  mod_df <- merge(mod_df, as.data.frame(mod_sign_m) %>%
                    rownames_to_column(var = "name"), by = "name", all.y = F)
  con_sum <- lapply(1:length(m_cor_list), function(x) {
    y <- as.data.frame(con_list[[x]]) %>%
      rownames_to_column(var = "lr_pair") %>% dplyr::select(lr_pair,kTotal)
    colnames(y) <- c("lr_pair",
                     paste(sample_names[x],
                           "connectivity",sep = "_"))
    return(y)
  }) %>% purrr::reduce(full_join, by = "lr_pair") %>%
    column_to_rownames(var = "lr_pair") %>%
    mutate_all(range01) %>% rownames_to_column(var = "lr_pair")
  mod_df <- merge(mod_df, con_sum, by = "lr_pair", all.y = F) %>%
    arrange(name)
  return(mod_df)

}

#' Score expression of single-cells by expression of discovered interaction programs
#'
#' @param seu A Seurat object
#' @param mods Interaction program data. Either the data.frame output of `InteractionProgramSignificance`, or a list of interaction program genes, as in the output of `InteractionPrograms` or `FindAllInteractionPrograms`
#'
#' @return Returns a Seurat object with interaction program scores as meta.data columns. Ligand and receptor expression of interaction programs are scored separately.
#' The scores for sender cells (ligand expression) are stored in columns named "ligands_[interaction program name]"
#' The scores for receiver cells (receptor expression) are stored in columns named "receptors_[interaction program name]"
#' @import dplyr stringr Seurat
#' @export
#'
#' @examples
ScoreInteractionPrograms <- function(seu, mods) {
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
  return(seu)
}

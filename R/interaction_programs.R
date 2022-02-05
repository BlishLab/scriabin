
#' Title
#'
#' @param object
#' @param assay
#' @param slot
#' @param database
#' @param ligands
#' @param recepts
#' @param iterate.threshold
#' @param n.iterate
#' @param specific
#' @param ranked_genes
#' @param return.mat
#' @param softPower
#' @param min.size
#' @param plot.mods
#' @param tree.cut.quantile
#'
#' @return
#' @import qlcMatrix WGCNA flashClust dynamicTreeCut reshape2
#' @export
#'
#' @examples
InteractionPrograms <- function(object, assay = "SCT", slot = "data",
                               database = "OmniPath", ligands = NULL,
                               recepts = NULL, iterate.threshold = 500,
                               n.iterate = NULL,
                               specific = F, ranked_genes = NULL,
                               return.mat = T, softPower = 1,
                               min.size = 5, plot.mods = F,
                               tree.cut.quantile = 0.4) {
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

  a <- as.matrix(cell.exprs.lig[,3:ncol(cell.exprs.lig)])
  a[is.na(a)] <- 0
  b <- as.matrix(cell.exprs.rec[,3:ncol(cell.exprs.rec)])
  b[is.na(b)] <- 0

  if(ncol(object)>iterate.threshold) {
    message("\nIteratively generating interaction matrix")
    if(is.null(n.iterate)) {
      n.rep <- round(sqrt(ncol(object)/iterate.threshold))
    }
    else {
      n.rep = n.iterate
    }
    mat_list <- pblapply(seq_along(1:n.rep), function(z) {
      x <- a[,sample(1:ncol(a),iterate.threshold)]
      y <- b[,colnames(x)]
      m <- sqrt(as.sparse((sapply(1:nrow(x), function(i) tcrossprod(x[i, ], y[i, ])))))
      colnames(m) <- paste(cell.exprs.lig$ligands, cell.exprs.rec$recepts, sep = "=")
      cna <- rep(colnames(x),ncol(x))
      cnb <- rep(colnames(x),each=ncol(x))
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



#' Title
#'
#' @param seu
#' @param group.by
#' @param n.replicate
#' @param min.members
#' @param return.mats
#' @param sim_threshold
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
FindAllInteractionPrograms <- function(seu, group.by = NULL,
                           return.mats = F, sim_threshold = 0.15,
                           ...) {
  if(is.null(group.by)) {
    message("Grouping by current idents")
    seu$mod_grouping <- Idents(seu)
    group.by = "mod_grouping"
  }

  seu_split <- SplitObject(qseu, split.by = group.by)
  q_mods <- lapply(seu_split, function(x) {
    InteractionPrograms(object = x, return.mat = T, tree.cut.quantile = 0.2)
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
    data.frame(a = y[1], b = y[2], ind = d[d>0.15&d<1][x])
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


#' Title
#'
#' @param ip_data
#' @param n.replicate
#' @param min.members
#'
#' @return
#' @import pbapply plyr dplyr purrr
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
  mod_list <- mod_list[lr_counts$ligands>1 & lr_counts$receptors>1]

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

#' Title
#'
#' @param seu
#' @param mods
#'
#' @return
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

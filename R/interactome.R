

#' Title
#'
#' @param seu
#' @param gene_rankings
#' @param interactome
#' @param min.pct
#' @param assay.use
#' @param slot.use
#' @param send_cells
#' @param rec_cells
#'
#' @return
#' @export
#'
#' @examples
PrioritizeInteractome <- function(seu, gene_rankings, interactome = NULL,
                                  min.pct = 0.05, assay.use = "RNA", slot.use = "counts",
                                  send_cells = NULL, rec_cells = NULL) {
  if(is.null(send_cells)) {
    send_cells = unique(interactome$source)
  }
  if(is.null(rec_cells)) {
    rec_cells = unique(interactome$receiver)
  }

  s.exprs <- GetAssayData(seu, assay = assay.use, slot = slot.use)[,send_cells]
  r.exprs <- GetAssayData(seu, assay = assay.use, slot = slot.use)[,rec_cells]

  s.expressed_genes <- rownames(s.exprs)[(Matrix::rowSums(s.exprs !=0)/ncol(s.exprs))>min.pct]
  r.expressed_genes <- rownames(r.exprs)[(Matrix::rowSums(r.exprs !=0)/ncol(r.exprs))>min.pct]

  background_expressed_genes <- r.expressed_genes %>% .[. %in% rownames(ligand_target_matrix)]
  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()
  expressed_ligands = intersect(ligands,s.expressed_genes)
  expressed_receptors = intersect(receptors,r.expressed_genes)
  potential_ligands = lr_network %>% dplyr::filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

  gene_rankings = gene_rankings[names(gene_rankings) %in% rec_cells]
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
#' @param ranked_genes
#' @param send_cells
#' @param rec_cells
#' @param clusters
#' @param send_cluster
#' @param rec_cluster
#' @param name
#' @param cell.type.calls
#' @param ident.label
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ExploreClusters <- function(seu, ranked_genes, send_cells = NULL, rec_cells = NULL,
                            clusters = NULL, send_cluster = NULL, rec_cluster = NULL,
                            name = NULL, cell.type.calls = "celltype.l2",
                            ident.label = NULL, ...) {
  if(is.null(clusters)) {
    clusteranalyze.interactome <- GenerateInteractome_new(seu = seu, ranked_genes = ranked_genes, send_cells = send_cells, rec_cells = rec_cells,
                                                      ident.label = ident.label, cell.type.calls = cell.type.calls, group.by = group.by, use_clusters = F)

  }
  else {
    clusteranalyze.interactome <- GenerateInteractome_new(seu = seu, ranked_genes = ranked_genes,
                                                      cluster_results = clusters, send_cluster = send_cluster, rec_cluster = rec_cluster,
                                                      ident.label = ident.label, cell.type.calls = cell.type.calls, group.by = group.by, use_clusters = T)

  }

  nichenet.results_analyze <- PrioritizeInteractome(seu, ranked_genes, clusteranalyze.interactome)
  nn.ligands <- bind_rows(lapply(nichenet.results_analyze, function(x) {x[[1]]}), .id = "cell") %>%
    filter(pearson>0.075)
  nn.ligands$cell_ligand <- paste(nn.ligands$cell,nn.ligands$test_ligand,sep = "_")
  nn.genes <- bind_rows(lapply(nichenet.results_analyze, function(x) {x[[2]]}), .id = "cell")
  nn.genes$cell_ligand <- paste(nn.genes$cell,nn.genes$ligand,sep = "_")
  nn.merge <- merge(nn.ligands,nn.genes[,c("cell_ligand","target","weight")],by = "cell_ligand", all.x = T, all.y = T)
  nn.merge <- nn.merge[!(is.na(nn.merge$pearson)),]
  colnames(nn.merge)[colnames(nn.merge)=="weight"] <- "geneweight"

  interactome = clusteranalyze.interactome
  interactome$perturbed <- ifelse(interactome$lig_geneset==T | interactome$rec_geneset==T,T,F)
  interactome_filtered <- interactome %>% dplyr::filter(perturbed==T)
  interactome_filtered$cell_ligand <- paste(interactome_filtered$receiver,interactome_filtered$ligand,sep = "_")

  interactome_merge <- merge(interactome_filtered,nn.merge,by = "cell_ligand")
  interactome_merge$weight <- interactome_merge$edgeweight * interactome_merge$geneweight
  interactome_merge <- interactome_merge[!is.na(interactome_merge$target),]

  # PlotAlluvium(interactome_merge)
  # saveit(format = "png", height = 8, width = 16, name = name)

  return(list(base_interactome=clusteranalyze.interactome,
              nichenet_results = nichenet.results_analyze,
              nichenet_interactome = interactome_merge))
}


#' Title
#'
#' @param seu
#' @param ranked_genes
#' @param cluster_results
#' @param send_cluster
#' @param rec_cluster
#' @param ident.label
#' @param use_clusters
#' @param send_cells
#' @param rec_cells
#' @param cell.type.calls
#' @param assay.use
#' @param slot.use
#' @param group.by
#' @param database
#'
#' @return
#' @export
#'
#' @examples
GenerateInteractome <- function(seu, ranked_genes = ranked_genes, cluster_results = clusters,
                                send_cluster = NULL, rec_cluster = NULL, ident.label = NULL, use_clusters = F,
                                send_cells = NULL, rec_cells = NULL, cell.type.calls = "celltype.l2",
                                assay.use = "SCT", slot.use = "data", group.by = "time.orig", database = "OmniPath") {
  if(use_clusters) {
    if(!is.null(send_cells) | !is.null(rec_cells)) {
      stop("Cell names are supplied both through cluster results and send_cells/rec_cells.
           If you wish to supply cell names explicitly, set use_clusters = F")
    }
    message("Using clusters")
    message(ident.label)
    send_nb = names(cluster_results$SCluster)[cluster_results$SCluster %in% send_cluster]
    rec_nb = names(cluster_results$RCluster)[cluster_results$RCluster %in% rec_cluster]

    send_cells = colnames(seu)[seu$bins %in% send_nb & (seu@meta.data[,group.by] %in% ident.label)]
    rec_cells = colnames(seu)[seu$bins %in% rec_nb & (seu@meta.data[,group.by] %in% ident.label)]
  }
  else {
    if(is.null(send_cells) | is.null(rec_cells)) {
      stop("Must supply cells to analyze either in cluster results, or as vectors send_cells and rec_cells")
    }
    message("Using explicitly supplied cell names")
  }
  # send_genes = ranked_genes[send_cells]
  # names(send_genes) <- NULL
  # ligands_use = unique(names(unlist(send_genes)))
  # names(send_genes) <- send_cells
  #
  # rec_genes = ranked_genes[rec_cells]
  # names(rec_genes) <- NULL
  # recept_use = unique(names(unlist(rec_genes)))
  # names(rec_genes) <- rec_cells

  if(database=="custom") {
    message("Using custom database")
    ligands <- ligands
    recepts <- recepts
    pairs <- data.frame(ligands = ligands, recepts = recepts)
  }
  else {
    all <- readRDS(system.file(package = "scriabin", "lr_resources.rds"))
    if(database %notin% names(all)) {
      stop("Database must be one of: OmniPath, CellChatDB, CellPhoneDB, Ramilowski2015, Baccin2019, LRdb, Kirouac2010, ICELLNET, iTALK, EMBRACE, HPMR, Guide2Pharma, connectomeDB2020, talklr, CellTalkDB")
    }
    message(paste("Using database",database))
    pairs <- as.data.frame(all[[database]][,c("source_genesymbol","target_genesymbol")] %>% dplyr::mutate_all(as.character))
    ligands <- as.character(lit.put[, "source_genesymbol"])
    recepts <- as.character(lit.put[, "target_genesymbol"])
  }
  colnames(pairs) <- c("ligand","receptor")

  cell.exprs <- GetAssayData(seu, assay = assay.use, slot = slot.use)[rownames(seu) %in% c(unique(pairs$ligand),unique(pairs$receptor)),c(send_cells,rec_cells)]
  cell.exprs <- cell.exprs[Matrix::rowSums(cell.exprs)>0,]
  # cell.exprs <- as.data.frame(cell.exprs) %>% rownames_to_column("gene")
  sources <- send_cells
  targets <- rec_cells

  message(paste0(length(sources)," sources"))
  message(paste0(length(targets)," targets"))

  message(paste("\nGenerating Interactome"))
  pb <- txtProgressBar(min = 0, max = length(sources), initial = 0,
                       style = 3)
  connectome <- data.frame()
  for (i in 1:length(sources)) {
    temp <- data.frame()
    for (j in 1:length(targets)) {
      # pairs_use <- pairs[pairs$ligand %in% names(send_genes[[sources[i]]]) & pairs$receptor %in% names(rec_genes[[targets[j]]]),]
      cell.exprs_use <- cell.exprs[,c(sources[i],targets[j])]
      cell.exprs_use <- as.data.frame(cell.exprs_use[Matrix::rowSums(cell.exprs_use)>0,]) %>% rownames_to_column("gene")
      pairs_use <- pairs %>%
        dplyr::filter(ligand %in% cell.exprs_use$gene) %>%
        dplyr::filter(receptor %in% cell.exprs_use$gene)
      if(nrow(cell.exprs_use)>0 & sources[i]!=targets[j] & nrow(pairs_use)>0) {
        lig_i <- pairs_use[,"ligand"]
        recept_j <- pairs_use[,"receptor"]
        ligands.df <- data.frame(ligands = lig_i)
        ligands.df$id <- 1:nrow(ligands.df)
        recepts.df <- data.frame(recepts = recept_j)
        recepts.df$id <- 1:nrow(recepts.df)
        cell.exprs_use.rec <- merge(recepts.df, cell.exprs_use[,c("gene",rec_cells[j])],
                                    by.x = "recepts", by.y = "gene", all.x = T)
        cell.exprs_use.rec <- cell.exprs_use.rec[order(cell.exprs_use.rec$id),
                                                 ]
        cell.exprs_use.lig <- merge(ligands.df, cell.exprs_use[,c("gene",send_cells[i])],
                                    by.x = "ligands", by.y = "gene", all.x = T)
        cell.exprs_use.lig <- cell.exprs_use.lig[order(cell.exprs_use.lig$id),
                                                 ]
        lig_geneset <- ifelse(lig_i %in% names(ranked_genes[[sources[i]]]),T,F)
        rec_geneset <- ifelse(recept_j %in% names(ranked_genes[[targets[j]]]),T,F)
        vector <- data.frame(source = sources[i], receiver = targets[j],
                             ligand = lig_i, receptor = recept_j,
                             pair = paste(lig_i, recept_j, sep = " - "),
                             ligand.expression = cell.exprs_use.lig[, sources[i]],
                             recept.expression = cell.exprs_use.rec[, targets[j]],
                             lig_geneset = lig_geneset, rec_geneset = rec_geneset#,
                             # ligand.scale = cluster.avgs.scale.df.lig[, sources[i]],
                             # recept.scale = cluster.avgs.scale.df.rec[, targets[j]],
                             # percent.source = cluster.pcts.df.lig[, sources[i]],
                             # percent.target = cluster.pcts.df.rec[, targets[j]]
        )
        temp <- rbind(temp, vector)
      }
    }
    connectome <- rbind(connectome, temp)
    Sys.sleep(0.5)
    setTxtProgressBar(pb, i)
  }
  connectome$ligand.expression <- log1p(connectome$ligand.expression)
  connectome$recept.expression <- log1p(connectome$recept.expression)
  connectome$edgeweight <- connectome$ligand.expression*connectome$recept.expression
  connectome <- connectome[connectome$edgeweight>0,]

  source_types <- data.frame(source_type = seu@meta.data[,cell.type.calls])
  source_types$source <- colnames(seu)
  target_types <- data.frame(receiver_type = seu@meta.data[,cell.type.calls])
  target_types$receiver <- colnames(seu)

  connectome <- merge(connectome,target_types,by = "receiver",all.x=T,all.y=F)
  connectome <- merge(connectome,source_types,by = "source",all.x=T,all.y=F)


  return(connectome)
}



#' Title
#'
#' @param seu
#' @param ranked_genes
#' @param send_cells
#' @param rec_cells
#' @param name
#'
#' @return
#' @export
#'
#' @examples
ExploreInteractome <- function(seu, ranked_genes, send_cells, rec_cells, name = NULL) {
  clusteranalyze.interactome <- GenerateInteractome(seu, ranked_genes, use_clusters = F, send_cells = send_cells, rec_cells = rec_cells)

  nichenet.results_analyze <- PrioritizeInteractome(seu, ranked_genes, clusteranalyze.interactome)
  nn.ligands <- bind_rows(lapply(nichenet.results_analyze, function(x) {x[[1]]}), .id = "cell") %>%
    filter(pearson>0.075)
  nn.ligands$cell_ligand <- paste(nn.ligands$cell,nn.ligands$test_ligand,sep = "_")
  nn.genes <- bind_rows(lapply(nichenet.results_analyze, function(x) {x[[2]]}), .id = "cell")
  nn.genes$cell_ligand <- paste(nn.genes$cell,nn.genes$ligand,sep = "_")
  nn.merge <- merge(nn.ligands,nn.genes[,c("cell_ligand","target","weight")],by = "cell_ligand", all.x = T, all.y = T)
  nn.merge <- nn.merge[!(is.na(nn.merge$pearson)),]
  colnames(nn.merge)[colnames(nn.merge)=="weight"] <- "geneweight"

  interactome = clusteranalyze.interactome
  interactome$perturbed <- ifelse(interactome$lig_geneset==T | interactome$rec_geneset==T,T,F)
  interactome_filtered <- interactome %>% dplyr::filter(perturbed==T)
  interactome_filtered$cell_ligand <- paste(interactome_filtered$receiver,interactome_filtered$ligand,sep = "_")

  interactome_merge <- merge(interactome_filtered,nn.merge,by = "cell_ligand")
  interactome_merge$weight <- interactome_merge$edgeweight * interactome_merge$geneweight
  interactome_merge <- interactome_merge[!is.na(interactome_merge$target),]

  PlotAlluvium(interactome_merge)
  saveit(format = "png", height = 8, width = 16, name = name)

  return(list(base_interactome=clusteranalyze.interactome,
              nichenet_results = nichenet.results_analyze,
              nichenet_interactome = interactome_merge))
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
#'
#' @return
#' @export
#'
#' @examples
BuildWeightedInteraction_test <- function (object, nichenet_results = late1.nnr, assay = "SCT", slot = "data",
                                           pearson.cutoff = 0.1, scale.factors = c(1.5,3),
                                           database = "fantom5", ligands = NULL, recepts = NULL) {
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(plotrix)
  if(database=="custom") {
    message("Using custom database")
    ligands <- ligands
    recepts <- recepts
    lit.put <- data.frame(pair.name = paste(ligands,recepts,sep="_"), ligands = ligands, ligand_name = ligands, recepts = recepts)
  }
  else {
    message("Using fantom5")
    fantom5 <- Connectome::ncomms8866_human
    # object <- subset(object, cells = colnames(object)[object$time.orig %in% tps])
    lit.put <- fantom5[fantom5$Pair.Evidence %in% c("literature supported",
                                                    "putative"), ]
    ligands <- as.character(lit.put[, 2])
    recepts <- as.character(lit.put[, 4])
  }
  ligands.use <- intersect(ligands, rownames(object@assays[[assay]]))
  recepts.use <- intersect(recepts, rownames(object@assays[[assay]]))
  genes.use = union(ligands.use, recepts.use)
  cell.exprs <- as.data.frame(GetAssayData(object, assay = assay, slot = slot)[genes.use,]) %>% rownames_to_column(var = "gene")
  ligands.df <- data.frame(lit.put[, c(1,2)]) %>% dplyr::mutate_all(as.character)
  colnames(ligands.df) <- c("pair","ligands")
  ligands.df$id <- 1:nrow(ligands.df)
  recepts.df <- data.frame(lit.put[, c(1,4)]) %>% dplyr::mutate_all(as.character)
  colnames(recepts.df) <- c("pair","recepts")
  recepts.df$id <- 1:nrow(recepts.df)

  ###
  cell.exprs_rec <- as.data.frame(GetAssayData(object, assay = assay, slot = slot)[genes.use,colnames(object) %in% names(nichenet_results)]) %>% rownames_to_column(var = "gene")
  ###

  cell.exprs.rec <- merge(recepts.df, cell.exprs_rec,
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
  int.to.merge <- lit.put[,c(1,2,4)]
  colnames(int.to.merge) <- c("pair","ligand","recepts")
  test <- merge(test,int.to.merge,by = "ligand", all.x=T)
  test <- test[!is.na(test$recepts),]

  if(nrow(test)==0) {
    warning("No weights found.")
    weighted.rec <- cell.exprs.rec
  }
  else {
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
  }

  message(paste("\nGenerating Interaction Matrix..."))
  m <- outer (
    cell.exprs.lig[,4:ncol(cell.exprs.lig)],     # First dimension:  the rows     (x)
    weighted.rec[,4:ncol(weighted.rec)],     # Second dimension: the columns  (y)
    Vectorize(function (x, y)   sum(x*y,na.rm=T))
    ## Rows represent "senders" and columns represent "receivers". Eg. m[1,2] is the interaction potential between cell #1 as a sender and cell #2 as a receiver.
  )
  results <- as.Graph(m)
  object@graphs$weighted_interaction <- results
  return(object)
}




#' Title
#'
#' @param seu
#' @param ranked_genes
#' @param cluster_results
#' @param send_cluster
#' @param rec_cluster
#' @param ident.label
#' @param use_clusters
#' @param send_cells
#' @param rec_cells
#' @param cell.type.calls
#' @param assay.use
#' @param slot.use
#' @param group.by
#' @param nichenet_results
#' @param pearson.cutoff
#' @param database
#'
#' @return
#' @export
#'
#' @examples
GeneratePrioritizedInteractome <- function(seu, ranked_genes = ranked_genes, cluster_results = clusters,
                                    send_cluster = NULL, rec_cluster = NULL, ident.label = NULL, use_clusters = T,
                                    send_cells = NULL, rec_cells = NULL, cell.type.calls = "celltype.l2",
                                    assay.use = "SCT", slot.use = "data", group.by = "time.orig",
                                    nichenet_results, pearson.cutoff = 0.075, database = "OmniPath") {
  if(use_clusters) {
    if(!is.null(send_cells) | !is.null(rec_cells)) {
      stop("Cell names are supplied both through cluster results and send_cells/rec_cells.
           If you wish to supply cell names explicitly, set use_clusters = F")
    }
    message("Using clusters")
    message(ident.label)
    send_nb = names(cluster_results$SCluster)[cluster_results$SCluster %in% send_cluster]
    rec_nb = names(cluster_results$RCluster)[cluster_results$RCluster %in% rec_cluster]

    send_cells = colnames(seu)[seu$bins %in% send_nb & (seu@meta.data[,group.by] %in% ident.label)]
    rec_cells = colnames(seu)[seu$bins %in% rec_nb & (seu@meta.data[,group.by] %in% ident.label)]
  }
  else {
    if(is.null(send_cells) | is.null(rec_cells)) {
      stop("Must supply cells to analyze either in cluster results, or as vectors send_cells and rec_cells")
    }
    message("Using explicitly supplied cell names")
  }

  if(database=="custom") {
    message("Using custom database")
    ligands <- ligands
    recepts <- recepts
    pairs <- data.frame(ligands = ligands, recepts = recepts)
  }
  else {
    all <- readRDS(system.file(package = "scriabin", "lr_resources.rds"))
    if(database %notin% names(all)) {
      stop("Database must be one of: OmniPath, CellChatDB, CellPhoneDB, Ramilowski2015, Baccin2019, LRdb, Kirouac2010, ICELLNET, iTALK, EMBRACE, HPMR, Guide2Pharma, connectomeDB2020, talklr, CellTalkDB")
    }
    message(paste("Using database",database))
    pairs <- as.data.frame(all[[database]][,c("source_genesymbol","target_genesymbol")] %>% dplyr::mutate_all(as.character))
    ligands <- as.character(pairs[, "source_genesymbol"])
    recepts <- as.character(pairs[, "target_genesymbol"])
  }
  colnames(pairs) <- c("ligand","receptor")

  cell.exprs <- GetAssayData(seu, assay = assay.use, slot = slot.use)[rownames(seu) %in% c(unique(pairs$ligand),unique(pairs$receptor)),c(send_cells,rec_cells)]
  cell.exprs <- cell.exprs[Matrix::rowSums(cell.exprs)>0,]
  # cell.exprs <- as.data.frame(cell.exprs) %>% rownames_to_column("gene")
  sources <- send_cells
  targets <- rec_cells

  message(paste0(length(sources)," sources"))
  message(paste0(length(targets)," targets"))

  nnr_oi_ligands <- bind_rows(lapply(nichenet_results, function(x) {x[[1]]}), .id = "cell") %>% dplyr::filter(pearson>pearson.cutoff)

  connectome <- bind_rows(pblapply(seq_along(1:length(sources)), function(i) {
    df <- bind_rows(lapply(seq_along(1:length(targets)), function(j) {
      cell.exprs_use <- cell.exprs[,c(sources[i],targets[j])]
      cell.exprs_use <- as.data.frame(cell.exprs_use[Matrix::rowSums(cell.exprs_use)>0,]) %>% rownames_to_column("gene")
      pairs_use <- pairs %>%
        dplyr::filter(ligand %in% cell.exprs_use$gene) %>%
        dplyr::filter(receptor %in% cell.exprs_use$gene) %>%
        dplyr::filter(ligand %in% (nnr_oi_ligands %>% dplyr::filter(cell==targets[j]) %>% pull(test_ligand)))
      if(nrow(cell.exprs_use)>0 & sources[i]!=targets[j] & nrow(pairs_use)>0) {
        lig_i <- pairs_use[,"ligand"]
        recept_j <- pairs_use[,"receptor"]
        ligands.df <- data.frame(ligands = lig_i)
        ligands.df$id <- 1:nrow(ligands.df)
        recepts.df <- data.frame(recepts = recept_j)
        recepts.df$id <- 1:nrow(recepts.df)
        cell.exprs_use.rec <- merge(recepts.df, cell.exprs_use[,c("gene",rec_cells[j])],
                                    by.x = "recepts", by.y = "gene", all.x = T)
        cell.exprs_use.rec <- cell.exprs_use.rec[order(cell.exprs_use.rec$id),
                                                 ]
        cell.exprs_use.lig <- merge(ligands.df, cell.exprs_use[,c("gene",send_cells[i])],
                                    by.x = "ligands", by.y = "gene", all.x = T)
        cell.exprs_use.lig <- cell.exprs_use.lig[order(cell.exprs_use.lig$id),
                                                 ]
        lig_geneset <- ifelse(lig_i %in% names(ranked_genes[[sources[i]]]),T,F)
        rec_geneset <- ifelse(recept_j %in% names(ranked_genes[[targets[j]]]),T,F)
        vector <- data.frame(source = sources[i], receiver = targets[j],
                             ligand = lig_i, receptor = recept_j,
                             pair = paste(lig_i, recept_j, sep = " - "),
                             ligand.expression = cell.exprs_use.lig[, sources[i]],
                             recept.expression = cell.exprs_use.rec[, targets[j]],
                             lig_geneset = lig_geneset, rec_geneset = rec_geneset)
        return(vector)
      }
    }))
    return(df)
  }))
  connectome$edgeweight <- connectome$ligand.expression*connectome$recept.expression
  connectome <- connectome[connectome$edgeweight>0,]

  connectome %<>% dplyr::filter(connectome$lig_geneset==T | connectome$rec_geneset==T) %>%
    dplyr::mutate(cell_ligand = paste(receiver,ligand,sep = "_"))

  connectome$source_type <- mapvalues(connectome$source,
                                      from = colnames(seu),
                                      to = as.character(seu@meta.data[,cell.type.calls]),
                                      warn_missing = F)
  connectome$receiver_type <- mapvalues(connectome$receiver,
                                      from = colnames(seu),
                                      to = as.character(seu@meta.data[,cell.type.calls]),
                                      warn_missing = F)

  nnr_oi_targets <- bind_rows(lapply(nichenet_results, function(x) {x[[2]]}), .id = "cell") %>%
    dplyr::mutate(cell_ligand = paste(cell,ligand,sep = "_"))

  connectome <- merge(connectome, nnr_oi_targets[,c("cell_ligand","target","weight")],by = "cell_ligand", all.x = T, all.y = F)

  colnames(connectome)[colnames(connectome)=="weight"] <- "target_weight"
  connectome$cell_ligand <- NULL

  return(connectome)
}

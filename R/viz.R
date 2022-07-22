
#' Title
#'
#' @param interaction_graphs
#' @param name
#'
#' @return
#' @export
#'
#' @examples
BuildRVHeatmap <- function(interaction_graphs, name = "rv_results") {
  interaction_graphs <- lapply(interaction_graphs, function(x) {x[[1]]})
  ig_ktab <- ktab.list.df(lapply(interaction_graphs,as.data.frame))
  rv_results <- statis(ig_ktab, scannf=F)
  mat_plot <- rv_results$RV
  colnames(mat_plot) <- translateTimes(colnames(mat_plot))
  rownames(mat_plot) <- translateTimes(rownames(mat_plot))
  max_value = quantile(mat_plot,0.99, na.rm=T)
  mid_value = quantile(mat_plot,0.75, na.rm=T)
  min_value = quantile(mat_plot,0.01, na.rm=T)
  col_fun.use = colorRamp2(c(min_value,mid_value,max_value), c("#08306B","grey80","white"))
  # col_fun.use = colorRamp2(breaks = c(min_value,mid_value,2), colors = c("yellow","purple","black"))
  # pdf(paste0("~/Downloads/",name,".pdf"), height = 5, width = 7)
  Heatmap(mat_plot, cluster_rows = F, cluster_columns = F, name = "RV\ncoef", col = col_fun.use)
}


#' Title
#'
#' @param interaction_graphs
#' @param seu
#' @param name
#' @param row.k
#' @param col.k
#' @param cell.type.calls
#' @param cell.cols
#'
#' @return
#' @export
#'
#' @examples
BuildVarianceHeatmap <- function(interaction_graphs, seu, name = "var_results",
                                 row.k = NULL, col.k = NULL, cell.type.calls = "celltype.l2",
                                 cell.cols = NULL) {
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  f1 <- function(lst) {
    n <- length(lst)
    rc <- dim(lst[[1]])
    ar1 <- array(unlist(lst), c(rc, n))
    var_results <- apply(ar1, c(1, 2), sd)
    median_results <- apply(ar1, c(1, 2), median)
    median_results <- rep(median_results,n)
    time_results <- abs(ar1-median_results)
    # pre_results <- rep(lst[[1]],n)
    # time_results <- abs(ar1-pre_results)
    time_results <- apply(time_results,c(1, 2), function(y) {resample(which(y==max(y)),1)})
    time_results <- mapvalues(time_results, from = 1:n, to = names(lst))
    return(list(var_results,time_results))
  }
  results <- f1(interaction_graphs)
  var_results <- results[[1]]
  rownames(var_results) <- rownames(interaction_graphs[[1]])
  colnames(var_results) <- colnames(interaction_graphs[[1]])

  time_results <- results[[2]]
  rownames(time_results) <- rownames(interaction_graphs[[1]])
  colnames(time_results) <- colnames(interaction_graphs[[1]])

  #cluster rows
  input = var_results
  d <- dist(input, method = "euclidean")
  hc1 <- hclust(d, method = "ward")
  test <- fviz_nbclust(input, FUN = hcut, method = "wss")
  nclust = elbow(test$data %>% mutate_all(as.numeric))$clusters_selected
  if(!is.null(row.k)) {
    nclust = row.k
  }
  row_grp <- as.character(cutree(hc1, k = nclust))
  names(row_grp) <- rownames(var_results)

  #cluster columns
  input = t(var_results)
  d <- dist(input, method = "euclidean")
  hc1 <- hclust(d, method = "ward")
  test <- fviz_nbclust(input, FUN = hcut, method = "wss")
  nclust = elbow(test$data %>% mutate_all(as.numeric))$clusters_selected
  if(!is.null(col.k)) {
    nclust = col.k
  }
  col_grp <- as.character(cutree(hc1, k = nclust))
  names(col_grp) <- colnames(var_results)

  nbct <- as.data.frame(table(seu@meta.data[,cell.type.calls],seu$bins))
  colnames(nbct) <- c("celltype","nb","freq")
  nbct.f <- nbct %>% group_by(nb) %>% top_n(1,freq) %>% sample_n(1)


  row.colors <- pal_igv()(length(unique(row_grp)))
  names(row.colors) <- unique(row_grp)

  col.colors <- pal_igv()(length(unique(col_grp)))
  names(col.colors) <- unique(col_grp)

  rowanno <- nbct.f[match(rownames(var_results),as.character(nbct.f$nb)),]
  colanno <- nbct.f[match(colnames(var_results),as.character(nbct.f$nb)),]
  if(is.null(cell.cols)) {
    cell.cols <- pal_igv()(length(unique(nbct.f$celltype)))
    names(cell.cols) <- unique(nbct.f$celltype)
  }
  cell.cols_use <- cell.cols[names(cell.cols) %in% nbct.f$celltype]

  longer = resample(which(c(length(col.colors),length(row.colors))==max(length(col.colors),length(row.colors))),1)
  if(longer==1) {
    la = rowAnnotation("Sender cell" = rowanno$celltype, "SCluster" = row_grp,
                       col = list("Sender cell" = cell.cols_use,
                                  "SCluster" = row.colors),
                       annotation_legend_param = list("Sender cell" = list(labels = names(cell.cols_use),
                                                                           legend_gp = gpar(fill = cell.cols_use),
                                                                           title = "Cell type", ncol=1)),
                       show_legend = c(T,F))

    ta = columnAnnotation("Receiver cell" = colanno$celltype, "RCluster" = col_grp,
                          col = list("Receiver cell" = cell.cols_use,
                                     "RCluster" = col.colors),
                          annotation_name_side = "left",
                          show_legend = c(F,T),
                          annotation_legend_param = list("RCluster" = list(labels = names(col.colors),
                                                                           legend_gp = gpar(fill = col.colors),
                                                                           title = "Cluster #")))
  }
  if(longer==2) {
    la = rowAnnotation("Sender cell" = rowanno$celltype, "SCluster" = row_grp,
                       col = list("Sender cell" = cell.cols_use,
                                  "SCluster" = row.colors),
                       annotation_legend_param = list("Sender cell" = list(labels = names(cell.cols_use),
                                                                           legend_gp = gpar(fill = cell.cols_use),
                                                                           title = "Cell type", ncol=1),
                                                      "SCluster" = list(labels = names(row.colors),
                                                                        legend_gp = gpar(fill = row.colors),
                                                                        title = "Cluster #")))

    ta = columnAnnotation("Receiver cell" = colanno$celltype, "RCluster" = col_grp,
                          col = list("Receiver cell" = cell.cols_use,
                                     "RCluster" = col.colors),
                          annotation_name_side = "left",
                          show_legend = c(F,F))
  }


  #
  max_value = quantile(var_results,0.99)
  col_fun.use = colorRamp2(c(0,max_value), c("white","#08306B"))
  #   ct_lgd = Legend(labels = names(cell.cols), legend_gp = gpar(fill = cell.cols),
  #                title = "Cell type", ncol=2)
  #
  #   cluster_lgd = Legend(labels = names(cl_cols.use), legend_gp = gpar(fill = cl_cols.use),
  #                        title = "Cluster #")


  # hm <- Heatmap(var_results, show_row_names = F, show_column_names = F,
  #              name = "SD", top_annotation = ta, left_annotation = la, col = col_fun.use)
  #
  # hm <- draw(hm)
  #
  # ro <- row_order(hm)
  # co <- column_order(hm)
  #
  # time_results <- time_results[ro,co]
  # time_results <- translateTimes(time_results)

  hc = hclust(dist(t(var_results), method="euclidean"), method = "complete")

  hm <- Heatmap(var_results, show_row_names = F, show_column_names = F,
                cluster_columns = hc, column_dend_reorder = F,
                name = "SD", top_annotation = ta, left_annotation = la, col = col_fun.use)

  thm <- Heatmap(time_results, name = "Timepoint",
                 # col = time.cols,
                 cluster_columns = hc, column_dend_reorder = F, show_column_dend = T,
                 show_row_names = F, show_column_names = F)

  return(list(SCluster=row_grp,RCluster=col_grp,var_hm=hm,time_hm=thm))
}




#' Title
#'
#' @param interaction_graphs
#' @param seu
#' @param name
#' @param ident.1
#' @param ident.2
#' @param row.k
#' @param col.k
#'
#' @return
#' @export
#'
#' @examples
BuildCompHeatmap <- function(interaction_graphs, seu, name = "var_results", ident.1 = NULL, ident.2 = "pre",
                             row.k = NULL, col.k = NULL) {
  var_results <- abs(interaction_graphs[[ident.1]]-interaction_graphs[[ident.2]])
  rownames(var_results) <- rownames(interaction_graphs[[1]])
  colnames(var_results) <- colnames(interaction_graphs[[1]])

  #cluster rows
  input = var_results
  d <- dist(input, method = "euclidean")
  hc1 <- hclust(d, method = "ward")
  test <- fviz_nbclust(input, FUN = hcut, method = "wss")
  nclust = elbow(test$data %>% mutate_all(as.numeric))$clusters_selected
  if(!is.null(row.k)) {
    nclust = row.k
  }
  row_grp <- as.character(cutree(hc1, k = nclust))
  names(row_grp) <- rownames(var_results)

  #cluster columns
  input = t(var_results)
  d <- dist(input, method = "euclidean")
  hc1 <- hclust(d, method = "ward")
  test <- fviz_nbclust(input, FUN = hcut, method = "wss")
  nclust = elbow(test$data %>% mutate_all(as.numeric))$clusters_selected
  if(!is.null(col.k)) {
    nclust = col.k
  }
  col_grp <- as.character(cutree(hc1, k = nclust))
  names(col_grp) <- colnames(var_results)

  nbct <- as.data.frame(table(seu$celltype.l2,seu$bins))
  colnames(nbct) <- c("celltype","nb","freq")
  nbct.f <- nbct %>% group_by(nb) %>% top_n(1,freq) %>% sample_n(1)

  cell.cols <- pal_igv()(length(unique(nbct.f$celltype)))
  names(cell.cols) <- unique(nbct.f$celltype)

  row.colors <- pal_igv()(length(unique(row_grp)))
  names(row.colors) <- unique(row_grp)

  col.colors <- pal_igv()(length(unique(col_grp)))
  names(col.colors) <- unique(col_grp)

  rowanno <- nbct.f[match(rownames(var_results),as.character(nbct.f$nb)),]
  la = rowAnnotation("Sender cell" = rowanno$celltype, "SCluster" = row_grp,
                     col = list("Sender cell" = cell.cols,
                                "SCluster" = row.colors))
  colanno <- nbct.f[match(colnames(var_results),as.character(nbct.f$nb)),]
  ta = columnAnnotation("Receiver cell" = colanno$celltype, "RCluster" = col_grp,
                        col = list("Receiver cell" = cell.cols,
                                   "RCluster" = col.colors),
                        annotation_name_side = "left")
  # pdf(paste0("~/Downloads/",name,".pdf"), height = 7, width = 9.5)
  p <- Heatmap(var_results, show_row_names = T, show_column_names = T,
               name = "Diff", top_annotation = ta, left_annotation = la,
               row_names_gp = gpar(fontsize=0), column_names_gp = gpar(fontsize=0))

  return(list(SCluster=row_grp,RCluster=col_grp,Heatmap=p))
}



#' Title
#'
#' @param seu
#' @param interactome
#' @param subset.time
#'
#' @return
#' @export
#'
#' @examples
HighlightSR <- function(seu, interactome, subset.time = NULL) {
  if(!is.null(subset.time)) {
    seu <- subset(seu, cells = colnames(seu)[seu$time.orig==subset.time])
    message("Running PCA")
    seu <- RunPCA(seu, verbose = F)
    message("Running UMAP")
    seu <- RunUMAP(seu, dims = 1:50, verbose = F)
  }
  cells = list(Senders = unique(interactome$source),
               Receivers = unique(interactome$receiver))
  DimPlot(seu, cells.highlight = cells, cols.highlight = list("blue","red")) + NoLegend() + theme_umap() +
    labs(x = "UMAP1", y = "UMAP2")
}



#' Title
#'
#' @param connectome
#'
#' @return
#' @export
#'
#' @examples
PlotAlluviumLigRec <- function(connectome) {
  connectome_plot <- connectome[,c("source_type","source","ligand","receptor","target","target_type","weight")]
  colnames(connectome_plot) <- c("Sender\ncelltype","Sender\ncell","Ligand",
                                 "Receptor","Receiver\ncell","Receiver\ncelltype","Weight")
  connectome_lodes <- to_lodes_form(connectome_plot, axes = 1:6, id = "Cohort")
  connectome_lodes$stratum <- as.factor(connectome_lodes$stratum)
  ggplot(connectome_lodes, aes(x = x, stratum = stratum, alluvium = Cohort,
                               y = Weight, label = stratum, fill = Weight)) +
    geom_flow(stat = "alluvium", lode.guidance = "frontback",
              color = "darkgray") +
    geom_stratum() +
    ggrepel::geom_text_repel(aes(label = ifelse(as.numeric(x) %in% c(1,3,4,6), as.character(stratum), NA)),
                             stat = "stratum", size = 4, direction = "y") +
    theme_cowplot() + scale_fill_gradient(low = "grey50", high = "red") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + labs(x = NULL)
}

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
PrioritizeInteractome <- function(seu, gene_rankings, interactome,
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
#' @param connectome
#' @param optimize.flows
#' @param cell.cols
#'
#' @return
#' @export
#'
#' @examples
PlotAlluvium <- function(connectome, optimize.flows = T, cell.cols = NULL) {
  # if(nrow(connectome)>1e5) {
  #   connectome <- connectome %>% group_by(source,receiver,pair) %>% top_n(n=10,wt=target_weight)
  # }

  connectome_plot <- as.data.frame(connectome[,c("source_type","source","ligand","receptor","receiver","receiver_type","target","edgeweight")])
  colnames(connectome_plot) <- c("Sender\ncelltype","Sender\ncell","Ligand",
                                 "Receptor","Receiver\ncell","Receiver\ncelltype","Target\ngene","Weight")
  connectome_lodes <- to_lodes_form(connectome_plot, axes = 1:7, id = "Cohort")
  connectome_lodes$stratum <- as.factor(connectome_lodes$stratum)

  #to distinguish senders from receivers, so coloring doesn't fail in case the same cell is in each (we want to treat them separately)
  connectome_lodes$stratum <- apply(connectome_lodes,1,function(x) {
    x[4] <- ifelse(x[3]=="Sender\ncell",paste0(x[4],"_send"),x[4])
  })

  if(is.null(cell.cols)) {
    ct <- unique(c(connectome$source_type,connectome$receiver_type))
    cell.cols <- pal_igv()(length(ct))
    names(cell.cols) <- ct
  }
  sender_fill <- cell.cols
  unique_genes <- unique(c(connectome_plot$Ligand,
                           connectome_plot$Receptor,
                           connectome_plot$`Target
gene`))
  gene_fill <- hue_pal()(length(unique_genes))
  names(gene_fill) <- unique_genes
  cells <- as.character(unique(connectome_lodes[connectome_lodes$x %in% c("Receiver\ncell","Sender\ncell"),"stratum"]))
  senders <- as.character(unique(connectome_lodes[connectome_lodes$x %in% c("Sender\ncell"),"stratum"]))
  receivers <- as.character(unique(connectome_lodes[connectome_lodes$x %in% c("Receiver\ncell"),"stratum"]))

  senders_weights <- aggregate(connectome_lodes[connectome_lodes$stratum %in% senders,"Weight"],by = list(connectome_lodes[connectome_lodes$stratum %in% senders,"stratum"]), sum) %>% arrange(x)
  senders_fill <- colorRampPalette(c("grey80","yellow","red"))(nrow(senders_weights))
  names(senders_fill) <- senders_weights$Group.1

  receivers_weights <- aggregate(connectome_lodes[connectome_lodes$stratum %in% receivers,"Weight"],by = list(connectome_lodes[connectome_lodes$stratum %in% receivers,"stratum"]), sum) %>% arrange(x)
  receivers_fill <- colorRampPalette(c("grey80","yellow","red"))(nrow(receivers_weights))
  names(receivers_fill) <- receivers_weights$Group.1

  cell_fill <- c(receivers_fill,senders_fill)

  connectome_lodes$stratum <- as.character(connectome_lodes$stratum)
  ##attempt to make all cells equal. But I'm not sure you can (easily) or should.
  ##Because it's not an even mix and match between senders and receivers.
  ##You can scale all the senders to be the same, but that doesn't guarantee that the all the receivers will be the same, or vice versa.
  # if(cells.equal){
  #   #calculate number of lodes per cell
  #   connectome_cells <- connectome_lodes[connectome_lodes$stratum %in% cells,]
  #   connectome_cells$stratum <- as.character(connectome_cells$stratum)
  #   lpc <- as.data.frame(table(connectome_cells$stratum))
  #   lpc$Freq <- (1/lpc$Freq)
  #   wpc <- unique(connectome_cells)[,c("Cohort","stratum")]
  #   wpc$Weight <- mapvalues(wpc$stratum, from = lpc$Var1, to = lpc$Freq)
  #   connectome_lodes$Weight <- mapvalues(connectome_lodes$Cohort, from = wpc$Cohort, to = wpc$Weight)
  # }

  # connectome_lodes$Weight <- 1

  if(optimize.flows) {
    message("Preparing data for networkD3 . . . ")
    send_to_cell <- merge(connectome_lodes[connectome_lodes$x=="Sender\ncelltype",c("Cohort","x","stratum")],
                          connectome_lodes[connectome_lodes$x=="Sender\ncell",c("Cohort","x","stratum")],
                          by = "Cohort")
    send_to_cell$stratum.x <- paste(send_to_cell$stratum.x,"send",sep = "_")

    cell_to_ligand <- merge(connectome_lodes[connectome_lodes$x=="Sender\ncell",c("Cohort","x","stratum")],
                            connectome_lodes[connectome_lodes$x=="Ligand",c("Cohort","x","stratum")],
                            by = "Cohort")
    cell_to_ligand$stratum.y <- paste(cell_to_ligand$stratum.y,"ligand", sep = "_")

    ligand_to_receptor <- merge(connectome_lodes[connectome_lodes$x=="Ligand",c("Cohort","x","stratum")],
                                connectome_lodes[connectome_lodes$x=="Receptor",c("Cohort","x","stratum")],
                                by = "Cohort")
    ligand_to_receptor$stratum.x <- paste(ligand_to_receptor$stratum.x,"ligand", sep = "_")
    ligand_to_receptor$stratum.y <- paste(ligand_to_receptor$stratum.y,"receptor", sep = "_")

    receptor_to_receiver <- merge(connectome_lodes[connectome_lodes$x=="Receptor",c("Cohort","x","stratum")],
                                  connectome_lodes[connectome_lodes$x=="Receiver\ncell",c("Cohort","x","stratum")],
                                  by = "Cohort")
    receptor_to_receiver$stratum.x <- paste(receptor_to_receiver$stratum.x,"receptor", sep = "_")

    receiver_to_type <- merge(connectome_lodes[connectome_lodes$x=="Receiver\ncell",c("Cohort","x","stratum")],
                              connectome_lodes[connectome_lodes$x=="Receiver\ncelltype",c("Cohort","x","stratum")],
                              by = "Cohort")
    receiver_to_type$stratum.y <- paste(receiver_to_type$stratum.y,"receiver", sep = "_")

    type_to_target <- merge(connectome_lodes[connectome_lodes$x=="Receiver\ncelltype",c("Cohort","x","stratum")],
                            connectome_lodes[connectome_lodes$x=="Target\ngene",c("Cohort","x","stratum")],
                            by = "Cohort")
    type_to_target$stratum.x <- paste(type_to_target$stratum.x,"receiver", sep = "_")
    type_to_target$stratum.y <- paste(type_to_target$stratum.y,"target", sep = "_")

    links <- rbind(data.frame(source=send_to_cell$stratum.x,target=send_to_cell$stratum.y,value=1),
                   data.frame(source=cell_to_ligand$stratum.x,target=cell_to_ligand$stratum.y,value=1),
                   data.frame(source=ligand_to_receptor$stratum.x,target=ligand_to_receptor$stratum.y,value=1),
                   data.frame(source=receptor_to_receiver$stratum.x,target=receptor_to_receiver$stratum.y,value=1),
                   data.frame(source=receiver_to_type$stratum.x,target=receiver_to_type$stratum.y,value=1),
                   data.frame(source=type_to_target$stratum.x,target=type_to_target$stratum.y,value=1))

    node.list <- unique(c(links$source,links$target))
    nodes <- data.frame(node = 1:length(unique(node.list)), name = unique(node.list))
    nodes$node <- nodes$node-1

    links$source <- as.numeric(mapvalues(links$source, from = nodes$name, to = nodes$node, warn_missing = F))
    links$target <- as.numeric(mapvalues(links$target, from = nodes$name, to = nodes$node, warn_missing = F))

    p1 <- sankeyNetwork(Links = links, Nodes = nodes,
                        Source = "source", Target = "target",
                        Value = "value", NodeID = "name",
                        fontSize = 12, nodeWidth = 30)

    customJS <- 'function() { console.log(this.sankey.nodes().map(d => [d.name, d.x, d.y])); }'
    p2 <- htmlwidgets::onRender(p1, customJS)
    saveNetwork(p2, "~/Downloads/sankey_networkD3.html")
    browseURL('file:///Users/aaronwilk/Downloads/sankey_networkD3.html')

    message("Please browse the code for the opened Sankey diagram (control-command-C)
            Navigate to the console tab and copy the contents of this tab
            Type 'yes' below when this has been completed.\n
            If you are prompted more than once, try copying again, you probably copied without the line breaks")

    status <- "no"
    while(status!="yes") {
      status <- readline(prompt = "Have node positions been copied?: ")
      nodePositions <- read_clip_tbl()
      nodePositions <- data.frame(x=nodePositions[!grepl("^\\[",nodePositions[,1]),])
      status <- ifelse(nrow(nodePositions)>1,"yes","no")
    }

    # nodePositions <- read.csv("~/Downloads/nodePositions.csv")

    colnames(nodePositions) <- "x_axis"
    nodePositions$x_axis <- sub("\\].*", "", sub(".*\\[", "", nodePositions$x_axis))

    nodePositions <- as.data.frame(t(as.data.frame(str_split(nodePositions$x_axis,","))))
    rownames(nodePositions) <- NULL
    colnames(nodePositions) <- c("stratum_full","x_axis","y")

    nodePositions <- nodePositions[1:(nrow(nodePositions)-2),]
    trim.leading <- function (x)  sub("^\\s+", "", x)

    nodePositions$x_axis <- as.numeric(trim.leading(nodePositions$x_axis))
    nodePositions$y <- as.numeric(trim.leading(nodePositions$y))
    nodePositions$stratum_full <- gsub("'","",nodePositions$stratum_full)

    nodePositions$x_axis <- scales::rescale(nodePositions$x_axis, newrange = c(0,10*max(nodePositions$y)))
    nodePositions$num_rank <- nodePositions$x_axis+nodePositions$y

    nodePositions %<>%
      mutate(rank = order(order(num_rank, decreasing = F)))
    new.levels <- nodePositions %>% arrange(num_rank) %>% pull(stratum_full)

    connectome_lodes_merge <- connectome_lodes

    connectome_lodes_merge$stratum <- ifelse(connectome_lodes_merge$x=="Sender\ncelltype",
                                             paste(connectome_lodes_merge$stratum,"send",sep = "_"),
                                             connectome_lodes_merge$stratum)
    connectome_lodes_merge$stratum <- ifelse(connectome_lodes_merge$x=="Ligand",
                                             paste(connectome_lodes_merge$stratum,"ligand",sep = "_"),
                                             connectome_lodes_merge$stratum)
    connectome_lodes_merge$stratum <- ifelse(connectome_lodes_merge$x=="Receptor",
                                             paste(connectome_lodes_merge$stratum,"receptor",sep = "_"),
                                             connectome_lodes_merge$stratum)
    connectome_lodes_merge$stratum <- ifelse(connectome_lodes_merge$x=="Receiver\ncelltype",
                                             paste(connectome_lodes_merge$stratum,"receiver",sep = "_"),
                                             connectome_lodes_merge$stratum)
    connectome_lodes_merge$stratum <- ifelse(connectome_lodes_merge$x=="Target\ngene",
                                             paste(connectome_lodes_merge$stratum,"target",sep = "_"),
                                             connectome_lodes_merge$stratum)

    connectome_lodes_merge$stratum <- factor(connectome_lodes_merge$stratum, levels = new.levels)

    send_fill <- sender_fill
    names(send_fill) <- paste(names(send_fill),"send", sep = "_")
    receiver_fill <- sender_fill
    names(receiver_fill) <- paste(names(receiver_fill),"receiver", sep = "_")
    ligand_fill <- gene_fill
    names(ligand_fill) <- paste(names(ligand_fill),"ligand", sep = "_")
    receptor_fill <- gene_fill
    names(receptor_fill) <- paste(names(receptor_fill),"receptor", sep = "_")
    target_fill <- gene_fill
    names(target_fill) <- paste(names(target_fill),"target", sep = "_")

    ggplot(connectome_lodes_merge, aes(x = x, stratum = stratum, alluvium = Cohort,
                                       y = Weight, label = stratum)) +
      geom_flow(lode.guidance = "frontback",color = "darkgray",cement.alluvia = T) +
      stat_stratum(aes(fill = stratum), width = 1/3) + scale_fill_manual(values = c(send_fill,
                                                                                    receiver_fill,
                                                                                    ligand_fill,
                                                                                    receptor_fill,
                                                                                    target_fill,
                                                                                    cell_fill)) +
      ggfittext::geom_fit_text(aes(label = ifelse(as.numeric(x) %in% c(1,3,4,6,7), sub("_[^_]+$", "", stratum), NA)),
                               stat = "stratum", width = 1/3, min.size = 3) +
      theme_cowplot() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + labs(x = NULL) + NoLegend()

  }

  else {
    ggplot(connectome_lodes, aes(x = x, stratum = stratum, alluvium = Cohort,
                                 y = Weight, label = stratum)) +
      geom_flow(lode.guidance = "frontback",color = "darkgray",cement.alluvia = T) +
      stat_stratum(aes(fill = stratum), width = 1/3) + scale_fill_manual(values = c(sender_fill,gene_fill,cell_fill)) +
      ggfittext::geom_fit_text(aes(label = ifelse(as.numeric(x) %in% c(1,3,4,6,7), as.character(stratum), NA)),
                               stat = "stratum", width = 1/3, min.size = 3) +
      theme_cowplot() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + labs(x = NULL) + NoLegend()
  }

}




#' Title
#'
#' @param interaction_graph
#' @param seu
#' @param name
#' @param row.k
#' @param col.k
#' @param cell.type.calls
#' @param cell.cols
#'
#' @return
#' @export
#'
#' @examples
BuildSingleHeatmap <- function(interaction_graph, seu, name = "var_results",
                               row.k = NULL, col.k = NULL, cell.type.calls = "celltype.l2",
                               cell.cols = NULL) {
  interaction_graph <- as.matrix(interaction_graph)
  #cluster rows
  message("Clustering rows")
  input = interaction_graph
  d <- dist(input, method = "euclidean")
  hc1 <- hclust(d, method = "ward.D")
  if(!is.null(row.k)) {
    nclust = row.k
  }
  else {
    test <- fviz_nbclust(input, FUN = hcut, method = "wss")
    nclust = elbow(test$data %>% mutate_all(as.numeric))$clusters_selected
  }
  row_grp <- as.character(cutree(hc1, k = nclust))
  names(row_grp) <- rownames(interaction_graph)

  message("Clustering columns")
  #cluster columns
  input = t(interaction_graph)
  d <- dist(input, method = "euclidean")
  hc1 <- hclust(d, method = "ward.D")
  if(!is.null(col.k)) {
    nclust = col.k
  }
  else {
    test <- fviz_nbclust(input, FUN = hcut, method = "wss")
    nclust = elbow(test$data %>% mutate_all(as.numeric))$clusters_selected
  }
  col_grp <- as.character(cutree(hc1, k = nclust))
  names(col_grp) <- colnames(interaction_graph)

  # nbct <- as.data.frame(table(seu@meta.data[,cell.type.calls],seu$bins))
  # colnames(nbct) <- c("celltype","nb","freq")
  # nbct.f <- nbct %>% group_by(nb) %>% top_n(1,freq) %>% sample_n(1)

  message("Building annotation")
  nbct.f <- data.frame(nb=colnames(seu),celltype=seu@meta.data[,cell.type.calls]) %>%
    dplyr::filter(nb %in% rownames(interaction_graph))


  row.colors <- pal_igv()(length(unique(row_grp)))
  names(row.colors) <- unique(row_grp)

  col.colors <- pal_igv()(length(unique(col_grp)))
  names(col.colors) <- unique(col_grp)

  rowanno <- nbct.f[match(rownames(interaction_graph),as.character(nbct.f$nb)),]
  colanno <- nbct.f[match(colnames(interaction_graph),as.character(nbct.f$nb)),]
  if(is.null(cell.cols)) {
    cell.cols <- pal_igv()(length(unique(seu@meta.data[,cell.type.calls])))
    names(cell.cols) <- unique(seu@meta.data[,cell.type.calls])
  }
  cell.cols_use <- cell.cols[names(cell.cols) %in% nbct.f$celltype]

  longer = resample(which(c(length(col.colors),length(row.colors))==max(length(col.colors),length(row.colors))),1)
  if(longer==1) {
    la = rowAnnotation("SCluster" = row_grp,
                       col = list("SCluster" = row.colors),
                       show_legend = c(F))

    ta = columnAnnotation("RCluster" = col_grp,
                          col = list("RCluster" = col.colors),
                          annotation_name_side = "left",
                          show_legend = c(T),
                          annotation_legend_param = list("RCluster" = list(labels = names(col.colors),
                                                                           legend_gp = gpar(fill = col.colors),
                                                                           title = "Cluster #")))
  }
  if(longer==2) {
    la = rowAnnotation("SCluster" = row_grp,
                       col = list("SCluster" = row.colors),
                       annotation_legend_param = list("SCluster" = list(labels = names(row.colors),
                                                                        legend_gp = gpar(fill = row.colors),
                                                                        title = "Cluster #")))

    ta = columnAnnotation("RCluster" = col_grp,
                          col = list("RCluster" = col.colors),
                          annotation_name_side = "left",
                          show_legend = c(F))
  }


  #
  max_value = quantile(interaction_graph,0.99)
  col_fun.use = colorRamp2(c(0,max_value), c("white","#08306B"))
  #   ct_lgd = Legend(labels = names(cell.cols), legend_gp = gpar(fill = cell.cols),
  #                title = "Cell type", ncol=2)
  #
  #   cluster_lgd = Legend(labels = names(cl_cols.use), legend_gp = gpar(fill = cl_cols.use),
  #                        title = "Cluster #")


  # hm <- Heatmap(interaction_graph, show_row_names = F, show_column_names = F,
  #              name = "SD", top_annotation = ta, left_annotation = la, col = col_fun.use)
  #
  # hm <- draw(hm)
  #
  # ro <- row_order(hm)
  # co <- column_order(hm)
  #
  # time_results <- time_results[ro,co]

  message("Building Heatmap")
  hm <- Heatmap(interaction_graph, show_row_names = F, show_column_names = F,
          name = "Interaction\nscore", top_annotation = ta, left_annotation = la, col = col_fun.use,
          row_split = rowanno$celltype, column_split = colanno$celltype,
          clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D",
          clustering_distance_columns = "euclidean", clustering_method_columns = "ward.D")

  return(list(SCluster=row_grp,RCluster=col_grp,hm=hm))

}


#' Title
#'
#' @param interaction_graphs
#' @param seu
#' @param name
#' @param row.k
#' @param col.k
#' @param cell.type.calls
#'
#' @return
#' @export
#'
#' @examples
BuildGranHeatmap <- function(interaction_graphs, seu, name = "var_results",
                             row.k = NULL, col.k = NULL, cell.type.calls = "celltype.l2") {
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  ig_rr <- interaction_graphs[6:9]
  ig_ll <- interaction_graphs[1:5]

  n <- length(ig_rr)
  rc <- dim(ig_rr[[1]])
  ar_rr <- array(unlist(ig_rr), c(rc, n))

  n <- length(ig_ll)
  rc <- dim(ig_ll[[1]])
  ar_ll <- array(unlist(ig_ll), c(rc, n))

  rr_mean <- apply(ar_rr, c(1, 2), mean)
  ll_mean <- apply(ar_ll, c(1, 2), mean)

  num <- rr_mean-ll_mean

  rr_sd <- apply(ar_rr, c(1, 2), sd)
  ll_sd <- apply(ar_ll, c(1, 2), sd)

  rr_sd <- (rr_sd^2)/4
  ll_sd <- (ll_sd^2)/5

  denom <- sqrt(rr_sd+ll_sd)

  var_results <- num/denom

  rownames(var_results) <- rownames(interaction_graphs[[1]])
  colnames(var_results) <- colnames(interaction_graphs[[1]])


  #cluster rows
  input = var_results
  d <- dist(input, method = "euclidean")
  hc1 <- hclust(d, method = "ward")
  test <- fviz_nbclust(input, FUN = hcut, method = "wss")
  nclust = elbow(test$data %>% mutate_all(as.numeric))$clusters_selected
  if(!is.null(row.k)) {
    nclust = row.k
  }
  row_grp <- as.character(cutree(hc1, k = nclust))
  names(row_grp) <- rownames(var_results)

  #cluster columns
  input = t(var_results)
  d <- dist(input, method = "euclidean")
  hc1 <- hclust(d, method = "ward")
  test <- fviz_nbclust(input, FUN = hcut, method = "wss")
  nclust = elbow(test$data %>% mutate_all(as.numeric))$clusters_selected
  if(!is.null(col.k)) {
    nclust = col.k
  }
  col_grp <- as.character(cutree(hc1, k = nclust))
  names(col_grp) <- colnames(var_results)

  nbct <- as.data.frame(table(seu@meta.data[,cell.type.calls],seu$bins))
  colnames(nbct) <- c("celltype","nb","freq")
  nbct.f <- nbct %>% group_by(nb) %>% top_n(1,freq) %>% sample_n(1)

  row.colors <- pal_d3()(length(unique(row_grp)))
  names(row.colors) <- unique(row_grp)

  col.colors <- pal_d3()(length(unique(col_grp)))
  names(col.colors) <- unique(col_grp)

  rowanno <- nbct.f[match(rownames(var_results),as.character(nbct.f$nb)),]
  colanno <- nbct.f[match(colnames(var_results),as.character(nbct.f$nb)),]

  rowanno$celltype <- factor(as.character(rowanno$celltype),
                             levels = unique(c(as.character(rowanno$celltype),
                                               as.character(colanno$celltype))))
  colanno$celltype <- factor(as.character(colanno$celltype),
                             levels = unique(c(as.character(rowanno$celltype),
                                               as.character(colanno$celltype))))

  cell.cols_use <- cell.cols[names(cell.cols) %in% levels(rowanno$celltype)]

  longer = resample(which(c(length(col.colors),length(row.colors))==max(length(col.colors),length(row.colors))),1)
  if(longer==1) {
    la = rowAnnotation("Sender cell" = rowanno$celltype, "SCluster" = row_grp,
                       col = list("Sender cell" = cell.cols_use,
                                  "SCluster" = row.colors),
                       annotation_legend_param = list("Sender cell" = list(labels = names(cell.cols_use),
                                                                           legend_gp = gpar(fill = cell.cols_use),
                                                                           title = "Cell type", ncol=1)),
                       show_legend = c(T,F))

    ta = columnAnnotation("Receiver cell" = colanno$celltype, "RCluster" = col_grp,
                          col = list("Receiver cell" = cell.cols_use,
                                     "RCluster" = col.colors),
                          annotation_name_side = "left",
                          show_legend = c(F,T),
                          annotation_legend_param = list("RCluster" = list(labels = names(col.colors),
                                                                           legend_gp = gpar(fill = col.colors),
                                                                           title = "Cluster #")))
  }
  if(longer==2) {
    la = rowAnnotation("Sender cell" = rowanno$celltype, "SCluster" = row_grp,
                       col = list("Sender cell" = cell.cols_use,
                                  "SCluster" = row.colors),
                       annotation_legend_param = list("Sender cell" = list(labels = names(cell.cols_use),
                                                                           legend_gp = gpar(fill = cell.cols_use),
                                                                           title = "Cell type", ncol=1),
                                                      "SCluster" = list(labels = names(row.colors),
                                                                        legend_gp = gpar(fill = row.colors),
                                                                        title = "Cluster #")))

    ta = columnAnnotation("Receiver cell" = colanno$celltype, "RCluster" = col_grp,
                          col = list("Receiver cell" = cell.cols_use,
                                     "RCluster" = col.colors),
                          annotation_name_side = "left",
                          show_legend = c(F,F))
  }


  #
  min_value = quantile(var_results,0.01)
  max_value = quantile(var_results,0.99)
  col_fun.use = colorRamp2(c(min_value,0,max_value), c("indianred","white","blue"))
  #   ct_lgd = Legend(labels = names(cell.cols), legend_gp = gpar(fill = cell.cols),
  #                title = "Cell type", ncol=2)
  #
  #   cluster_lgd = Legend(labels = names(cl_cols.use), legend_gp = gpar(fill = cl_cols.use),
  #                        title = "Cluster #")


  # hm <- Heatmap(var_results, show_row_names = F, show_column_names = F,
  #              name = "SD", top_annotation = ta, left_annotation = la, col = col_fun.use)
  #
  # hm <- draw(hm)
  #
  # ro <- row_order(hm)
  # co <- column_order(hm)
  #
  # time_results <- time_results[ro,co]
  # time_results <- translateTimes(time_results)

  hc = hclust(dist(t(var_results), method="euclidean"), method = "ward")
  hc1 = hclust(dist((var_results), method="euclidean"), method = "ward")

  hm <- Heatmap(var_results, show_row_names = F, show_column_names = F,
                cluster_columns = hc, column_dend_reorder = F,
                cluster_rows = hc1, row_dend_reorder = F,
                name = "T-statistic", top_annotation = ta, left_annotation = la, col = col_fun.use)

  # thm <- Heatmap(time_results,col = time.cols, name = "Timepoint",
  # cluster_columns = hc, column_dend_reorder = F, show_column_dend = T,
  # show_row_names = F, show_column_names = F)

  return(list(SCluster=row_grp,RCluster=col_grp,var_hm=hm))
}



#' Title
#'
#' @param connectome
#' @param optimize.flows
#'
#' @return
#' @export
#'
#' @examples
PlotAlluvium_nocells <- function(connectome, optimize.flows = T) {
  connectome_plot <- connectome[,c("source_type","ligand","receptor","receiver_type","target","weight")]

  connectome_plot$unique <- paste(connectome_plot$source_type,
                                  connectome_plot$ligand,
                                  connectome_plot$receptor,
                                  connectome_plot$receiver_type,
                                  connectome_plot$target,
                                  sep = "_")

  connectome_plot <- aggregate(connectome_plot$weight, by = list(connectome_plot$unique), FUN = sum)
  colnames(connectome_plot) <- c("unique","Weight")
  connectome_plot %<>% separate(unique, into = c("source_type","ligand","receptor","receiver_type","target"), sep = "_")

  colnames(connectome_plot) <- c("Sender\ncelltype","Ligand",
                                 "Receptor","Receiver\ncelltype","Target\ngene","Weight")
  connectome_lodes <- to_lodes_form(connectome_plot, axes = 1:5, id = "Cohort")
  connectome_lodes$stratum <- as.factor(connectome_lodes$stratum)

  sender_fill <- cell.cols
  unique_genes <- unique(c(connectome_plot$Ligand,
                           connectome_plot$Receptor,
                           connectome_plot$`Target
                           gene`))
  gene_fill <- hue_pal()(length(unique_genes))
  names(gene_fill) <- unique_genes

  connectome_lodes$stratum <- as.character(connectome_lodes$stratum)


  if(optimize.flows) {
    message("Preparing data for networkD3 . . . ")
    send_to_ligand <- merge(connectome_lodes[connectome_lodes$x=="Sender\ncelltype",c("Cohort","x","stratum","Weight")],
                            connectome_lodes[connectome_lodes$x=="Ligand",c("Cohort","x","stratum")],
                            by = "Cohort")
    send_to_ligand$stratum.x <- paste(send_to_ligand$stratum.x,"send",sep = "_")
    send_to_ligand$stratum.y <- paste(send_to_ligand$stratum.y,"send",sep = "_")

    ligand_to_receptor <- merge(connectome_lodes[connectome_lodes$x=="Ligand",c("Cohort","x","stratum","Weight")],
                                connectome_lodes[connectome_lodes$x=="Receptor",c("Cohort","x","stratum")],
                                by = "Cohort")
    ligand_to_receptor$stratum.x <- paste(ligand_to_receptor$stratum.x,"ligand", sep = "_")
    ligand_to_receptor$stratum.y <- paste(ligand_to_receptor$stratum.y,"receptor", sep = "_")

    receptor_to_type <- merge(connectome_lodes[connectome_lodes$x=="Receptor",c("Cohort","x","stratum","Weight")],
                              connectome_lodes[connectome_lodes$x=="Receiver\ncelltype",c("Cohort","x","stratum")],
                              by = "Cohort")
    receptor_to_type$stratum.x <- paste(receptor_to_type$stratum.x,"receiver", sep = "_")
    receptor_to_type$stratum.y <- paste(receptor_to_type$stratum.y,"receiver", sep = "_")

    type_to_target <- merge(connectome_lodes[connectome_lodes$x=="Receiver\ncelltype",c("Cohort","x","stratum","Weight")],
                            connectome_lodes[connectome_lodes$x=="Target\ngene",c("Cohort","x","stratum")],
                            by = "Cohort")
    type_to_target$stratum.x <- paste(type_to_target$stratum.x,"receiver", sep = "_")
    type_to_target$stratum.y <- paste(type_to_target$stratum.y,"target", sep = "_")

    links <- rbind(data.frame(source=send_to_ligand$stratum.x,target=send_to_ligand$stratum.y,value=send_to_ligand$Weight),
                   data.frame(source=ligand_to_receptor$stratum.x,target=ligand_to_receptor$stratum.y,value=ligand_to_receptor$Weight),
                   data.frame(source=receptor_to_type$stratum.x,target=receptor_to_type$stratum.y,value=receptor_to_type$Weight),
                   data.frame(source=type_to_target$stratum.x,target=type_to_target$stratum.y,value=type_to_target$Weight))

    node.list <- unique(c(links$source,links$target))
    nodes <- data.frame(node = 1:length(unique(node.list)), name = unique(node.list))
    nodes$node <- nodes$node-1

    links$source <- as.numeric(plyr::mapvalues(links$source, from = nodes$name, to = nodes$node, warn_missing = F))
    links$target <- as.numeric(plyr::mapvalues(links$target, from = nodes$name, to = nodes$node, warn_missing = F))

    p1 <- sankeyNetwork(Links = links, Nodes = nodes,
                        Source = "source", Target = "target",
                        Value = "value", NodeID = "name",
                        fontSize = 12, nodeWidth = 30)

    customJS <- 'function() { console.log(this.sankey.nodes().map(d => [d.name, d.x, d.y])); }'
    p2 <- htmlwidgets::onRender(p1, customJS)
    saveNetwork(p2, "~/Downloads/sankey_networkD3.html")
    browseURL('file:///Users/aaronwilk/Downloads/sankey_networkD3.html')

    message("Please browse the code for the opened Sankey diagram (control-command-C)
            Navigate to the console tab and copy the contents of this tab
            Type 'yes' below when this has been completed.\n
            If you are prompted more than once, try copying again, you probably copied without the line breaks")

    status <- "no"
    while(status!="yes") {
      status <- readline(prompt = "Have node positions been copied?: ")
      nodePositions <- read_clip_tbl()
      nodePositions <- data.frame(x=nodePositions[!grepl("^\\[",nodePositions[,1]),])
      status <- ifelse(nrow(nodePositions)>1,"yes","no")
    }

    # nodePositions <- read.csv("~/Downloads/nodePositions.csv")

    colnames(nodePositions) <- "x_axis"
    nodePositions$x_axis <- sub("\\].*", "", sub(".*\\[", "", nodePositions$x_axis))

    nodePositions <- as.data.frame(t(as.data.frame(str_split(nodePositions$x_axis,","))))
    rownames(nodePositions) <- NULL
    colnames(nodePositions) <- c("stratum_full","x_axis","y")

    nodePositions <- nodePositions[1:(nrow(nodePositions)-2),]
    trim.leading <- function (x)  sub("^\\s+", "", x)

    nodePositions$x_axis <- as.numeric(trim.leading(nodePositions$x_axis))
    nodePositions$y <- as.numeric(trim.leading(nodePositions$y))
    nodePositions$stratum_full <- gsub("'","",nodePositions$stratum_full)

    nodePositions$x_axis <- scales::rescale(nodePositions$x_axis, newrange = c(0,10*max(nodePositions$y)))
    nodePositions$num_rank <- nodePositions$x_axis+nodePositions$y

    nodePositions %<>%
      mutate(rank = order(order(num_rank, decreasing = F)))
    new.levels <- nodePositions %>% arrange(num_rank) %>% pull(stratum_full)

    connectome_lodes_merge <- connectome_lodes

    connectome_lodes_merge$stratum <- ifelse(connectome_lodes_merge$x=="Sender\ncelltype",
                                             paste(connectome_lodes_merge$stratum,"send",sep = "_"),
                                             connectome_lodes_merge$stratum)
    connectome_lodes_merge$stratum <- ifelse(connectome_lodes_merge$x=="Ligand",
                                             paste(connectome_lodes_merge$stratum,"ligand",sep = "_"),
                                             connectome_lodes_merge$stratum)
    connectome_lodes_merge$stratum <- ifelse(connectome_lodes_merge$x=="Receptor",
                                             paste(connectome_lodes_merge$stratum,"receptor",sep = "_"),
                                             connectome_lodes_merge$stratum)
    connectome_lodes_merge$stratum <- ifelse(connectome_lodes_merge$x=="Receiver\ncelltype",
                                             paste(connectome_lodes_merge$stratum,"receiver",sep = "_"),
                                             connectome_lodes_merge$stratum)
    connectome_lodes_merge$stratum <- ifelse(connectome_lodes_merge$x=="Target\ngene",
                                             paste(connectome_lodes_merge$stratum,"target",sep = "_"),
                                             connectome_lodes_merge$stratum)

    connectome_lodes_merge$stratum <- factor(connectome_lodes_merge$stratum, levels = new.levels)

    send_fill <- sender_fill
    names(send_fill) <- paste(names(send_fill),"send", sep = "_")
    receiver_fill <- sender_fill
    names(receiver_fill) <- paste(names(receiver_fill),"receiver", sep = "_")
    ligand_fill <- gene_fill
    names(ligand_fill) <- paste(names(ligand_fill),"ligand", sep = "_")
    receptor_fill <- gene_fill
    names(receptor_fill) <- paste(names(receptor_fill),"receptor", sep = "_")
    target_fill <- gene_fill
    names(target_fill) <- paste(names(target_fill),"target", sep = "_")

    ggplot(connectome_lodes_merge, aes(x = x, stratum = stratum, alluvium = Cohort,
                                       y = Weight, label = stratum)) +
      geom_flow(stat = "alluvium", lode.guidance = "frontback",color = "darkgray",cement.alluvia = T) +
      geom_stratum(aes(fill = stratum), width = 1/3) + scale_fill_manual(values = c(send_fill,
                                                                                    receiver_fill,
                                                                                    ligand_fill,
                                                                                    receptor_fill,
                                                                                    target_fill)) +
      ggfittext::geom_fit_text(aes(label = ifelse(as.numeric(x) %in% c(1:5), sub("_[^_]+$", "", stratum), NA)),
                               stat = "stratum", width = 1/3, min.size = 3) +
      theme_cowplot() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + labs(x = NULL) + NoLegend()

  }

  else {
    ggplot(connectome_lodes, aes(x = x, stratum = stratum, alluvium = Cohort,
                                 y = Weight, label = stratum)) +
      geom_flow(lode.guidance = "frontback",color = "darkgray",cement.alluvia = T) +
      stat_stratum(aes(fill = stratum), width = 1/3) + scale_fill_manual(values = c(sender_fill,gene_fill,cell_fill)) +
      ggfittext::geom_fit_text(aes(label = ifelse(as.numeric(x) %in% c(1,3,4,6,7), as.character(stratum), NA)),
                               stat = "stratum", width = 1/3, min.size = 3) +
      theme_cowplot() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + labs(x = NULL) + NoLegend()
  }

}



#' Title
#'
#' @param ccim
#' @param seu
#' @param features
#' @param type_plot
#' @param ...
#'
#' @return
#' @import magrittr
#' @export
#'
#' @examples
CCIMFeaturePlot <- function(ccim, seu,
                        features = NULL,
                        type_plot = c("sender","receiver"), ...) {
  #find cells of interest in the original object
  type_plot <- match.arg(arg = type_plot, several.ok = F)
  p <- FeaturePlot(seu, features = features, combine = F)
  if(length(features)>1) {
    pdata <- as.data.frame(t(do.call(rbind,lapply(p,function(x) {x$data[,4]}))))
    colnames(pdata) <- features
    rownames(pdata) <- colnames(seu)
  }
  else {
    pdata <- as.data.frame(x=p[[1]]$data[,4])
    colnames(pdata) <- features
    rownames(pdata) <- colnames(seu)
  }
  ccim_meta <- merge((ccim@meta.data %>% rownames_to_column("cell")),
                     (pdata %>% rownames_to_column(type_plot)),
                     by = type_plot, all.y = F)
  ccim_meta <- ccim_meta[match(colnames(ccim),ccim_meta$cell),]
  rownames(ccim_meta) <- NULL
  ccim_meta %<>%
    column_to_rownames("cell")
  ccim@meta.data <- ccim_meta
  FeaturePlot(ccim, features = features, ... )
}


#' Plot expression of interaction programs
#'
#' @param seu A seurat object that has interaction program expression information stored in the meta.data slot as the output of `ScoreInteractionPrograms`
#' @param ip The name of the interaction program to plot
#' @param cols Character vector of colors to use for plotting
#' @param order Logical. Order points so highest expressers are on top?
#'
#' @return
#' @import ggplot2
#' @export
#'
#' @examples
IPFeaturePlot <- function(seu, ip, cols = c("grey90","blue","orangered3"), order = T) {
  seu$lig_feature <- seu[["IPligands"]]@data[ip,]
  seu$rec_feature <- seu[["IPreceptors"]]@data[ip,]
  # p <- FeaturePlot(seu,
  #             features = c(paste0("ligands_",ip),
  #                          paste0("receptors_",ip)),
  #             blend = T, combine = F,
  #             cols = cols, order = order)
  p <- FeaturePlot(seu,
                   features = c("lig_feature","rec_feature"),
                   blend = T, combine = F,
                   cols = cols, order = order)
  p[[3]] + NoLegend() + labs(x = NULL, y = NULL, title = NULL) +
    theme(aspect.ratio = 1,axis.ticks = element_blank(),
          axis.text = element_blank())
}


#' Plot Alluvial of ligand-target links
#'
#' @param connectome
#' @param optimize.flows
#' @param colors
#'
#' @return
#' @export
#'
#' @examples
PlotAlluviumLT <- function(connectome, optimize.flows = T, colors = NULL) {
  connectome_plot <- as.data.frame(connectome[,c("ligand","target","mean")])
  colnames(connectome_plot) <- c("Ligand","Target","Weight")
  connectome_lodes <- to_lodes_form(connectome_plot, axes = 1:2, id = "Cohort")
  connectome_lodes$stratum <- as.character(connectome_lodes$stratum)

  if(is.null(colors)) {
    unique_genes <- unique(connectome_lodes$stratum)
    gene_fill <- hue_pal()(length(unique_genes))
    names(gene_fill) <- unique_genes
  }
  else {
    gene_fill <- colors
  }

  if(optimize.flows) {
    message("Preparing data for networkD3 . . . ")

    ligand_to_receptor <- merge(connectome_lodes[connectome_lodes$x=="Ligand",c("Cohort","x","stratum")],
                                connectome_lodes[connectome_lodes$x=="Target",c("Cohort","x","stratum")],
                                by = "Cohort")
    ligand_to_receptor$stratum.x <- paste(ligand_to_receptor$stratum.x,"ligand", sep = "_")
    ligand_to_receptor$stratum.y <- paste(ligand_to_receptor$stratum.y,"target", sep = "_")

    links <- rbind(data.frame(source=ligand_to_receptor$stratum.x,
                              target=ligand_to_receptor$stratum.y,value=1))

    node.list <- unique(c(links$source,links$target))
    nodes <- data.frame(node = 1:length(unique(node.list)), name = unique(node.list))
    nodes$node <- nodes$node-1

    links$source <- as.numeric(mapvalues(links$source, from = nodes$name, to = nodes$node, warn_missing = F))
    links$target <- as.numeric(mapvalues(links$target, from = nodes$name, to = nodes$node, warn_missing = F))

    p1 <- sankeyNetwork(Links = links, Nodes = nodes,
                        Source = "source", Target = "target",
                        Value = "value", NodeID = "name",
                        fontSize = 12, nodeWidth = 30)

    customJS <- 'function() { console.log(this.sankey.nodes().map(d => [d.name, d.x, d.y])); }'
    p2 <- htmlwidgets::onRender(p1, customJS)
    saveNetwork(p2, "~/Downloads/sankey_networkD3.html")
    browseURL('file:///Users/aaronwilk/Downloads/sankey_networkD3.html')

    message("Please browse the code for the opened Sankey diagram (control-command-C)
            Navigate to the console tab and copy the contents of this tab
            Type 'yes' below when this has been completed.\n
            If you are prompted more than once, try copying again, you probably copied without the line breaks")

    status <- "no"
    while(status!="yes") {
      status <- readline(prompt = "Have node positions been copied?: ")
      nodePositions <- read_clip_tbl()
      nodePositions <- data.frame(x=nodePositions[!grepl("^\\[",nodePositions[,1]),])
      status <- ifelse(nrow(nodePositions)>1,"yes","no")
    }

    # nodePositions <- read.csv("~/Downloads/nodePositions.csv")

    colnames(nodePositions) <- "x_axis"
    nodePositions$x_axis <- sub("\\].*", "", sub(".*\\[", "", nodePositions$x_axis))

    nodePositions <- as.data.frame(t(as.data.frame(str_split(nodePositions$x_axis,","))))
    rownames(nodePositions) <- NULL
    colnames(nodePositions) <- c("stratum_full","x_axis","y")

    nodePositions <- nodePositions[1:(nrow(nodePositions)-1),]
    trim.leading <- function (x)  sub("^\\s+", "", x)

    nodePositions$x_axis <- as.numeric(trim.leading(nodePositions$x_axis))
    nodePositions$y <- as.numeric(trim.leading(nodePositions$y))
    nodePositions$stratum_full <- gsub("'","",nodePositions$stratum_full)

    nodePositions$x_axis <- scales::rescale(nodePositions$x_axis, newrange = c(0,10*max(nodePositions$y)))
    nodePositions$num_rank <- nodePositions$x_axis+nodePositions$y

    nodePositions %<>%
      mutate(rank = order(order(num_rank, decreasing = F)))
    new.levels <- nodePositions %>% arrange(num_rank) %>% pull(stratum_full)

    connectome_lodes_merge <- connectome_lodes
    connectome_lodes_merge$stratum <- ifelse(connectome_lodes_merge$x=="Ligand",
                                             paste(connectome_lodes_merge$stratum,"ligand",sep = "_"),
                                             connectome_lodes_merge$stratum)
    connectome_lodes_merge$stratum <- ifelse(connectome_lodes_merge$x=="Target",
                                             paste(connectome_lodes_merge$stratum,"target",sep = "_"),
                                             connectome_lodes_merge$stratum)

    connectome_lodes_merge$stratum <- factor(connectome_lodes_merge$stratum, levels = new.levels)

    ligand_fill <- gene_fill
    names(ligand_fill) <- paste(names(ligand_fill),"ligand", sep = "_")
    target_fill <- gene_fill
    names(target_fill) <- paste(names(target_fill),"target", sep = "_")

    ggplot(connectome_lodes_merge, aes(x = x, stratum = stratum, alluvium = Cohort,
                                       y = Weight, label = stratum)) +
      geom_flow(lode.guidance = "frontback",color = "darkgray",cement.alluvia = T) +
      stat_stratum(aes(fill = stratum), width = 1/3) + scale_fill_manual(values = c(ligand_fill,
                                                                                    target_fill)) +
      ggfittext::geom_fit_text(aes(label = ifelse(as.numeric(x) %in% c(1,2), sub("_[^_]+$", "", stratum), NA)),
                               stat = "stratum", width = 1/3, min.size = 3) +
      theme_cowplot() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + labs(x = NULL) + NoLegend()

  }

  else {
    ggplot(connectome_lodes, aes(x = x, stratum = stratum, alluvium = Cohort,
                                 y = Weight, label = stratum)) +
      geom_flow(lode.guidance = "frontback",color = "darkgray",cement.alluvia = T) +
      stat_stratum(aes(fill = stratum), width = 1/3) + scale_fill_manual(values = c(gene_fill)) +
      ggfittext::geom_fit_text(aes(label = ifelse(as.numeric(x) %in% c(1,2), as.character(stratum), NA)),
                               stat = "stratum", width = 1/3, min.size = 3) +
      theme_cowplot() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + labs(x = NULL) + NoLegend()
  }

}



#' Plot Pearson distributions for ligand activities
#'
#' @param ligand_activities the matrix output of `RankActiveLigands`
#' @param pearson_cutoff determines where to draw a vertical line illustrating the chosen Pearson threshold cutoff (default: 0.075)
#'
#' @return Two plots of all Pearson coefficients (left) and median Pearson coefficient per cell (right)
#' @export
#'
#' @examples
PearsonDist <- function(ligand_activities, pearson_cutoff = 0.075) {
  ral_df <- data.frame(x = as.vector(ligand_activities))
  minx <- min(ral_df$x)
  maxx <- max(ral_df$x)
  p1 <- ggplot(ral_df, aes(x = x)) + geom_density(fill = "steelblue") + theme_cowplot() +
    labs(x = "Pearson coefficient", y = "Density", title = "All Pearson coefficients") +
    xlim(c(minx,maxx)) + geom_vline(xintercept = pearson_cutoff) +
    theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1)
  ral_df <- data.frame(x = as.vector(apply(ligand_activities,2,median)))
  p2 <- ggplot(ral_df, aes(x = x)) + geom_density(fill = "steelblue") + theme_cowplot() +
    labs(x = "Pearson coefficient", y = "Density", title = "Median Pearson coefficient\nper cell") +
    xlim(c(minx,maxx)) + geom_vline(xintercept = pearson_cutoff) +
    theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1)
  p1|p2
}

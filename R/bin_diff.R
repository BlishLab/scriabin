

#' Identify most perturbed bins across >2 samples
#'
#' @param seu A binned Seurat object with bin identities in the "bin" column of meta.data
#' @param interaction_graphs List of summarized interaction graphs built from `BuildPriorInteraction` or `BuiltWeightedInteraction`
#'
#' @return For each sender bin-receiver bin combination, returns a data.frame of Kruskal-Wallis p-value and test statistic
#' @export
#'
#' @examples
PerturbedBins <- function(seu, interaction_graphs = NULL) {
  nbin = length(unique(seu$bins))
  bin_combos <- expand.grid(1:nbin,1:nbin)
  bins <- unique(seu$bins)

  bin_pvals <- pblapply(seq_along(1:nrow(bin_combos)), function(x) {
    sbin <- as.character(bin_combos[x,1])
    rbin <- as.character(bin_combos[x,2])
    scells <- colnames(seu)[seu$bins==sbin]
    rcells <- colnames(seu)[seu$bins==rbin]

    bin_ig <- reshape2::melt(lapply(interaction_graphs, function(x) {
      as.vector(x[rownames(x) %in% scells,colnames(x) %in% rcells])
    }))

    res <- kruskal.test(value~L1, data = bin_ig)
    return(c(res$p.value,
             res$statistic,
             nrow(bin_ig),
             sd(aggregate(bin_ig$value,by = list(bin_ig$L1), FUN=mean)$x)))
  })
  return(bin_pvals)
}

#' Summary plots of perturbed bins
#'
#' @param seu A binned Seurat object with bin identities in the "bin" column of meta.data
#' @param bin_pvals Bin perturbation signifance data, ie. the output of `PerturbedBins`
#' @param cell.type.calls Meta.data slot column corresponding to cell type annotations for summarization
#'
#' @return Returns plots of bin-bin perturbation p-values, cross sample communication standard deviation, and summary by cell type or other meta.data of interest
#' @import ggplot2 dplyr cowplot scales
#' @export
#'
#' @examples
PerturbedBinSummary <- function(seu, bin_pvals, cell.type.calls = "celltype.l2") {
  bin_pvals_p <- unlist(lapply(bin_pvals, function(x) {x[1]}))
  bin_pvals_stat <- unlist(lapply(bin_pvals, function(x) {x[2]}))
  bin_pvals_size <- unlist(lapply(bin_pvals, function(x) {x[3]}))
  bin_pvals_sd <- unlist(lapply(bin_pvals, function(x) {x[4]}))

  ggplot(data.frame(x = bin_pvals_size, y = -log(bin_pvals_p)), aes(x=x,y=y)) +
    geom_bin2d(bins=100) + theme_cowplot() + scale_fill_gradient(trans="log", low = "blue", high = "red")

  ggplot(data.frame(x = bin_pvals_size, y = -log(bin_pvals_p), color = bin_pvals_sd), aes(x=x,y=y,color=color)) +
    geom_point(size=0) + theme_cowplot() + scale_color_gradient(trans="log", low = "blue", high = "red")

  #Correct for number of bins as a function of cell type proportion
  nbct <- as.data.frame(table(as.character(seu@meta.data[,cell.type.calls]),as.character(seu$bins)))
  colnames(nbct) <- c("celltype","nb","freq")
  nbct.f <- nbct %>% group_by(nb) %>% top_n(1,freq) %>% sample_n(1) %>% ungroup() %>% mutate_all(as.character)

  expected_sizes <- as.data.frame(table(seu@meta.data[,cell.type.calls]))
  colnames(expected_sizes) <- c("celltype","celln")
  bin_types <- as.data.frame(table(nbct.f$celltype))
  colnames(bin_types) <- c("celltype","binn")
  expected_sizes <- merge(expected_sizes,bin_types, by = "celltype")
  expected_sizes$corrected <- lm(binn~celln, data = expected_sizes)$fitted.values
  expected_sizes$adjustment_factor <- expected_sizes$corrected/expected_sizes$binn

  pbins_plot <- as.data.frame(table(sig_bin_combos$send_celltype,sig_bin_combos$rec_celltype)) %>% dplyr::filter(Freq>0)
  colnames(pbins_plot) <- c("send_celltype","rec_celltype","Freq")
  pbins_plot$send_adjustment <- scriabin::mapvalues(pbins_plot$send_celltype,
                                                from = expected_sizes$celltype,
                                                to = expected_sizes$adjustment_factor)
  pbins_plot$rec_adjustment <- scriabin::mapvalues(pbins_plot$rec_celltype,
                                               from = expected_sizes$celltype,
                                               to = expected_sizes$adjustment_factor)
  pbins_plot$freq_adjusted <- as.numeric(as.character(pbins_plot$Freq))*
    as.numeric(as.character(pbins_plot$send_adjustment))*
    as.numeric(as.character(pbins_plot$rec_adjustment))

  ggplot(pbins_plot, aes(x = rec_celltype, y = send_celltype)) + geom_point(aes(size=freq_adjusted,color=freq_adjusted)) +
    theme_cowplot() + ggpubr::rotate_x_text() +
    scale_color_gradient(low = "lightgrey",high = "red", limits = c(10,100), oob = scales::squish) +
    labs(x = "Receiving celltype", y = "Sender celltype", color = "Number of\nperturbed\nbin pairs") + guides(size=F)
}

#' Perform post-hoc analysis of perturbed bins
#'
#' @param seu A binned Seurat object with bin identities in the "bin" column of meta.data
#' @param bin_pvals Bin perturbation signifance data, ie. the output of `PerturbedBins`
#' @param interaction_graphs List of summarized interaction graphs built from `BuildPriorInteraction` or `BuiltWeightedInteraction`
#' @param split.by Meta.data column name indicating how data was split for interaction graph generation
#' @param cell.type.calls Meta.data slot column corresponding to cell type annotations for summarization
#' @param kw_p.value Bin-bin combinations with a KW p value above this threshold will be discarded. Default: 0.001
#' @param bin_sd.quantile Bin-bin combinations with a summarized interaction standard deviation below this quantile will be discarded. Ensures that bin-bin combinations displaying both statistical significance and effect size of perturbation are analyzed. Default: 0.9.
#'
#' @return Performs Dunn's Kruskal-Wallis multiple comparison post-hoc test to evaluate which samples within perturbed bins are significantly perturbed
#' @import pbapply FSA dplyr
#' @export
#'
#' @examples
BinPostHoc <- function(seu, bin_pvals, interaction_graphs,
                       split.by = "orig.ident", cell.type.calls = "celltype",
                       kw_p.value = 0.001, bin_sd.quantile = 0) {
  nbin = length(unique(seu$bins))
  bin_combos <- expand.grid(1:nbin,1:nbin)
  bins <- unique(seu$bins)

  bin_pvals_p <- unlist(lapply(bin_pvals, function(x) {x[1]}))
  bin_pvals_stat <- unlist(lapply(bin_pvals, function(x) {x[2]}))
  bin_pvals_size <- unlist(lapply(bin_pvals, function(x) {x[3]}))
  bin_pvals_sd <- unlist(lapply(bin_pvals, function(x) {x[4]}))

  #Perform a post-hoc test to identify which groups differ most in each significant bin
  bin_combos$p_vals <- bin_pvals_p
  bin_combos$kw_stat <- bin_pvals_stat
  bin_combos$size <- bin_pvals_size
  bin_combos$sd <- bin_pvals_sd
  sig_bin_combos <- bin_combos %>%
    dplyr::filter(p_vals<kw_p.value) %>%
    dplyr::filter(sd>quantile(bin_combos$sd,bin_sd.quantile))

  bin_posthoc <- pblapply(seq_along(1:nrow(sig_bin_combos)), function(x) {
    sbin <- bins[sig_bin_combos[x,1]]
    rbin <- bins[sig_bin_combos[x,2]]
    scells <- colnames(seu)[seu$bins==sbin]
    rcells <- colnames(seu)[seu$bins==rbin]

    bin_ig <- lapply(interaction_graphs, function(x) {
      as.vector(x[rownames(x) %in% scells,colnames(x) %in% rcells])
    })
    bin_ig.df <- data.frame(name = as.factor(unlist(lapply(seq_along(1:length(bin_ig)),
                                                           function(x) {rep(names(bin_ig)[x],lapply(bin_ig,length)[[x]])}))),
                            score=unlist(bin_ig))

    p <- dunnTest(score~name, data = bin_ig.df, method = "bh")$res %>%
      separate(Comparison, into = c("tp1","tp2"), sep = " - ") %>%
      dplyr::filter(P.adj<0.05)

    bin_means <- aggregate(bin_ig.df$score, by = list(bin_ig.df$name), FUN=mean)
    # bin mean that is furthest from group mean (mean of all timepoints)
    # bin_means$x <- abs(bin_means$x-mean(bin_means$x))

    # bin mean that is furthest from baseline
    bin_means$x <- abs(bin_means$x-bin_means[1,"x"])
    # bin_means %<>% dplyr::filter(Group.1 %in% unique(c(p$tp1,p$tp2))) %>%
    #   top_n(1,x) %>% pull(Group.1)

    return(as.character(bin_means %>% dplyr::filter(Group.1 %in% unique(c(p$tp1,p$tp2))) %>%
                          top_n(1,x) %>% pull(Group.1)))
  })


  sig_bin_combos$tp <- (bin_posthoc)
  #some of these bin-bin combos will return zero statistically significant results between timepoints. Remove those.
  sig_bin_combos <- sig_bin_combos[sig_bin_combos$tp %in% unique(seu@meta.data[,split.by]),]


  #now let's transfer these bin numbers to cell type calls
  sig_bin_combos$Var1 <- scriabin::mapvalues(sig_bin_combos$Var1, from = 1:length(bins), to = bins,warn_missing = F)
  sig_bin_combos$Var2 <- scriabin::mapvalues(sig_bin_combos$Var2, from = 1:length(bins), to = bins,warn_missing = F)

  nbct <- as.data.frame(table(as.character(seu@meta.data[,cell.type.calls]),as.character(seu$bins)))
  colnames(nbct) <- c("celltype","nb","freq")
  nbct.f <- nbct %>% group_by(nb) %>% top_n(1,freq) %>% sample_n(1) %>% ungroup() %>% mutate_all(as.character)

  sig_bin_combos$send_celltype <- as.character(scriabin::mapvalues(sig_bin_combos$Var1,
                                                               from = nbct.f$nb,
                                                               to = nbct.f$celltype, warn_missing = F))
  sig_bin_combos$rec_celltype <- as.character(scriabin::mapvalues(sig_bin_combos$Var2,
                                                              from = nbct.f$nb, to = nbct.f$celltype,
                                                              warn_missing = F))
  return(sig_bin_combos)
}

#' Assign quality scores to each bin
#'
#' @param seu A binned Seurat object with bin identities in the "bin" column of meta.data
#' @param celltype.calls Name of meta.data column containing the cell type labels to use for connectivity testing. Must be a finer (more granular) cell type label than used for coarse_cell_types in `BinDatasets`.
#' @param split.by Meta.data column name indicating how data was split for interaction graph generation
#'
#' @return A named numeric vector of bin quality scores
#' @import pbapply dplyr
#' @export
#'
#' @examples
ScoreBins <- function(seu, celltype.calls = NULL, split.by = NULL) {
  bins <- as.character(unique(seu$bins))
  scores <- unlist(pblapply(1:length(bins), function(x) {
    meta <- seu@meta.data[seu$bins==bins[x],c(celltype.calls,"bins",split.by)]
    colnames(meta) <- c("celltype","bins","split")
    mp_total <- meta %>% dplyr::count(split)
    mp_count <- meta %>% dplyr::count(split,celltype)
    n_unique <- length(unique(mp_count$celltype))
    mp_count <- mp_count %>%
      group_by(split) %>% dplyr::slice(which.max(n)) %>%
      left_join(.,mp_total,by = "split") %>% dplyr::mutate(prop = n.x/n.y) %>%
      pull(prop)
    score <- (1-sd(mp_count))/n_unique
  }))
  names(scores) <- bins
  return(scores)
}











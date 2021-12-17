


#' Title
#'
#' @param seu
#' @param nsample
#' @param interaction_graphs
#'
#' @return
#' @export
#'
#' @examples
LongPerturbedBins <- function(seu, nsample = 500, interaction_graphs = NULL) {
  nbin = length(unique(seu$bins))
  bin_combos <- expand.grid(1:nbin,1:nbin)
  bins <- unique(seu$bins)
  bin_pvals <- pblapply(seq_along(1:(nbin*nbin)), function(x) {
    sbin <- bins[bin_combos[x,1]]
    rbin <- bins[bin_combos[x,2]]
    scells <- colnames(seu)[seu$bins==sbin]
    rcells <- colnames(seu)[seu$bins==rbin]

    bin_ig <- lapply(ogig, function(x) {
      as.vector(x[rownames(x) %in% scells,colnames(x) %in% rcells])
    })
    bin_ig.df <- data.frame(name = unlist(lapply(seq_along(1:length(bin_ig)),
                                                 function(x) {rep(names(bin_ig)[x],lapply(bin_ig,length)[[x]])})),
                            score=unlist(bin_ig))

    res <- kruskal.test(score~name, data = bin_ig.df)
    return(c(res$p.value,
             res$statistic,
             nrow(bin_ig.df),
             sd(aggregate(bin_ig.df$score,by = list(bin_ig.df$name), FUN=mean)$x)))
  })
  return(bin_pvals)
}

#' Title
#'
#' @param seu
#' @param bin_pvals
#' @param cell.type.calls
#'
#' @return
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

  ggplot(data.frame(x = bin_pvals_size, y = ((bin_pvals_p)), color = bin_pvals_sd), aes(x=x,y=y,color=color)) +
    geom_point(size=0) + theme_cowplot() + scale_color_gradient(trans="log", low = "blue", high = "red")

  #need an appropriate way of correcting for the number of bins. Eg. even though there aren't many B cells there are a ton of B intermediate bins. This will make things look weird.
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
  pbins_plot$send_adjustment <- plyr::mapvalues(pbins_plot$send_celltype,
                                                from = expected_sizes$celltype,
                                                to = expected_sizes$adjustment_factor)
  pbins_plot$rec_adjustment <- plyr::mapvalues(pbins_plot$rec_celltype,
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

#' Title
#'
#' @param seu
#' @param bin_pvals
#' @param interaction_graphs
#' @param split.by
#' @param cell.type.calls
#'
#' @return
#' @export
#'
#' @examples
BinPostHoc <- function(seu, bin_pvals, interaction_graphs,
                       split.by = "time.orig", cell.type.calls = "celltype.l2") {
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
    dplyr::filter(p_vals<0.001) %>%
    dplyr::filter(sd>quantile(bin_combos$sd,0.9))

  bin_posthoc <- pblapply(seq_along(1:nrow(sig_bin_combos)), function(x) {
    sbin <- bins[sig_bin_combos[x,1]]
    rbin <- bins[sig_bin_combos[x,2]]
    scells <- colnames(seu)[seu$bins==sbin]
    rcells <- colnames(seu)[seu$bins==rbin]

    bin_ig <- lapply(ogig, function(x) {
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
  sig_bin_combos <- sig_bin_combos[sig_bin_combos$tp %in% unique(seu$time.orig),]


  #now let's transfer these bin numbers to cell type calls
  sig_bin_combos$Var1 <- plyr::mapvalues(sig_bin_combos$Var1, from = 1:length(bins), to = bins,warn_missing = F)
  sig_bin_combos$Var2 <- plyr::mapvalues(sig_bin_combos$Var2, from = 1:length(bins), to = bins,warn_missing = F)

  nbct <- as.data.frame(table(as.character(seu@meta.data[,cell.type.calls]),as.character(seu$bins)))
  colnames(nbct) <- c("celltype","nb","freq")
  nbct.f <- nbct %>% group_by(nb) %>% top_n(1,freq) %>% sample_n(1) %>% ungroup() %>% mutate_all(as.character)

  sig_bin_combos$send_celltype <- as.character(plyr::mapvalues(sig_bin_combos$Var1,
                                                               from = nbct.f$nb,
                                                               to = nbct.f$celltype, warn_missing = F))
  sig_bin_combos$rec_celltype <- as.character(plyr::mapvalues(sig_bin_combos$Var2,
                                                              from = nbct.f$nb, to = nbct.f$celltype,
                                                              warn_missing = F))
  return(sig_bin_combos)
}












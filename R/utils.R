


#' Load ligand-receptor database
#'
#' @param species character. Name of species from which to load ligand-receptor databases. One of: "human", "mouse", "rat". Default: "human"
#' @param database Name of ligand-receptor database to use. Default: "OmniPath"
#' When species is "human", one of: OmniPath, CellChatDB, CellPhoneDB, Ramilowski2015, Baccin2019, LRdb, Kirouac2010, ICELLNET, iTALK, EMBRACE, HPMR, Guide2Pharma, connectomeDB2020, talklr, CellTalkDB
#' When species is "mouse" or "rat", only "OmniPath" is supported.
#' To pass a custom ligand-receptor database to this function, set database = "custom"
#' @param ligands Character vector of custom ligands to use for interaction graph generation. Ignored unless database = "custom"
#' When ligands is supplied, recepts must also be supplied and equidimensional.
#' @param recepts Character vector of custom receptors to use for interaction graph generation. Ignored unless database = "custom"
#' When recepts is supplied, ligands must also be supplied and equidimensional.
#'
#' @return Returns a dataframe containing gene symbols for ligand and receptor pairs. Column "pair" contains underscore-separated ligand-receptor pairs. Ligand and receptor gene names are stored in "source_genesymbol" and "target_genesymbol", respectively.
#' @import dplyr
#' @references Turei, et al. Molecular Systems Biology (2021); Raredon, et al. bioRxiv (2021)
#' @export
#'
#' @examples
LoadLR <- function(species = "human", database = "OmniPath", ligands = NULL, recepts = NULL) {
  if(database=="custom") {
    message("Using custom database")
    ligands <- ligands
    recepts <- recepts
    lit.put <- data.frame(pair.name = paste(ligands,recepts,sep="_"), ligands = ligands, recepts = recepts)
  }
  else {
    if(species %notin% c("human","mouse","rat")) {
      stop("Only human, mouse, and rat supported as species")
    }
    all <- readRDS(system.file(package = "scriabin", "lr_resources.rds"))[[species]]
    if(database %notin% names(all)) {
      stop("Database must be one of: OmniPath, CellChatDB, CellPhoneDB, Ramilowski2015, Baccin2019, LRdb, Kirouac2010, ICELLNET, iTALK, EMBRACE, HPMR, Guide2Pharma, connectomeDB2020, talklr, CellTalkDB\nFor rat or mouse, only OmniPath is supported")
    }
    message(paste("Using database",database))
    pairs <- as.data.frame(all[[database]][,c("source_genesymbol","target_genesymbol")] %>% mutate_all(as.character))
    lit.put <- pairs %>% dplyr::mutate(pair = paste(source_genesymbol,target_genesymbol, sep = "_"))
    lit.put <- as.data.frame(lit.put[,c("pair","source_genesymbol","target_genesymbol")])
  }
  return(lit.put)
}

#' Identify potential ligands for ligand activity prediction
#'
#' @param seu A Seurat object
#' @param assay Assay in Seurat object from which to pull expression values
#' @param slot Slot within assay from which to pull expression values
#' @param min.pct Minimum percentage of cells in which a gene must be detected in order to be considered an "expressed" gene. Default 0.025 (ie. 2.5%)
#' @param species character. Name of species from which to load ligand-receptor databases. One of: "human", "mouse", "rat". Default: "human"
#' @param database Name of ligand-receptor database to use. Default: "OmniPath"
#' When species is "human", one of: OmniPath, CellChatDB, CellPhoneDB, Ramilowski2015, Baccin2019, LRdb, Kirouac2010, ICELLNET, iTALK, EMBRACE, HPMR, Guide2Pharma, connectomeDB2020, talklr, CellTalkDB
#' When species is "mouse" or "rat", only "OmniPath" is supported.
#' To pass a custom ligand-receptor database to this function, set database = "custom"
#' @param ligands Character vector of custom ligands to use for interaction graph generation. Ignored unless database = "custom"
#' When ligands is supplied, recepts must also be supplied and equidimensional.
#' @param recepts Character vector of custom receptors to use for interaction graph generation. Ignored unless database = "custom"
#' When recepts is supplied, ligands must also be supplied and equidimensional.
#'
#' @return Returns a list of length 2: 1) a character vector of potential ligands, 2) a character vector of background expressed genes
#' @import dplyr
#' @export
#'
#' @examples
IDPotentialLigands <- function(seu, assay = "SCT", slot = "data", min.pct = 0.025,
                               species = "human", database = "OmniPath",
                               ligands = NULL, recepts = NULL) {
  if(!exists("ligand_target_matrix")) {
    stop("Error: Please load NicheNet database into environment via scriabin::load_nichenet_database()")
  }
  lr_network <- LoadLR(database = database, species = species, ligands = ligands, recepts = recepts)
  exprs <- GetAssayData(seu, assay = assay, slot = slot)
  expressed_genes <- rownames(exprs)[(Matrix::rowSums(exprs !=0)/ncol(exprs))>min.pct]
  background_expressed_genes <- expressed_genes %>% .[. %in% rownames(ligand_target_matrix)]
  ligands = lr_network %>% pull(source_genesymbol) %>% unique()
  receptors = lr_network %>% pull(target_genesymbol) %>% unique()
  expressed_ligands = intersect(ligands,expressed_genes)
  expressed_receptors = intersect(receptors,expressed_genes)
  potential_ligands = lr_network %>% dplyr::filter(source_genesymbol %in% expressed_ligands & target_genesymbol %in% expressed_receptors) %>% pull(source_genesymbol) %>% unique()
  return(list(potential_ligands,background_expressed_genes))
}




#' Convenient negation
#'
#' @return
#' @export
#'
#' @examples
`%notin%` <- Negate(`%in%`)

#' Resampling
#'
#' @return
#' @export
#'
#' @examples
resample <- function(x, ...) x[sample.int(length(x), ...)]


#' Load NicheNet Database
#'
#' @return
#' @export
#'
#' @examples
load_nichenet_database <- function() {
  f = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
  assign("ligand_target_matrix", f, envir = .GlobalEnv)
  lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
  assign("lr_network", lr_network, envir = .GlobalEnv)
  weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
  assign("weighted_networks", weighted_networks, envir = .GlobalEnv)
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
  assign("weighted_networks_lr", weighted_networks_lr, envir = .GlobalEnv)
}



#' Mapvalues (ripped from plyr)
#'
#' @param x the factor or vector to modify
#' @param from a vector of the items to replace
#' @param to a vector of replacement values
#' @param warn_missing print a message if any of the old values are not present in x
#'
#' @return
#' @export
#'
#' @examples
mapvalues <- function (x, from, to, warn_missing = TRUE)
{
  if (length(from) != length(to)) {
    stop("`from` and `to` vectors are not the same length.")
  }
  if (!is.atomic(x)) {
    stop("`x` must be an atomic vector.")
  }
  if (is.factor(x)) {
    levels(x) <- mapvalues(levels(x), from, to, warn_missing)
    return(x)
  }
  mapidx <- match(x, from)
  mapidxNA <- is.na(mapidx)
  from_found <- sort(unique(mapidx))
  if (warn_missing && length(from_found) != length(from)) {
    message("The following `from` values were not present in `x`: ",
            paste(from[!(1:length(from) %in% from_found)], collapse = ", "))
  }
  x[!mapidxNA] <- to[mapidx[!mapidxNA]]
  x
}


#' Detect Inflection Point of a Concave Curve (Elbox Method). Ripped from:
#' \url{https://github.com/ahasverus/elbow}
#'
#' This function implements the Elbow (or knee of a curve) method to detect
#'   the inflection point of a concave curve. More information on this method:
#'   \url{https://en.wikipedia.org/wiki/Elbow_method_(clustering)}.
#'
#' @param data a two-columns data frame (x and y respectively).
#' @param plot a boolean. If TRUE (default) curves are plotted.
#'
#' @return This function returns a 2-elements list with:
#'   - the value on x-Axis corresponding to the inflection point
#'   - a data frame with the original data and two additional columns used in
#'     the graphic.
#'
#' @details This function detects unique inflection point in a simple concave
#'   curve. The curve can be concave down/up with a positive/negative slope.
#'
#' @import stats graphics
#'
#' @export
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.com}
#'
#' @examples
#' ## Loading dataset ----
#' data(profiles)
#' head(profiles)
#'
#' ## Object returned ----
#' x <- elbow(profiles[ , c("x", "concave_down_pos_slo")], plot = FALSE)
#' class(x)
#' names(x)
#' x$"x_selected"
#'
#' ## Graphical usage ----
#' x <- elbow(profiles[ , c("x", "concave_down_pos_slo")])
#'
#' ## The four implemented profiles ----
#' curves <- colnames(profiles)[-1]
#' par(mfrow = c(2, 2))
#' for (i in curves) {
#'   elbow(profiles[ , c("x", i)])
#'   title(i)
#' }
#'

elbow <- function(data, plot = TRUE) {
  ## Argument checks ----

  if (missing(data)) {
    stop("Please provide a two-columns data frame.")
  }

  if (!is.list(data)) {
    stop("`data` must be a two-columns data frame.")
  }

  if (ncol(data) > 2) {
    warning("Only the first two columns will be considered.")
  }

  if (!is.numeric(data[ , 1]) || !is.numeric(data[ , 2])) {
    stop("Non-numeric data detected.")
  }

  if (sum(is.na(data[ , 1])) + sum(is.na(data[ , 2]))) {
    stop("Missing values detected.")
  }

  if (!is.logical(plot)) {
    stop("`plot` must be a boolean.")
  }


  ## Data transformation ----

  data <- data[ , 1:2]
  data <- data[order(data[ , 1]), ]


  ## Get constant increase/decrease in y ----

  constant <- data[c(1, nrow(data)), ]
  colnames(constant) <- c("x", "y")

  mod <- stats::lm(y ~ x, data = constant)

  data[ , "constant"] <- round(mod$coef[[1]] + mod$coef[[2]] * data[ , 1], 3)


  ## Detect inflection point ----

  pos <- round(nrow(data) / 2)

  if (data[pos, "constant"] < data[pos, 2]) { # Concave Down

    ymin <- min(data[ , 2])
    data[ , "benefits"] <- ymin + round(data[ , 2] - data[ , "constant"], 3)
    maxi <- data[which.max(data[ , "benefits"]), ]

  } else { # Concave Up

    ymax <- max(data[ , 2])
    data[ , "benefits"] <- ymax - round(data[ , "constant"] - data[ , 2], 3)
    maxi <- data[which.min(data[ , "benefits"]), ]
  }


  ## Store results ----

  xxx <- list()
  xxx[[1]] <- maxi[1, 1]
  xxx[[2]] <- data
  names(xxx) <- c(paste(colnames(data)[1], "selected", sep = "_"), "data")


  ## Plot ----

  if (plot) {

    xlims <- range(data[ , 1])
    ylims <- c(min(c(data[ , 2], data[ , 3], data[ , 4])), max(data[ , 2]))


    ## Graphical parameters ----

    graphics::par(
      mar      = c(2.5, 2.5, 1.5, 1.5),
      family   = "serif",
      cex.axis = 0.85,
      mgp      = c(2, .15, 0),
      tcl      = -0.25,
      fg       = "#666666",
      col      = "#666666",
      col.axis = "#666666"
    )


    ## Background plot ----

    graphics::plot(
      x    = data[ , 1],
      y    = data[ , 2],
      xlim = xlims,
      ylim = ylims,
      ann  = FALSE,
      axes = FALSE,
      type = "n"
    )

    graphics::grid()
    graphics::box()


    ## Add axes ----

    graphics::par(mgp = c(2, 0.00, 0))
    graphics::axis(1, lwd = 0)
    graphics::axis(1, maxi[ , 1], lwd = 0, font = 2, col.axis = "black")

    graphics::par(mgp = c(2, 0.25, 0))
    graphics::axis(2, lwd = 0, las = 1)
    at <- round(maxi[ , 2], 3)
    graphics::axis(2, at, lwd = 0, font = 2, col.axis = "black", las = 1)

    graphics::mtext(side = 1, cex = 1, line = 1.25, text = expression("x"))
    graphics::mtext(side = 2, cex = 1, line = 1.45, text = expression("f(x)"))


    ## Real gains/losses in y while x increases ----

    graphics::polygon(
      x      = c(data[ , 1], data[1, 1]),
      y      = c(data[ , "benefits"], data[1, "benefits"]),
      col    = "#aaaaaa66",
      border = "#aaaaaa"
    )


    ## Data serie ----

    graphics::points(
      x    = data[ , 1],
      y    = data[ , 2],
      type = "b",
      pch  = 19,
      col  = "black"
    )


    ## Inflection point informations ----

    graphics::lines(
      x   = rep(maxi[1, 1], 2),
      y   = c(par()$usr[3], maxi[1, 2]),
      col = "black",
      lwd = 0.5
    )

    graphics::lines(
      x   = c(par()$usr[1], maxi[1, 1]),
      y   = rep(maxi[1, 2], 2),
      col = "black",
      lwd = 0.5
    )

    graphics::points(
      x    = maxi[ , 1],
      y    = maxi[ , 2],
      type = "b",
      pch  = 19,
      cex  = 1.5,
      col  = "black"
    )

    graphics::points(
      x    = maxi[ , 1],
      y    = maxi[ , 4],
      type = "b",
      pch  = 19,
      cex  = 1,
      col  = "#666666"
    )
  }

  return(xxx)
}










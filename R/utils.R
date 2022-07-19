


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
#' @param lit_support Numeric. Only entries with curation_effort (number of unique database - citation pairs per interaction) greater than or equal to this value are retained. Default 7.
#'
#' @return Returns a dataframe containing gene symbols for ligand and receptor pairs. Column "pair" contains underscore-separated ligand-receptor pairs. Ligand and receptor gene names are stored in "source_genesymbol" and "target_genesymbol", respectively.
#' @import dplyr
#' @references Turei, et al. Molecular Systems Biology (2021); Raredon, et al. bioRxiv (2021)
#' @export
#'
#' @examples
LoadLR <- function(species = "human", database = "OmniPath", ligands = NULL, recepts = NULL, lit_support = 7) {
  if(database=="custom") {
    message("Using custom database")
    ligands <- ligands
    recepts <- recepts
    lit.put <- data.frame(pair.name = paste(ligands,recepts,sep="_"), source_genesymbol = ligands, target_genesymbol = recepts)
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
    pairs <- all[[database]] %>% dplyr::filter(curation_effort>=lit_support) %>%
      dplyr::select(source_genesymbol,target_genesymbol) %>% mutate_all(as.character) %>% as.data.frame()
    lit.put <- pairs %>% dplyr::mutate(pair = paste(source_genesymbol,target_genesymbol, sep = "_"))
    lit.put <- as.data.frame(lit.put[,c("pair","source_genesymbol","target_genesymbol")])
  }
  return(lit.put)
}

#' Load unformatted ligand-receptor database
#'
#' @param species character. Name of species from which to load ligand-receptor databases. One of: "human", "mouse", "rat". Default: "human"
#' @param database Name of ligand-receptor database to use. Default: "OmniPath"
#' When species is "human", one of: OmniPath, CellChatDB, CellPhoneDB, Ramilowski2015, Baccin2019, LRdb, Kirouac2010, ICELLNET, iTALK, EMBRACE, HPMR, Guide2Pharma, connectomeDB2020, talklr, CellTalkDB
#' When species is "mouse" or "rat", only "OmniPath" is supported.
#' To pass a custom ligand-receptor database to this function, set database = "custom"
#'
#' @return Returns a dataframe containing full information for the specified LR database. If species==NULL, will return all LR databases in the package. If database==NULL, will return all databases for the given species.
#' @import dplyr
#' @references Turei, et al. Molecular Systems Biology (2021); Raredon, et al. bioRxiv (2021)
#' @export
#'
#' @examples
LoadRawLR <- function(species = "human", database = "OmniPath") {
  if(is.null(species)) {
    all <- readRDS(system.file(package = "scriabin", "lr_resources.rds"))
    message("Returning all LR databases in package")
    return(all)
  } else if(species %notin% c("human","mouse","rat")) {
    stop("Only human, mouse, and rat supported as species")
  } else {
    all <- readRDS(system.file(package = "scriabin", "lr_resources.rds"))[[species]]
    if(is.null(database)) {
      message(paste0("Returning all LR databases for species: ",species))
      return(all)
    } else if(database %notin% names(all)) {
      stop("Database must be one of: OmniPath, CellChatDB, CellPhoneDB, Ramilowski2015, Baccin2019, LRdb, Kirouac2010, ICELLNET, iTALK, EMBRACE, HPMR, Guide2Pharma, connectomeDB2020, talklr, CellTalkDB\nFor rat or mouse, only OmniPath is supported")
    } else {
      return(all[[database]])
    }
  }
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
#' @import dplyr nichenetr
#' @export
#'
#' @examples
IDPotentialLigands <- function(seu, assay = "SCT", slot = "data", min.pct = 0.025,
                               species = "human", database = "OmniPath",
                               ligands = NULL, recepts = NULL) {
  if(!exists("ligand_target_matrix")) {
    stop("Error: Please load NicheNet database into environment via scriabin::load_nichenet_database()")
  }
  library(nichenetr)
  lr_network <- LoadLR(database = database, species = species, ligands = ligands, recepts = recepts)
  exprs <- GetAssayData(seu, assay = assay, slot = slot)
  if(species != "human") {
    rownames(exprs) <- nichenetr::convert_mouse_to_human_symbols(rownames(exprs))
    lr_network <- lr_network %>% dplyr::mutate(source_genesymbol = nichenetr::convert_mouse_to_human_symbols(source_genesymbol)) %>%
      dplyr::mutate(target_genesymbol = nichenetr::convert_mouse_to_human_symbols(target_genesymbol))
  }
  expressed_genes <- rownames(exprs)[(Matrix::rowSums(exprs !=0)/ncol(exprs))>min.pct]
  background_expressed_genes <- expressed_genes %>% .[. %in% rownames(ligand_target_matrix)]
  ligands = lr_network %>% pull(source_genesymbol) %>% unique()
  receptors = lr_network %>% pull(target_genesymbol) %>% unique()
  expressed_ligands = intersect(ligands,expressed_genes)
  expressed_receptors = intersect(receptors,expressed_genes)
  potential_ligands = lr_network %>% dplyr::filter(source_genesymbol %in% expressed_ligands) %>% pull(source_genesymbol) %>% unique()
  potential_ligands <- potential_ligands[!is.na(potential_ligands)]
  potential_ligands <- potential_ligands[potential_ligands %in% colnames(ligand_target_matrix)]
  background_expressed_genes <- background_expressed_genes[!is.na(background_expressed_genes)]

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
#'
#' @return This function returns a 2-elements list with:
#'   - the value on x-Axis corresponding to the inflection point
#'   - a data frame with the original data and two additional columns used in
#'     the graphic.
#'
#' @details This function detects unique inflection point in a simple concave
#'   curve. The curve can be concave down/up with a positive/negative slope.
#'
#' @export
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.com}
#'
#' @examples

elbow <- function(data) {
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

  return(xxx)
}










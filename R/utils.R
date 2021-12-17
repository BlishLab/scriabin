
#' Title
#'
#' @return
#' @export
#'
#' @examples
`%notin%` <- Negate(`%in%`)

#' Title
#'
#' @return
#' @export
#'
#' @examples
resample <- function(x, ...) x[sample.int(length(x), ...)]

#' Title
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

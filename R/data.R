#' @docType data
#' @name DataSet
#'
#' @aliases psp
#' @aliases ubi
#' @aliases ace
#' @aliases input_psp_example
#' @aliases input_ubi_example
#' @aliases input_ace_example
#' @aliases wholeProteome_psp_example
#' @aliases motif_example
#' @keywords datasets
NULL


#' @export psp_data
#' @export ubi_data
#' @export ace_data
#' @export inputExample_PSPdata
#' @export inputExample_UBIdata
#' @export inputExample_ACEdata
#' @export wholeProteomeExample_data
#' @export motifExample_data

#'
#'

psp_data <- function() {
  utils::data(list="PhosphositePlus", package="JUMPsem")
  get("psp", envir = .GlobalEnv)
}

ubi_data <- function() {
  utils::data(list="UbiquitinDatabase", package="JUMPsem")
  get("ubi", envir = .GlobalEnv)
}


ace_data <- function(){
  utils::data(list="AcetyDatabase", package="JUMPsem")
  get("ace", envir = .GlobalEnv)
}




inputExample_PSPdata <- function() {
  utils::data(list="input_psp_example", package="JUMPsem")
  get("input_psp_example", envir = .GlobalEnv)
}

inputExample_UBIdata <- function() {
  utils::data(list="input_ubi_example", package="JUMPsem")
  get("input_ubi_example", envir = .GlobalEnv)
}


inputExample_ACEdata <- function() {
  utils::data(list="input_ace_example", package="JUMPsem")
  get("input_ace_example", envir = .GlobalEnv)
}



wholeProteomeExample_data <- function() {
  utils::data(list="wholeProteome_psp_example", package="JUMPsem")
  get("wholeProteome_psp_example", envir = .GlobalEnv)
}


motifExample_data <- function() {
  utils::data(list="motif_example", package="JUMPsem")
  get("motif_example", envir = .GlobalEnv)
}









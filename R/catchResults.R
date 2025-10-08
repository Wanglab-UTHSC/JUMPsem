#' @import dplyr
#' @import stringr
#' @param x The results from modeling fitting.
#' @param rawData  Input raw data (quantitative phospho-proteomics / ubiqutin-proteomics data)
#' @title catchResults
#' @description Catch Activity and Affinity. Activity contains Raw-Activity, Mean-Center, Z-scaled.


#' @export catchActRaw
#' @export catchActMeanCenter
#' @export catchActZscore
#' @export centerApply
#' @rdname catchResults
#'
#'


### Raw activity ###
catchActRaw <- function(x){
  x.nonull <- x[sapply(x, typeof) == "list"]
  xdataframe <- do.call(cbind, lapply(x.nonull, data.frame))
  # transpose the matrix
  xdataframe_t <- data.frame(t(xdataframe))
  rownames(xdataframe_t) <- names(x.nonull)
  xdataframe_t$ENZYME <- rownames(xdataframe_t)
  final <- xdataframe_t[, c(ncol(xdataframe_t), 1:(ncol(xdataframe_t) - 1))]
  return(final)
}


### Mean-Center normalization ###
centerApply <- function(i) {
  apply(i, 2, function(y) y / mean(abs(y))) #raw_value / mean_of_absolute_values
}

catchActMeanCenter <- function(x){
  x.nonull <- x[sapply(x, typeof) == "list"]
  xdataframe <- do.call(cbind,lapply(x.nonull, data.frame))
  # Mean-center transform based on each enzyme
  xdataframe_m <- centerApply(xdataframe)
  # transpose the matrix
  xdataframe_t <- data.frame(t(xdataframe_m))
  rownames(xdataframe_t) <- names(x.nonull)
  xdataframe_t$ENZYME <- rownames(xdataframe_t)
  final <- xdataframe_t[, c(ncol(xdataframe_t), 1:(ncol(xdataframe_t) - 1))]
  return(final)
}


### Z-score normalization ###
catchActZscore <- function(x){
  x.nonull <- x[sapply(x, typeof) == "list"]
  xdataframe <- do.call(cbind,lapply(x.nonull, data.frame))
  # zscore transform based on each enzyme
  xdataframe_s <- scale(xdataframe)
  # transpose the matrix
  xdataframe_t <- data.frame(t(xdataframe_s))
  rownames(xdataframe_t) <- names(x.nonull)
  xdataframe_t$ENZYME <- rownames(xdataframe_t)
  final <- xdataframe_t[, c(ncol(xdataframe_t), 1:(ncol(xdataframe_t) - 1))]
  return(final)
}



#' @export catchAff
#' @rdname catchResult

### Catch Affinity ###
catchAff <- function(x, rawData){

  x.nonull <- x[sapply(x, typeof) == "list"]
  xdataframe <- do.call(rbind,lapply(x.nonull, data.frame))
  xdataframe$ENZYME <- rownames(xdataframe)
  xdataframe$ENZYME <- str_split(xdataframe$ENZYME, "[.]", simplify = TRUE)[, 1]
  # Join with substrate gene name
  xgenename <- separate(xdataframe, col = 1, into = c("Substrate_ACC", "sites"), sep = "_", remove = F)
  subgene <- unique(rawData[, 1:2])
  subgene <- data.frame(subgene)
  subgene$Substrate_ACC <- gsub("-", ".", subgene$Substrate_ACC)

  xfull <- left_join(xgenename, subgene, by = "Substrate_ACC")
  final <- xfull[, c(7, 2, 8, 3, 4:6)]
  # Rename the colnames
  colnames(final) <- c("ENZYME_GN", "SUB_ACC_ID", "SUBSTRATE_GN", "SUB_sites", "AFFINITY", "z", "Pvalue")
  return(final)
}











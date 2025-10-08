#' @description Get single HAT activity.
#' @export singleHatAct
#' @title singleHatAct
#' @import lavaan
#' @import EFAtools
#' @importFrom EFAtools KMO
#' @param HAT Single HAT name.
#' @param input Normalized and transformed acetyl-proteomics data.
#' @param adj Adjacency matrix.
#' @param cor.off Set up correlation cutoff value 0-1 to remove high collinear variables. Default is 0.95.
#' @param kmo.off Set up KMO cutoff value 0-1. Default is 0.




singleHatAct <- function(HAT,
                         input,
                         adj,
                         cor.off,
                         kmo.off)
{
  try({
    activity_list <- list()
    HAT_substrates <- rownames(adj[which(adj[, HAT] == 1), ])# extract "1" from special HAT column
    rownames(input) <- gsub("-", ".", rownames(input))

    intersect_substrates <- intersect(x = HAT_substrates, y = rownames(input))# intersection of HAT_substrates between adj_matrix and input

    if (length(intersect_substrates) > 1){
      # correlation between each HAT_substrate
      data <- as.matrix(input)
      HAT_table_raw <- t(data[which(rownames(data) %in% intersect_substrates),])

      # transfer matrix into numeric matrix (just make sure)
      HAT_table_raw_c <- as.numeric(HAT_table_raw)
      dim(HAT_table_raw_c) <- dim(HAT_table_raw)
      colnames(HAT_table_raw_c) <- colnames(HAT_table_raw)
      rownames(HAT_table_raw_c) <- rownames(HAT_table_raw)
      HAT_table_raw <- HAT_table_raw_c

      cor_matrix <- cor(HAT_table_raw) # Correlation matrix

      # Remove highly&lowly correlated variables
      # Modify correlation matrix
      cor_matrix[upper.tri(cor_matrix)] <- 0
      diag(cor_matrix) <- 0
      HAT_table_new <- HAT_table_raw[ ,!apply(cor_matrix, 2, function(x) any(x > cor.off | x < -cor.off ))]


      # filter with KMO() cutoff
      if(length(colnames(HAT_table_new)) > 0){
        HAT_table_md <- HAT_table_new[, which(as.matrix(KMO(cor(HAT_table_new))[[2]]) > kmo.off)]
      }else{
        HAT_table_md <- HAT_table_new
      }

      # SEM model
      if (length(colnames(HAT_table_md)) > 1){
        HAT_new_t <- as.matrix(t(HAT_table_md))
        HAT_new_t_rowname <- rownames(HAT_new_t)

        substrate <- noquote(paste(HAT_new_t_rowname, collapse = " + "))
        HAT_model <- paste(HAT, " =~", substrate)

        # fit & summary SEM model
        fit_HAT = suppressWarnings(sem(HAT_model, data = t(data)))
        summary_PE = quiet(summary(fit_HAT)[[5]])

        sub_num <- length(colnames(HAT_table_md))
        estimate <- as.vector(as.matrix(summary_PE$est)[1:sub_num])
        activity_list$HAT <- HAT_table_md %*% estimate
        print(HAT)

        return(activity_list)
      }
    }
  }, silent= T)}

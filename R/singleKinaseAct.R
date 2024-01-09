#' @description Get single kinase activity.
#' @export singleKinaseAct
#' @title singleKinaseAct
#' @import lavaan
#' @import EFAtools
#' @importFrom EFAtools KMO
#' @param kinase Single kinase name.
#' @param input Normalized and transformed phosphoproteomics  data.
#' @param adj Adjacency matrix.
#' @param cor.off Set up correlation cutoff value 0-1 to remove high collinear variables. Default is 0.95.
#' @param kmo.off Set up KMO cutoff value 0-1. Default is 0.


singleKinaseAct <- function(kinase,
                            input,
                            adj,
                            cor.off,
                            kmo.off)
{
  try({
    activity_list <- list()
    kinase_substrates <- rownames(adj[which(adj[, kinase] == 1), ])# extract "1" from special kinase column
    intersect_substrates <- intersect(x = kinase_substrates, y = rownames(input))# intersection of kinase_substrates between adj_matrix and input

    if (length(intersect_substrates) > 1){
      # correlation between each kinase_substrate
      data <- as.matrix(input)
      kinase_table_raw <- t(data[which(rownames(data) %in% intersect_substrates),])
      cor_matrix <- cor(kinase_table_raw) # Correlation matrix

      # Remove highly&lowly correlated variables
      # Modify correlation matrix
      cor_matrix[upper.tri(cor_matrix)] <- 0
      diag(cor_matrix) <- 0
      kinase_table_new <- kinase_table_raw[ ,!apply(cor_matrix, 2, function(x) any(x > cor.off | x < -cor.off ))]


      # filter with KMO() cutoff
      if(length(colnames(kinase_table_new)) > 0){
        kinase_table_md <- kinase_table_new[, which(as.matrix(KMO(cor(kinase_table_new))[[2]]) > kmo.off)]
      }else{
        kinase_table_md <- kinase_table_new
      }

      # SEM model
      if (length(colnames(kinase_table_md)) > 1){
        kinase_new_t <- as.matrix(t(kinase_table_md))
        kinase_new_t_rowname <- rownames(kinase_new_t)

        substrate <- noquote(paste(kinase_new_t_rowname, collapse = " + "))
        kinase_model <- paste(kinase, " =~", substrate)

        # fit & summary SEM model
        fit_kinase = suppressWarnings(sem(kinase_model, data = t(data)))
        summary_PE = quiet(summary(fit_kinase)[[5]])

        sub_num <- length(colnames(kinase_table_md))
        estimate <- as.vector(as.matrix(summary_PE$est)[1:sub_num])
        activity_list$kinase <- kinase_table_md %*% estimate
        print(kinase)

        return(activity_list)
      }
    }
  }, silent= T)}

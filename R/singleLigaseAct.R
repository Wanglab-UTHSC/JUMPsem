#' @description Get single ligase activity.
#' @export singleLigaseAct
#' @title singleLigaseAct
#' @import lavaan
#' @import EFAtools
#' @importFrom EFAtools KMO
#' @param ligase Single ligase name.
#' @param input Normalized and transformed ubiquiti-proteomics data.
#' @param adj Adjacency matrix.
#' @param cor.off Set up correlation cutoff value 0-1 to remove high collinear variables. Default is 0.95.
#' @param kmo.off Set up KMO cutoff value 0-1. Default is 0.

#ligase = "BRCA1"
#HERC2


singleLigaseAct <- function(ligase,
                            input,
                            adj,
                            cor.off,
                            kmo.off)
{
  try({
    activity_list <- list()
    ligase_substrates <- rownames(adj[which(adj[, ligase] == 1), ])# extract "1" from special ligase column
    intersect_substrates <- intersect(x = ligase_substrates, y = rownames(input))# intersection of ligase_substrates between adj_matrix and input

    if (length(intersect_substrates) > 1){
      # correlation between each ligase_substrate
      data <- as.matrix(input)
      ligase_table_raw <- t(data[which(rownames(data) %in% intersect_substrates),])
      cor_matrix <- cor(ligase_table_raw) # Correlation matrix

      # Remove highly&lowly correlated variables
      # Modify correlation matrix
      cor_matrix[upper.tri(cor_matrix)] <- 0
      diag(cor_matrix) <- 0
      ligase_table_new <- ligase_table_raw[ ,!apply(cor_matrix, 2, function(x) any(x > cor.off | x < -cor.off ))]


      # filter with KMO() cutoff
      if(length(colnames(ligase_table_new)) > 0){
        ligase_table_md <- ligase_table_new[, which(as.matrix(KMO(cor(ligase_table_new))[[2]]) > kmo.off)]
      }else{
        ligase_table_md <- ligase_table_new
      }

      # SEM model
      if (length(colnames(ligase_table_md)) > 1){
        ligase_new_t <- as.matrix(t(ligase_table_md))
        ligase_new_t_rowname <- rownames(ligase_new_t)

        substrate <- noquote(paste(ligase_new_t_rowname, collapse = " + "))
        ligase_model <- paste(ligase, " =~", substrate)

        # fit & summary SEM model
        fit_ligase = suppressWarnings(sem(ligase_model, data = t(data)))
        summary_PE = quiet(summary(fit_ligase)[[5]])

        sub_num <- length(colnames(ligase_table_md))
        estimate <- as.vector(as.matrix(summary_PE$est)[1:sub_num])
        activity_list$ligase <- ligase_table_md %*% estimate
        print(ligase)

        return(activity_list)
      }
    }
  }, silent= T)}

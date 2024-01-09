#' @description Get single enzyme affinity.
#' @export
#' @import lavaan
#' @import EFAtools
#' @importFrom EFAtools KMO
#' @param enzyme Single enzyme name.
#' @param input Normalized and transformed phospho-proteomics / ubiqutin-proteomics data.
#' @param adj Adjacency matrix.
#' @param cor.off Set up correlation cutoff value 0-1 to remove high collinear variables. Default is 0.95.
#' @param kmo.off Set up KMO cutoff value 0-1. Default is 0.
#' @title singleEnzymeAff


singleEnzymeAff <- function(enzyme,
                            input,
                            adj,
                            cor.off,
                            kmo.off)
{
  try({
    affinitiy_list <- list()
    enzyme_substrates <- rownames(adj[which(adj[, enzyme] == 1), ])# extract "1" from special enzyme column
    intersect_substrates <- intersect(x = enzyme_substrates, y = rownames(input))# intersection of enzyme_substrates between adj_matrix and input

    if (length(intersect_substrates) > 1){
      # correlation between each enzyme_substrate
      data <- as.matrix(input)
      enzyme_table_raw <- t(data[which(rownames(data) %in% intersect_substrates),])

      # transfer matrix into numeric matrix (just make sure)
      enzyme_table_raw_c <- as.numeric(enzyme_table_raw)
      dim(enzyme_table_raw_c) <- dim(enzyme_table_raw)
      colnames(enzyme_table_raw_c) <- colnames(enzyme_table_raw)
      rownames(enzyme_table_raw_c) <- rownames(enzyme_table_raw)
      enzyme_table_raw <- enzyme_table_raw_c


      cor_matrix <- cor(enzyme_table_raw) # Correlation matrix

      # Remove highly&lowly correlated variables
      # Modify correlation matrix
      cor_matrix[upper.tri(cor_matrix)] <- 0
      diag(cor_matrix) <- 0
      enzyme_table_new <- enzyme_table_raw[ ,!apply(cor_matrix, 2, function(x) any(x > cor.off | x < -cor.off ))]

      # filter with KMO() cutoff
      if(length(colnames(enzyme_table_new)) > 0){
        enzyme_table_md <- enzyme_table_new[, which(as.matrix(KMO(cor(enzyme_table_new))[[2]]) > kmo.off)]
      }else{
        enzyme_table_md <- enzyme_table_new
      }

      # SEM model
      if (length(colnames(enzyme_table_md)) > 1){
        enzyme_new_t <- as.matrix(t(enzyme_table_md))
        enzyme_new_t_rowname <- rownames(enzyme_new_t)

        substrate <- noquote(paste(enzyme_new_t_rowname, collapse = " + "))
        enzyme_model <- paste(enzyme, " =~", substrate)

        # fit & summary SEM model
        fit_enzyme = suppressWarnings(sem(enzyme_model, data = t(data)))
        summary_PE = quiet(summary(fit_enzyme)[[5]])
        main_table <- summary_PE[which(summary_PE$op == "=~"),]
        main_table <- main_table[, c(3, 5, 7, 8)]
        colnames(main_table) <- c("substrate-sites", "affinity", "z", "Pvalue")
        affinitiy_list$enzyme <- main_table
        print(enzyme)

        return(affinitiy_list)
      }
    }
  }, silent= T)}


###
# adj_matrix <- pspAdjacency(sp = "human")
# enzyme = "EP300"
# load("C:\\Project_Dehui_2022\\KSEM\\02_KSEM_code\\KSEM_v1\\data\\input_example.RData")
# input_trans <- inputTransform(input_raw = input,
#                               input_log2_trans = F)
#
# enzylist <- as.list(colnames(adj_matrix))
# names(enzylist) <- colnames(adj_matrix)
#
# test <- lapply(enzylist,
#                singleKinaseAff,
#                input= input_trans,
#                adj = adj_matrix,
#                cor.off = 0.85,
#                kmo.off = 0)






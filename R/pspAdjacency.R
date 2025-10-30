#' @title pspAdjacency
#' @export pspAdjacency
#' @import dplyr
#' @import tidyr
#' @description  This script is used to build kinase-substrate adjacency matrix with only PSP database.
#' @param sp Substrate species; ("mouse", "rat", "human").
#' @param ep Enzyme species; ("mouse", "rat", "human").
#' @param databasePSP Default database or customized database input.
#' @param mdsite Logical. Mapping with precise phosphorylated modification sites or not. Default is TRUE.


#'

pspAdjacency <- function(sp, ep, databasePSP, mdsite){
  # Decide using default database or user customized database
  if(is.null(databasePSP)){ # default
    psp <- psp_data()
    psp_rel <- psp
    rm(psp)
  }else{
    psp_rel <- databasePSP # customized
  }

  # Select PSP matrix
  if(is.null(ep)){
    ep <- sp
  }

  psp_set <- psp_rel[which(psp_rel$KIN_ORGANISM %in% ep & psp_rel$SUB_ORGANISM %in% sp), ]

  # Build psite name and replace "-" with "."
  psp_set$psite <- paste(psp_set$SUB_ACC_ID, "_", psp_set$SUB_MOD_RSD, sep = "")
  psp_set$psite <- as.matrix(gsub("-", ".", as.matrix(psp_set$psite)))
  psp_set$SUB_ACC_ID <- as.matrix(gsub("-", ".", as.matrix(psp_set$SUB_ACC_ID)))

  # Change all kinase names into toupper
  psp_set$GENE <- toupper(psp_set$GENE)

  if (mdsite == TRUE){
    # Duplicate and Add value=1
    psp_set <- psp_set[, c(1, 17)]

  }else{
    # Duplicate and Add value=1
    psp_set <- psp_set[, c(1, 7)]
  }

  psp_set <- psp_set[!duplicated(psp_set), ]
  psp_set$value <- 1
  rel_adj_matrix <- psp_set %>%
    pivot_wider(names_from ='GENE',
                values_from ='value',
                values_fill = list(value = 0))
  rel_adj_matrix <- as.matrix(rel_adj_matrix)
  rownames(rel_adj_matrix) <- rel_adj_matrix[, 1]
  rel_adj_matrix <- rel_adj_matrix[, -1]
  #row_n <- dim(rel_adj_matrix)[1]
  #col_n <- dim(rel_adj_matrix)[2]
  rel_adj <- apply(rel_adj_matrix, 2, as.numeric)
  rownames(rel_adj) <- rownames(rel_adj_matrix)
  return(rel_adj)
}

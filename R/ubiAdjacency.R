#' @title ubiAdjacency
#' @export ubiAdjacency
#' @import dplyr
#' @import tidyr
#' @description  This script is used to build ligase-substrate adjacency matrix with only PSP database.
#' @param sp Substrate species; ("mouse", "rat", "human").
#' @param databaseUBI Default database or customized database input.
#'


ubiAdjacency <- function(sp, databaseUBI){
  if(is.null(databaseUBI)){
    ubi <- ubi_data()
    ubi_rel <- ubi
    rm(ubi)
  }else{
    ubi_rel <- databaseUBI
  }

  # Select ubiquitination matrix
  ubi_set <- ubi_rel[which(ubi_rel$organism %in% sp), ]
  # Rebuild E3 name and replace "/" with "."
  ubi_set$Gene.Symbol..E3. <- as.matrix(gsub("/", ".", as.matrix(ubi_set$Gene.Symbol..E3.)))
  # Change all kinase names into toupper
  ubi_set$Gene.Symbol..E3. <- toupper(ubi_set$Gene.Symbol..E3)

  # split ligase and separate into multiple rows with "#"
  ubi_set <- ubi_set[, c(5, 4)]
  ubi_set <- ubi_set %>%
    separate_rows(Gene.Symbol..E3., sep = "#")

  # Duplicate and Add value=1
  ubi_set <- ubi_set[!duplicated(ubi_set), ]
  ubi_set$value <- 1

  rel_adj_matrix <- ubi_set %>%
    pivot_wider(names_from = "Gene.Symbol..E3.",
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



#' @title aceAdjacency
#' @export aceAdjacency
#' @import dplyr
#' @import tidyr
#' @param databaseACE Default database or customized database input.
#' @description  This script is used to build HUMAN HAT-substrate adjacency matrix.




aceAdjacency <- function(databaseACE){
  if(is.null(databaseACE)){
    ace <- ace_data()
    ace <- as.data.frame(ace)
    ace_data <- ace
    rm(ace)
  }else{
    ace_data <- databaseACE
  }


  # extract enzyme gene name and substract uniprot id
  ace_set <- as.data.frame(ace_data[, c(1, 5)])

  # Duplicate and Add value=1
  ace_set <- ace_set[!duplicated(ace_set), ]
  ace_set$value <- 1

  rel_adj_matrix <- ace_set %>%
    pivot_wider(names_from = "HAT",
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



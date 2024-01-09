#' @title pspMotifAdjacency
#' @export pspMotifAdjacency
#' @import dplyr
#' @import tidyr
#' @description This script is used to build enzyme-substrate adjacency matrix with PSP database + Motif.
#' @param sp Substrate species; ("mouse", "rat", "human").
#' @param ep Enzyme species; ("mouse", "rat", "human").
#' @param databasePSP Default database or customized database input.

#' @param motif.ref Added enzyme-subatrate relationships from motif discovery. Same as KSEM -motif parameter.


pspMotifAdjacency <- function(sp, ep, databasePSP, motif.ref){
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
  psp_set$psite <- as.matrix(gsub("-",".",as.matrix(psp_set$psite)))
  # Change all kinase names into toupper
  psp_set$GENE <- toupper(psp_set$GENE)
  # Duplicate and Add value=1
  psp_set <- psp_set[, c(1, 17)]
  psp_set <- psp_set[!duplicated(psp_set), ]
  psp_set$value <- 1

  # add motif
  ### MOTIF ###
  # Duplicate and Add value=1
  motif.ref <- motif.ref[!duplicated(motif.ref), ]
  motif.ref$value = 1
  # combine PSP and Motif
  rel_adj_pspmotif <- rbind(psp_set, motif.ref)
  # Remove dulicated pairs
  rel_adj_pspmotif <- rel_adj_pspmotif[!duplicated(rel_adj_pspmotif), ]


  rel_adj_matrix <- rel_adj_pspmotif %>%
    pivot_wider(names_from ='GENE',
                values_from ='value',
                values_fill = list(value = 0))

  rel_adj_matrix <- as.matrix(rel_adj_matrix)
  rownames(rel_adj_matrix) <- rel_adj_matrix[, 1]
  rel_adj_matrix <- rel_adj_matrix[, -1]

  rel_adj <- apply(rel_adj_matrix, 2, as.numeric)
  rownames(rel_adj) <- rownames(rel_adj_matrix)
  return(rel_adj)
}


##
# ksp = "human"
# ssp = "human"
# motif.ref <- read.csv("C:\\Project_Dehui_2022\\KSEM\\04_Example\\motif_example.csv")










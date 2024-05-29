#' @title inputTransform
#' @export inputTransformPSP
#' @export inputTransformUBI
#' @description This script is used to log2 transform AND/OR whole proteome normalize input data.
#' @param input_raw Input raw data (quantitative phospho-proteomics/ubiquitin-proteomics data)
#' @param input_log2_trans FALSE or TRUE. Need program to do log2 transforming of the input file or not. Default is FALSE.
#' @param whole_log2_trans Need program to do log2 transforming of the whole proteome or not. Ignore if -whole.proteome is missing. Default is FALSE.
#' @param whole_proteome Set up whole proteome used to normalize phospho-proteomics/ubiquitin-proteomics data by whole proteome. Default is NULL.



####### phospho-proteomics #########
inputTransformPSP <- function(input_raw,
                              input_log2_trans = FALSE,
                              whole_log2_trans = FALSE,
                              whole_proteome = NULL){
  ## extract "Substrate_ACC", "psites" and samples
  input_1 <- input_raw[, c(1, 3, 5:ncol(input_raw))]
  ## separate rows with "," in "psites"
  input_2 <- input_1 %>%
    separate_rows(sites, sep = ",")
  ## delete rows with empty psites
  input_3 <- input_2[!(input_2$sites == "" | is.na(input_2$sites)), ]
  ## transform id into "SUB_ACC_psite" format
  input_3$psite <- paste(input_3$Substrate_ACC, "_", input_3$sites, sep = "")
  ## re-organize the colume order
  input_4 <- input_3[,c(ncol(input_3), 3:(ncol(input_3) - 1))]
  ## unique table
  index <- duplicated(input_4[, 1])
  data_uniq <- input_4[!index, ]
  data <- as.matrix(data_uniq[, 2:ncol(data_uniq)])
  rownames(data) <- data_uniq$psite

  if (input_log2_trans == TRUE){
    input_log2 <- log2(data)
  }else{
    input_log2 <- data
  }

  input_log2 <- data.frame(input_log2, check.names = F)

  if (!is.null(whole_proteome)){
    if (whole_log2_trans == TRUE) {
      whole_proteome_log2 <- log2(whole_proteome[, 3:ncol(whole_proteome)])
      whole_proteome_log2$Substrate_ACC <- whole_proteome[, 1]
      whole_proteome_log2 <- whole_proteome_log2[, c(ncol(whole_proteome_log2), 1:(ncol(whole_proteome_log2)-1))]

      # add Substrate_ACC
      input_log2$Substrate_ACC <- str_split(rownames(input_log2), "_", simplify = TRUE)[, 1]
      # order the colnames with whole_proteome_log2 colnames and add psite_id
      input_log2 <- input_log2[, colnames(whole_proteome_log2)]
      input_log2$id <- rownames(input_log2)
      # left join and remove rows with "NA" whole_proteome
      input_whole_join_na <- left_join(input_log2, whole_proteome_log2, by = "Substrate_ACC")
      input_whole_join <- input_whole_join_na[complete.cases(input_whole_join_na), ]
      # sample number and sample names
      smp_num <- ncol(input_log2) - 2
      smp_name <- colnames(whole_proteome_log2)[2:(smp_num + 1)]
      # log2 transform with whole proteome data
      input_whole_trans <- cbind(input_whole_join$Substrate_ACC, input_whole_join$id, input_whole_join[, 2:(smp_num + 1)] - input_whole_join[, (smp_num + 3):ncol(input_whole_join)])
      # row.name and col.name
      rownames(input_whole_trans) <- input_whole_trans[,2]
      input_whole_trans <- input_whole_trans[, 3:(smp_num + 2)]
      colnames(input_whole_trans) <- smp_name
      return(input_whole_trans)

    }else{
      whole_proteome_log2 <- whole_proteome[, -2]

      # add Substrate_ACC
      input_log2$Substrate_ACC <- str_split(rownames(input_log2), "_", simplify = TRUE)[, 1]
      # order the colnames with whole_proteome_log2 colnames and add psite_id
      input_log2 <- input_log2[, colnames(whole_proteome_log2)]
      input_log2$id <- rownames(input_log2)
      # left join and remove rows with "NA" whole_proteome
      input_whole_join_na <- left_join(input_log2, whole_proteome_log2, by = "Substrate_ACC")
      input_whole_join <- input_whole_join_na[complete.cases(input_whole_join_na), ]
      # sample number and sample names
      smp_num <- ncol(input_log2)-2
      smp_name <- colnames(whole_proteome_log2)[2:(smp_num + 1)]
      # log2 transform with whole proteome data
      input_whole_trans <- cbind(input_whole_join$Substrate_ACC, input_whole_join$id, input_whole_join[, 2:(smp_num + 1)] - input_whole_join[, (smp_num + 3):ncol(input_whole_join)])
      # row.name and col.name
      rownames(input_whole_trans) <- input_whole_trans[,2]
      input_whole_trans <- input_whole_trans[, 3:(smp_num + 2)]
      colnames(input_whole_trans) <- smp_name
      return(input_whole_trans)
    }
  }else{
    print("Running without whole proteome normalized!")
    return(input_log2)
  }

}



######### ubiquitin-proteomics ###########
inputTransformUBI <- function(input_raw,
                              input_log2_trans = FALSE,
                              whole_log2_trans = FALSE,
                              whole_proteome = NULL){
  ## extract "Substrate_ACC", and samples
  input_1 <- input_raw[, c(1, 5:ncol(input_raw))]

  data <- as.matrix(input_1[, 2:ncol(input_1)])
  rownames(data) <- input_1[, 1]

  if (input_log2_trans == TRUE){
    input_log2 <- log2(data)
  }else{
    input_log2 <- data
  }

  input_log2 <- data.frame(input_log2)

  if (!is.null(whole_proteome)){
    if (whole_log2_trans == TRUE) {
      whole_proteome_log2 <- log2(whole_proteome[, 3:ncol(whole_proteome)])
      whole_proteome_log2$Substrate_ACC <- whole_proteome[, 1]
      whole_proteome_log2 <- whole_proteome_log2[, c(ncol(whole_proteome_log2), 1:(ncol(whole_proteome_log2)-1))]

      # add Substrate_ACC
      input_log2$Substrate_ACC <- str_split(rownames(input_log2), "_", simplify = TRUE)[, 1]
      # order the colnames with whole_proteome_log2 colnames and add psite_id
      input_log2 <- input_log2[, colnames(whole_proteome_log2)]
      input_log2$id <- rownames(input_log2)
      # left join and remove rows with "NA" whole_proteome
      input_whole_join_na <- left_join(input_log2, whole_proteome_log2, by = "Substrate_ACC")
      input_whole_join <- input_whole_join_na[complete.cases(input_whole_join_na), ]
      # sample number and sample names
      smp_num <- ncol(input_log2)-2
      smp_name <- colnames(whole_proteome_log2)[2:(smp_num + 1)]
      # log2 transform with whole proteome data
      input_whole_trans <- cbind(input_whole_join$Substrate_ACC, input_whole_join$id, input_whole_join[, 2:(smp_num + 1)] - input_whole_join[, (smp_num + 3):ncol(input_whole_join)])
      # row.name and col.name
      rownames(input_whole_trans) <- input_whole_trans[,2]
      input_whole_trans <- input_whole_trans[, 3:(smp_num + 2)]
      colnames(input_whole_trans) <- smp_name
      return(input_whole_trans)

    }else{
      whole_proteome_log2 <- whole_proteome[, -2]

      # add Substrate_ACC
      input_log2$Substrate_ACC <- str_split(rownames(input_log2), "_", simplify = TRUE)[, 1]
      # order the colnames with whole_proteome_log2 colnames and add psite_id
      input_log2 <- input_log2[, colnames(whole_proteome_log2)]
      input_log2$id <- rownames(input_log2)
      # left join and remove rows with "NA" whole_proteome
      input_whole_join_na <- left_join(input_log2, whole_proteome_log2, by = "Substrate_ACC")
      input_whole_join <- input_whole_join_na[complete.cases(input_whole_join_na), ]
      # sample number and sample names
      smp_num <- ncol(input_log2)-2
      smp_name <- colnames(whole_proteome_log2)[2:(smp_num + 1)]
      # log2 transform with whole proteome data
      input_whole_trans <- cbind(input_whole_join$Substrate_ACC, input_whole_join$id, input_whole_join[, 2:(smp_num + 1)] - input_whole_join[, (smp_num + 3):ncol(input_whole_join)])
      # row.name and col.name
      rownames(input_whole_trans) <- input_whole_trans[,2]
      input_whole_trans <- input_whole_trans[, 3:(smp_num + 2)]
      colnames(input_whole_trans) <- smp_name
      return(input_whole_trans)
    }
  }else{
    print("Running without whole proteome normalized!")
    return(input_log2)
  }

}




######### acetyl-proteomics ###########
inputTransformACE <- function(input_raw,
                              input_log2_trans = FALSE,
                              whole_log2_trans = FALSE,
                              whole_proteome = NULL){
  ## extract "Substrate_ACC", and samples
  input_1 <- input_raw[, c(1, 5:ncol(input_raw))]

  data <- as.matrix(input_1[, 2:ncol(input_1)])
  rownames(data) <- input_1[, 1]

  if (input_log2_trans == TRUE){
    input_log2 <- log2(data)
  }else{
    input_log2 <- data
  }

  input_log2 <- data.frame(input_log2)

  if (!is.null(whole_proteome)){
    if (whole_log2_trans == TRUE) {
      whole_proteome_log2 <- log2(whole_proteome[, 3:ncol(whole_proteome)])
      whole_proteome_log2$Substrate_ACC <- whole_proteome[, 1]
      whole_proteome_log2 <- whole_proteome_log2[, c(ncol(whole_proteome_log2), 1:(ncol(whole_proteome_log2)-1))]

      # add Substrate_ACC
      input_log2$Substrate_ACC <- str_split(rownames(input_log2), "_", simplify = TRUE)[, 1]
      # order the colnames with whole_proteome_log2 colnames and add psite_id
      input_log2 <- input_log2[, colnames(whole_proteome_log2)]
      input_log2$id <- rownames(input_log2)
      # left join and remove rows with "NA" whole_proteome
      input_whole_join_na <- left_join(input_log2, whole_proteome_log2, by = "Substrate_ACC")
      input_whole_join <- input_whole_join_na[complete.cases(input_whole_join_na), ]
      # sample number and sample names
      smp_num <- ncol(input_log2)-2
      smp_name <- colnames(whole_proteome_log2)[2:(smp_num + 1)]
      # log2 transform with whole proteome data
      input_whole_trans <- cbind(input_whole_join$Substrate_ACC, input_whole_join$id, input_whole_join[, 2:(smp_num + 1)] - input_whole_join[, (smp_num + 3):ncol(input_whole_join)])
      # row.name and col.name
      rownames(input_whole_trans) <- input_whole_trans[,2]
      input_whole_trans <- input_whole_trans[, 3:(smp_num + 2)]
      colnames(input_whole_trans) <- smp_name
      return(input_whole_trans)

    }else{
      whole_proteome_log2 <- whole_proteome[, -2]

      # add Substrate_ACC
      input_log2$Substrate_ACC <- str_split(rownames(input_log2), "_", simplify = TRUE)[, 1]
      # order the colnames with whole_proteome_log2 colnames and add psite_id
      input_log2 <- input_log2[, colnames(whole_proteome_log2)]
      input_log2$id <- rownames(input_log2)
      # left join and remove rows with "NA" whole_proteome
      input_whole_join_na <- left_join(input_log2, whole_proteome_log2, by = "Substrate_ACC")
      input_whole_join <- input_whole_join_na[complete.cases(input_whole_join_na), ]
      # sample number and sample names
      smp_num <- ncol(input_log2)-2
      smp_name <- colnames(whole_proteome_log2)[2:(smp_num + 1)]
      # log2 transform with whole proteome data
      input_whole_trans <- cbind(input_whole_join$Substrate_ACC, input_whole_join$id, input_whole_join[, 2:(smp_num + 1)] - input_whole_join[, (smp_num + 3):ncol(input_whole_join)])
      # row.name and col.name
      rownames(input_whole_trans) <- input_whole_trans[,2]
      input_whole_trans <- input_whole_trans[, 3:(smp_num + 2)]
      colnames(input_whole_trans) <- smp_name
      return(input_whole_trans)
    }
  }else{
    print("Running without whole proteome normalized!")
    return(input_log2)
  }

}

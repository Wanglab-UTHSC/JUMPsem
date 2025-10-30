
#' @export JUMPsem
#' @title JUMPsem
#' @description Infer enzyme activity from phosphoralation, ubiquitination, acetylation
#' @details You can use JUMPsem package to infer enzyme activity.
#'
#' @param datatype Factor. Data type: ubiquitination or phosphoralation or acetylation; Corresponding values are "ubi", "psp" and "ace".
#' @param input Dataframe. Your quantitative ubiquiti-proteomics, phospho-proteomic and acetyl-proteomics data.
#' @param organism Vector. Substrate species. Corresponding values are “human”, “mouse” and “rat”.
#' @param enzyme.organism Vector. Enzyme species; Corresponding values are “human”, “mouse” and “rat”.
#' @param database Dataframe. Customized database input. Default is NULL using internal database.
#' @param kmo.off Numeric. Set up KMO cutoff value 0-1. Default is 0.
#' @param mdsite Logical. Mapping with precise phosphorylated modification sites or not. Default is TRUE.
#' @param enzyList Vector. Program only calculate the enzyme in the enzyList. Default is to output ALL enzyme activities and affinities.
#' @param input.log2.norm Logical. FALSE or TRUE. Need program to do log2 transforming of the input file or not. Default is FALSE.
#' @param relative.norm.p Logical. FALSE or TRUE. Need program to do relative normalization of the PTM input file or not. Default is TRUE.
#' @param whole.log2.trans Logical. FALSE or TRUE. Need program to do log2 transforming of the whole proteome or not. (Ignore if -whole.proteome is missing). Default is FALSE.
#' @param whole.proteome Dataframe. Set up whole proteome used to normalize phosphor-proteome or ubiquitin-proteome whole proteome. Default is NULL.
#' @param relative.norm.w Logical. FALSE or TRUE. Need program to do relative normalization of the whole proteomic input file or not. Default is TRUE.
#' @param motif Matrix. Added kinase-substrate relationships from motif discovery.
#' @param output.folder Character. Character vector of location to save files if desired. Default is current directory.
#' @return List
#'
#' @import lavaan
#' @import dplyr
#' @import tidyr
#' @import devtools
#' @import psych
#' @import MASS
#' @import Matrix
#' @import tidyverse
#' @import data.table


#' @importFrom stats complete.cases
#' @importFrom stats cor
#' @importFrom utils write.table
#' @examples
#'
#' result <- JUMPsem(input = input_ubi_example,
#'                   datatype = "ubi",
#'                   organism = "human",
#'                   input.log2.norm = T)
#'
#'
#' result <- JUMPsem(input = input_psp_example,
#'                   datatype = "psp",
#'                   organism = "mouse",
#'                   enzyme.organism = c("human", "mouse", "rat"),
#'                   input.log2.norm = TRUE,
#'                   whole.log2.trans = TRUE,
#'                   motif = motif_example,
#'                   whole.proteome = wholeProteome_example)
#'
#'
#' result <- JUMPsem(input = input_ace_example,
#'                   datatype = "ace",
#'                   organism = "human",
#'                   input.log2.norm = FALSE)
#'
#'
#'
#'



JUMPsem = function(
    input = NULL,
    datatype = NULL,
    organism = NULL,
    enzyme.organism = NULL,
    database = NULL,
    kmo.off = 0,
    mdsite = TRUE,
    enzyList = NULL,
    input.log2.norm = FALSE,
    relative.norm.p = TRUE,
    whole.log2.trans = FALSE,
    whole.proteome = NULL,
    relative.norm.w = TRUE,
    motif = NULL,
    output.folder = getwd()){


  # Check input file
  if (is.null(input)){
    stop("Please set your input file!")
  }else{
    message(date()) # print current time
  }

  # Check data type
  if (is.null(datatype)){
    stop("Please set your datatype!")
  }else{
    if (datatype == "ubi"){
      message("Your analysis is based on ubiquitination database!")
    }else if(datatype == "psp"){
      message("Your analysis is based on phosphoralation database!")
    }else if(datatype == "ace"){
      message("Your analysis is based on acetylation database!")
    }
  }

  # Check organism
  if (is.null(organism)){
    stop("Please set your data organism!")
  }else if(datatype == "ace" && organism != "human"){
    stop("Sorry but acetylation datatype is only available for human organism! \n Please change correct your datatype or organism")
  }else {
    message("Start running...")
  }

  # Check enzyme organism
  if (is.null(enzyme.organism)){
    enzyme.organism <- organism
  }else if(datatype == "ace" & all(enzyme.organism != "human")){
    stop("Sorry but currently enzyme organism can only be customerised for phosphorlation data! \n Please change correct your datatype or enzyme organism \n NOTE: You may need to leave 'enzyme.organism' as default or set as 'human' based on your 'ace' datatype.")
  }else if(datatype == "ubi" && enzyme.organism != organism){
    stop("Sorry but currently enzyme organism can only be customerised for phosphorlation data! \n Please change correct your datatype or enzyme organism \n NOTE: You may need to leave 'enzyme.organism' as default or set 'enzyme.organism' same as 'organism' based on your 'ubi' datatype.")
  }


  # enzyme organism info
  message("\n", "Enzyme organism: ")
  message(cat(enzyme.organism))

  # enzyme list
  if (is.null(enzyList)) {
    message("\n", "JUMPsem will calculate: ")
    message(cat("All enzymes."))
  }else{
    message("\n", "JUMPsem will calculate: ")
    message(cat(enzyList))
  }

  try({
    # transform the raw input into program-fitted input
    if (datatype == "psp"){
      input_trans <- inputTransformPSP(input_raw = input,
                                       input_log2_trans = input.log2.norm,
                                       relative.norm.p = relative.norm.p,
                                       whole_log2_trans = whole.log2.trans,
                                       whole_proteome = whole.proteome,
                                       relative.norm.w = relative.norm.w)
    }else if(datatype == "ubi"){
      input_trans <- inputTransformUBI(input_raw = input,
                                       input_log2_trans = input.log2.norm,
                                       relative.norm.p = relative.norm.p,
                                       whole_log2_trans = whole.log2.trans,
                                       whole_proteome = whole.proteome,
                                       relative.norm.w = relative.norm.w)
    }else if(datatype == "ace"){
      input_trans <- inputTransformACE(input_raw = input,
                                       input_log2_trans = input.log2.norm,
                                       relative.norm.p = relative.norm.p,
                                       whole_log2_trans = whole.log2.trans,
                                       whole_proteome = whole.proteome,
                                       relative.norm.w = relative.norm.w)
    }


    # print kmo.off info
    message("Your kmo.off: ", kmo.off)


    # check and build adjacency matrix
    if (is.null(motif)) {
      message("\n", "Running without motif added!")
      if (datatype == "psp"){
        adj_matrix <- pspAdjacency(sp = organism, ep = enzyme.organism, databasePSP = database, mdsite = mdsite)

      }else if(datatype == "ubi"){
        adj_matrix <- ubiAdjacency(sp = organism, databaseUBI = database)

      }else if(datatype == "ace"){
        adj_matrix <- aceAdjacency(databaseACE = database)
      }

    }else{
      message("\n", "Running with motif added!")
      adj_matrix <- pspMotifAdjacency(sp = organism, ep = enzyme.organism, databasePSP = database, motif.ref = motif, mdsite = mdsite)
    }
    message("\n", "Adjacency matrix completes!")
    # check and get kinase.list
    if (is.null(enzyList)){
      enzyList <- as.list(colnames(adj_matrix))
      names(enzyList) <- colnames(adj_matrix)
    }else{
      enzyme_list <- as.list(enzyList)
      names(enzyme_list) <- enzyList
      enzyList <- enzyme_list
    }
    # get activity and affinity
    message("\n", "Calculating activity begins!", "\n", "Valid Enzyme with activities: ")
    if (datatype == "psp"){
      result <- lapply(enzyList, singleEnzymeSEM, input= input_trans, adj = adj_matrix, kmo.off = kmo.off, mdsite = mdsite)

    }else if(datatype == "ubi"){
      result <- lapply(enzyList, singleEnzymeSEM, input= input_trans, adj = adj_matrix, kmo.off = kmo.off, mdsite = FALSE)

    }else if(datatype == "ace"){
      result <- lapply(enzyList, singleEnzymeSEM, input= input_trans, adj = adj_matrix, kmo.off = kmo.off, mdsite = FALSE)
    }


    ### Raw activity score
    timestamp <- format(Sys.time(), "%Y%m%d.%H%M%S")

    activityRaw <- catchActRaw(result)
    write.table(activityRaw, file = file.path(output.folder, paste0("Activity_", timestamp, ".txt")),
                quote = FALSE, sep = "\t", row.names = FALSE)
    ### Mean-Center score
    # activityMean <- catchActMeanCenter(result)
    # write.table(activityMean, file = file.path(output.folder, paste0("Activity_MeanCenter_", timestamp, ".txt")),
    #             quote = FALSE, sep = "\t", row.names = FALSE)

    ### Z-score
    # activityZscore <- catchActZscore(result)
    # write.table(activityZscore, file = file.path(output.folder, paste0("Activity_Zscore_", timestamp, ".txt")),
    #             quote = FALSE, sep = "\t", row.names = FALSE)


    ### AFFINITY
    affinity <- catchAff(result_list = result, rawData = input)
    ##
    write.table(affinity, file = file.path(output.folder, paste0("Affinity_", timestamp, ".txt")),
                quote = FALSE, sep = "\t", row.names = FALSE)

    ### Evaluations
    eval_table <- catchEvaluations(result_list = result)
    eval_table <- addEvalDescription(eval_table = eval_table)
    message("\n", "Calculation finished! Start gererating...")
    ##
    # Define output path
    eval_path <- file.path(output.folder, paste0("Eval_table_", timestamp, ".txt"))

    # Write header description first
    writeLines(attr(eval_table, "Description"), eval_path)

    # Then append the actual evaluation table
    suppressWarnings(write.table(eval_table,
                                 file = eval_path,
                                 quote = FALSE,
                                 sep = "\t",
                                 row.names = FALSE,
                                 append = TRUE))
    message("\n", "Done!!!")

    # combine activity and affinity
    out <- list(activityRaw, affinity, eval_table)
    names(out) <- c("Activity", "Affinity", "Evaluations")

    return(out)
  }, silent= T)
}









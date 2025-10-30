#' @export getFa
#' @export subKinase
#' @export getSeq
#' @title Build MEME input
#' @description script to extract MEME input .fa file
#' @param ss Vector; Substrate species; ("mouse", "rat", "human")
#' @param es Vector; Kinase species; ("mouse", "rat", "human")
#' @param filePath The path is used to save .fa files. Recommend to create a new folder to save all .fa.files.
#' @examples getFa(ss = "human", es = c("human", "mouse", "rat"), filePath = getwd())
#'
#'
#'
#'
#'
#'



getFa <- function(ss = NULL,
                  es = NULL,
                  filePath){
  if (is.null(ss)) {
    stop("Please set your substrate organism!")
  }else{
    if (is.null(es)){
      es <- ss
    }else{
      species_kinase <- getSeq(ens = es, subs = ss)
      kinase_list <- levels(factor(species_kinase$GENE))
      lapply(kinase_list, subKinase, psp_table = species_kinase, savePath = filePath)

    }
      }
}


getSeq <- function(ens, subs){
  # Load psp
  psp <- psp_data()
  # Toupper gene name
  psp$GENE <- toupper(psp$GENE)
  # Replace "-" with "." in SUB_ACC_ID
  psp$SUB_ACC_ID <- as.matrix(gsub("-",".",as.matrix(psp$SUB_ACC_ID)))
  # Extract ORGANISM psp
  psp_species <- psp[which(psp$KIN_ORGANISM %in% ens & psp$SUB_ORGANISM %in% subs), c(1, 7, 10, 12)]

  # Replace "_" with "" in peptide (remove empty amino acid sites)
  psp_species$SITE_...7_AA <- as.matrix(gsub("_", "",as.matrix(psp_species$SITE_...7_AA)))

  # Modify PSP_human
  psp_species$psite <- paste(psp_species$SUB_ACC_ID, psp_species$SUB_MOD_RSD, sep = "_")
  psp_species$nchar <- nchar(psp_species$SITE_...7_AA)

  psp_species$site_name <- paste(">", psp_species$psite, sep = "")
  psp_species$seq <- paste(psp_species$site_name, "\n", psp_species$SITE_...7_AA)
  return(psp_species)
}


subKinase <- function(kinase, psp_table, savePath){
  kinase_table <- psp_table[which(psp_table$GENE == kinase), ]
  kinase_file <- paste(kinase_table$seq)
  kinase_file <- gsub(" ","", kinase_file)
  write.table(kinase_file, file = file.path(savePath, paste0(kinase, ".fa", sep = "")),
              row.names = F, col.names = F, quote = F)
}

##

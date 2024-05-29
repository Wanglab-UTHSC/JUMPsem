#' @docType data
#' @name DataSet
#'
#' @aliases psp
#' @aliases ubi
#' @aliases ace
#' @aliases input_psp_example
#' @aliases input_ubi_example
#' @aliases input_ace_example
#' @aliases wholeProteome_psp_example
#' @aliases motif_example
#' @keywords datasets
NULL


#' @export psp_data
#' @export ubi_data
#' @export ace_data
#' @export inputExample_PSPdata
#' @export inputExample_UBIdata
#' @export inputExample_ACEdata
#' @export wholeProteomeExample_data
#' @export motifExample_data
#'
#'
#'
#'

psp_data <- function() {
  utils::data(list="PhosphositePlus", package="JUMPsem")
  get("psp", envir = .GlobalEnv)
}

ubi_data <- function() {
  utils::data(list="UbiquitinDatabase", package="JUMPsem")
  get("ubi", envir = .GlobalEnv)
}


ace_data <- function(){
  utils::data(list="AcetyDatabase", package="JUMPsem")
  get("ace", envir = .GlobalEnv)
}





inputExample_PSPdata <- function() {
  utils::data(list="input_psp_example", package="JUMPsem")
  get("input_psp_example", envir = .GlobalEnv)
}

inputExample_UBIdata <- function() {
  utils::data(list="input_ubi_example", package="JUMPsem")
  get("input_ubi_example", envir = .GlobalEnv)
}


inputExample_ACEdata <- function() {
  utils::data(list="input_ace_example", package="JUMPsem")
  get("input_ace_example", envir = .GlobalEnv)
}



wholeProteomeExample_data <- function() {
  utils::data(list="wholeProteome_psp_example", package="JUMPsem")
  get("wholeProteome_psp_example", envir = .GlobalEnv)
}


motifExample_data <- function() {
  utils::data(list="motif_example", package="JUMPsem")
  get("motif_example", envir = .GlobalEnv)
}


# psp <- read.csv("/Volumes/LaCie/Project_Dehui_2022/KSEM/02_KSEM_code/eSEM_v3/eSEM/data/PhosphositePlus.csv", fileEncoding = "UTF-8-BOM")
# save_path <- file.path(".", "data", "PhosphositePlus.RData")
# Save the object using the file path
# save(psp, file = save_path)


#ubi <- read.csv("/Volumes/LaCie/Project_Dehui_2022/KSEM/02_KSEM_code/eSEM_v1/eSEM/data/UbiquitinDatabase.csv")
#save_path <- file.path(".", "data", "UbiquitinDatabase.RData")
#save(ubi, file = save_path)

#ace <- read.csv("/Volumes/LaCie/Project_Dehui_2022/KSEM/02_KSEM_code/eSEM_v1/eSEM/data/AcetyDatabase.csv")
#save_path <- file.path(".", "data", "AcetyDatabase.RData")
#save(ace, file = save_path)


#input_ubi_example <- read.csv("/Volumes/LaCie/Project_Dehui_2022/KSEM/02_KSEM_code/eSEM_v1/eSEM/data/input_ubi_example.csv")
#input_ubi_example <- input_ubi_example[, 1:15]
#save_path <- file.path(".", "data", "input_ubi_example.RData")
#save(input_ubi_example, file = save_path)


# input_psp_example <- read.csv("/Volumes/LaCie/Project_Dehui_2022/KSEM/08_paper/Hong/input_deep_new.csv")
# input_psp_example <- read.csv("/Volumes/LaCie/Project_Dehui_2022/KSEM/13_RShiny/01_exampleData/input_psp_example_trim.csv")
# colnames(input_psp_example)[3] <- "sites"
# save_path <- file.path(".", "data", "input_psp_example.RData")
# save(input_psp_example, file = save_path)


# input_ubi_example <- read.csv("/Volumes/LaCie/Project_Dehui_2022/KSEM/02_KSEM_code/eSEM_v1/eSEM/data/input_ubi_example.csv")
# input_ubi_example <- input_ubi_example[, 1:15]
# save_path <- file.path(".", "data", "input_ubi_example.RData")
# save(input_ubi_example, file = save_path)


# input_ace_example <- read.csv("/Volumes/LaCie/Project_Dehui_2022/KSEM/02_KSEM_code/eSEM_v1/eSEM/data/input_ace_example_raw.csv", fileEncoding = "UTF-8-BOM")
# ace_set <- input_ace_example %>%
#  separate_rows(accession_numbers, sep = "\\|")
# write.csv(ace_set, file = "/Volumes/LaCie/Project_Dehui_2022/KSEM/02_KSEM_code/eSEM_v1/eSEM/data/input_ace_example.csv", row.names = F)
# aceMap <- read.csv("/Volumes/LaCie/Project_Dehui_2022/KSEM/02_KSEM_code/eSEM_v1/eSEM/data/ace_mappingLibrary.csv")
# ace <- left_join(ace_set, aceMap, by = "accession_numbers")
# ace$sites <- "NA"
# ace$Peptides <- "NA"
# ace_example <- ace[, c(11, 13, 14, 15, 2:10)]
# input_ace_example <- data.frame(ace_example[which(ace_example$Substrate_ACC != "NA"), ])
# write.csv(input_ace_example, file = "/Volumes/LaCie/Project_Dehui_2022/KSEM/02_KSEM_code/eSEM_v1/eSEM/data/input_ace_example.csv", row.names = F)
# save(input_ace_example, file = "./data/input_ace_example.RData")


# wholeProteome_example <- read.csv("/Volumes/LaCie/Project_Dehui_2022/KSEM/02_KSEM_code/eSEM_v1/eSEM/data/wholeProteome_psp_example.csv")
# save_path <- file.path(".", "data", "wholeProteome_psp_example.RData")
# save(wholeProteome_example, file = save_path)

# mastFile <- read.table("C:\\Project_Dehui_2022\\KSEM\\02_KSEM_code\\KSEM_v2\\KSEM\\data_raw\\extract_result.txt")
# motif_example <- read.csv("/Volumes/LaCie/Project_Dehui_2022/KSEM/08_paper/Hong/motif_input_15.csv")
# colnames(motif_example) <- c("GENE", "psite")
# save(motif_example, file = "./data/motif_example.RData")



######## Trim data for RShiny
# input_psp_example <- read.csv("/Volumes/LaCie/Project_Dehui_2022/KSEM/13_RShiny/01_exampleData/input_psp_example_trim.csv")
# pspDB <- read.csv("/Volumes/LaCie/Project_Dehui_2022/KSEM/02_KSEM_code/eSEM_v3/eSEM/data/PhosphositePlus.csv")
# input_psp_example$psite <- paste0(input_psp_example$Substrate_ACC, "_", input_psp_example$sites, sep = "")
# pspDB$psite <- paste0(pspDB$SUB_ACC_ID, "_", pspDB$SUB_MOD_RSD, sep = "")
# pspDBf <- subset(pspDB,
#                 SUB_ORGANISM %in% "mouse" &
#                    KIN_ORGANISM %in% c("mouse", "human", "rat") &
#                    KINASE %in% c("PRKG2", "PIM1", "RPS6KA5", "AURKB",
#                                  "GSK3A", "CAMK2A","DYRK1A", "PKN1",
#                                  "LYN", "CDK1"))
# input_psp_example_T <- subset(input_psp_example, psite %in% pspDBf$psite)
# input_psp_example_T <- input_psp_example_T[, 1: (ncol(input_psp_example_T) - 1)]
# save(input_psp_example_T, file = "./data/input_psp_example_T.RData")
# write.csv(input_psp_example_T, file = "/Volumes/LaCie/Project_Dehui_2022/KSEM/02_KSEM_code/eSEM_v3/eSEM/data/input_psp_example_T.csv", row.names = F)

# input_ubi_example <- read.csv("/Volumes/LaCie/Project_Dehui_2022/KSEM/02_KSEM_code/eSEM_v3/eSEM/data/input_ubi_example.csv")
# ubiDB <- read.csv("/Volumes/LaCie/Project_Dehui_2022/KSEM/02_KSEM_code/eSEM_v3/eSEM/data/UbiquitinDatabase.csv")
# ubiDBf <- subset(ubiDB, Gene.Symbol..E3. %in% c("PRPF19", "CCNF", "PPARG", "FBXO32", "CBLC",
#                                                 "AMFR", "GAN", "RLIM", "BIRC2", "MARCHF5"))
# input_ubi_example_T <-subset(input_ubi_example, Substrate_ACC %in% ubiDBf$SwissProt.AC..Substrate.)
# input_ubi_example_T <- input_ubi_example_T[, 1:15]
# colnames(input_ubi_example_T)[5:15] <- sapply(strsplit(colnames(input_ubi_example_T)[5:15], "_"), function(x) paste(x[-1], collapse = "_"))
# colnames(input_ubi_example_T)[5:15] <- sapply(strsplit(colnames(input_ubi_example_T)[5:15], "\\."), function(x) x[1])
# colnames(input_ubi_example_T)[14:15] <- c("Recovery1_Btz", "Recovery2_Btz")
# save(input_ubi_example_T, file = "./data/input_ubi_example_T.RData")
# write.csv(input_ubi_example_T, file = "/Volumes/LaCie/Project_Dehui_2022/KSEM/02_KSEM_code/eSEM_v3/eSEM/data/input_ubi_example_T.csv", row.names = F)


# input_ace_example_T <- read.csv("/Volumes/LaCie/Project_Dehui_2022/KSEM/02_KSEM_code/eSEM_v3/eSEM/data/input_ace_example_T.csv")
# aceDB <- read.csv("/Volumes/LaCie/Project_Dehui_2022/KSEM/02_KSEM_code/eSEM_v3/eSEM/data/AcetyDatabase.csv")
# overlap <- intersect(input_ace_example_T$Substrate_ACC, aceDB$UniProt)
# input_ace_example_T <- subset(input_ace_example_T, Substrate_ACC %in% overlap)
# save(input_ace_example_T, file = "./data/input_ace_example_T.RData")
# write.csv(input_ace_example_T, file = "/Volumes/LaCie/Project_Dehui_2022/KSEM/02_KSEM_code/eSEM_v3/eSEM/data/input_ace_example_T.csv", row.names = F)






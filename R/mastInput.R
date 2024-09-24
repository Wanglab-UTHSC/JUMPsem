#' @export mastInput
#' @title Build MAST input
#' @description Script to extract MAST input FASTA file.
#' @param input The same quantitative phosphoproteomics data as eSEM() input file.
#' @param savePath The path to save built MAST input file. Default is the current directory.
#' @importFrom utils write.table
#' @examples mastInput(input = input_psp_example)
#'
#'



mastInput <- function(input,
                      savePath = getwd()){
  input$sequence_name <- paste(">", "sites.", input$Substrate_ACC, "_", input$sites, sep = "")
  input$seq <- paste(input$sequence_name, "\n", input$Peptide, sep = "")
  mast_input <- as.matrix(input$seq)
  file_path <- file.path(savePath, "MAST_input.fa")
  #file_path <- paste(savePath, "\\", "MAST_input", ".fa", sep = "")
  write.table(mast_input, file = file_path,
              row.names = F, col.names = F, quote = F)
  return(mast_input)
}

#input <- read_excel("C:\\Project_Dehui_2022\\KSEM\\02_KSEM_code\\KSEM_v2\\KSEM\\data_raw\\input_raw.xlsx")
#mastInput(input = input, savePath = "C:\\Project_Dehui_2022\\KSEM\\02_KSEM_code\\KSEM_v2\\KSEM\\data_raw")


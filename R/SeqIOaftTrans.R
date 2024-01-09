#' @export SeqIOaftTrans
#' @description
#' @param SeqIOoutput output from SeqIO.py; Without header
#' @param input_trans transformed file from SeqIOpreTrans; could use input_raw also if there is NO comma in "sites" coloumn.
#' @param filePath Character. Character vector of location to save files if desired. Default is current directory.
#' @import stringr
#' @import dplyr



SeqIOaftTrans <- function(SeqIOoutput,
                          input_trans,
                          filePath = getwd())
{
  SeqIOoutput$Substrate_ACC <- (str_split(SeqIOoutput$V1, '_', simplify = TRUE))[, 1]
  SeqIOoutput$sites <- (str_split(SeqIOoutput$V1, '_', simplify = TRUE))[, 2]
  SeqIOoutput <- SeqIOoutput[, c(3, 4, 2)]

  results <- left_join(SeqIOpreTrans, SeqIOoutput, by = c("Substrate_ACC", "sites"))
  results <- results[, c(1:3, ncol(results), 5:(ncol(results)-1))]
  colnames(results)[4] <- "Peptide"
  results <- results[!results$Peptide == "", ]
  write.csv(results, file =  file.path(filePath, "Newinput_SeqIO.csv"),
            quote = FALSE, row.names = F)
  return(results)
}

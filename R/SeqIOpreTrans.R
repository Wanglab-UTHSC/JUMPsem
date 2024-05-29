#' @export SeqIOpreTrans
#' @description SeqIOpreTrans
#' @param input_raw Input raw data (quantitative phosphoproteomics data)
#' @param filePath Character. Character vector of location to save files if desired. Default is current directory.
#' @import dplyr
#' @import tidyr




SeqIOpreTrans <- function(input_raw,
                          filePath = getwd())
{
  ## separate rows with "," in "sites"
  input_1 <- input_raw %>%
    separate_rows(sites, sep = ",")
  ## delete rows with empty sites
  input_2 <- input_1[!(input_1$sites == "" | is.na(input_1$sites)), ]
  write.csv(input_2, file = file.path(filePath, "input_SeqIOpreTrans.csv"),
            quote = FALSE, row.names = F)
  return(input_2)
}


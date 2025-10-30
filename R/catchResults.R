#' @import dplyr
#' @import stringr
#' @import data.table
#' @import tidyr
#' @param result_list The result list from modeling fitting.
#' @param rawData  Input raw data (quantitative phospho-proteomics / ubiqutin-proteomics data)
#' @param eval_table eval_table A tibble or data.frame produced by `catchEvaluations()`.
#' @title catchResults
#' @description Catch Activity and Affinity.
#' @rdname catchResults
#'
#' @export catchActRaw

#' @export catchAff

#' @export catchEvaluations
#' @export addEvalDescription


# -------------------------------------------------------
# 1. Mean-Center Normalization
# -------------------------------------------------------

# centerApply <- function(mat) {
#   apply(mat, 2, function(y) y / mean(abs(y), na.rm = TRUE))
# }
#
# catchActMeanCenter <- function(result_list) {
#   valid_list <- result_list[sapply(result_list, function(x) {
#     is.list(x) && length(x) > 0 && !is.null(x[[1]])
#   })]
#
#   if (length(valid_list) == 0) {
#     stop("No valid enzyme activity entries found in input list.")
#   }
#
#   activity_mats <- lapply(valid_list, function(x) as.data.frame(x[[1]]))
#   activity_df <- do.call(cbind, activity_mats)
#   activity_centered <- centerApply(activity_df)
#
#   activity_t <- data.frame(t(activity_centered))
#   rownames(activity_t) <- names(valid_list)
#
#   activity_t$ENZYME <- rownames(activity_t)
#   final <- activity_t[, c(ncol(activity_t), 1:(ncol(activity_t) - 1))]
#
#   return(final)
# }

# -------------------------------------------------------
# 2. Z-score Normalization
# -------------------------------------------------------

# catchActZscore <- function(result_list) {
#   valid_list <- result_list[sapply(result_list, function(x) {
#     is.list(x) && length(x) > 0 && !is.null(x[[1]])
#   })]
#
#   if (length(valid_list) == 0) {
#     stop("No valid enzyme activity entries found in input list.")
#   }
#
#   activity_mats <- lapply(valid_list, function(x) as.data.frame(x[[1]]))
#   activity_df <- do.call(cbind, activity_mats)
#   activity_scaled <- scale(activity_df)
#
#   activity_t <- data.frame(t(activity_scaled))
#   rownames(activity_t) <- names(valid_list)
#
#   activity_t$ENZYME <- rownames(activity_t)
#   final <- activity_t[, c(ncol(activity_t), 1:(ncol(activity_t) - 1))]
#
#   return(final)
# }

# -------------------------------------------------------
# 3. Raw Activity (for consistency)
# -------------------------------------------------------

catchActRaw <- function(result_list) {
  # Keep only valid enzyme entries that contain a list with ≥1 non-null element
  valid_list <- result_list[sapply(result_list, function(x) {
    is.list(x) && length(x) > 0 && !is.null(x[[1]])
  })]

  if (length(valid_list) == 0) {
    stop("No valid enzyme activity entries found in input list.")
  }

  # Extract first element (kinase, ligase, etc.) automatically
  activity_mats <- lapply(valid_list, function(x) as.data.frame(x[[1]]))
  activity_df <- do.call(cbind, activity_mats)

  # Transpose (enzymes as rows)
  activity_t <- data.frame(t(activity_df))
  rownames(activity_t) <- names(valid_list)

  # Add ENZYME column
  activity_t$ENZYME <- rownames(activity_t)
  final <- activity_t[, c(ncol(activity_t), 1:(ncol(activity_t) - 1))]

  return(final)
}




catchAff <- function(result_list, rawData) {
  library(data.table)
  library(dplyr)
  library(tidyr)

  # 1. Filter valid enzyme entries
  valid_list <- result_list[sapply(result_list, function(x) {
    is.list(x) && length(x) > 0 && !is.null(x[[1]])
  })]

  if (length(valid_list) == 0) {
    warning("No valid enzyme affinity data found in input list.")
    return(data.frame())
  }

  # 2. Combine all affinity data safely
  aff_df <- tryCatch(
    data.table::rbindlist(
      lapply(valid_list, function(x) {
        df2 <- x[[2]]
        # Case 1: nested structure like list($PRKCD = data.frame(...))
        if (is.list(df2) && length(df2) == 1 && is.data.frame(df2[[1]])) {
          return(as.data.frame(df2[[1]]))
        }
        # Case 2: direct data frame
        if (is.data.frame(df2) && nrow(df2) > 0) {
          return(as.data.frame(df2))
        }
        # Otherwise skip
        return(NULL)
      }),
      use.names = FALSE,
      fill = TRUE,
      idcol = "ENZYME"
    ),
    error = function(e) {
      warning("Affinity data combination failed: ", conditionMessage(e))
      return(data.frame())
    }
  )


  # 3. If empty, return blank result
  if (nrow(aff_df) == 0 || ncol(aff_df) < 2) {
    warning("No valid affinity data to process.")
    return(data.frame())
  }

  # 4. Split the substrate column safely
  aff_df <- tryCatch({
    separate(
      aff_df,
      col = 2,
      into = c("Substrate_ACC", "sites"),
      sep = "_",
      remove = FALSE
    )
  }, error = function(e) {
    warning("Substrate split failed: ", conditionMessage(e))
    aff_df$Substrate_ACC <- NA
    aff_df$sites <- NA
    return(aff_df)
  })

  # 5. Ensure rawData has usable structure
  if (ncol(rawData) >= 2) {
    subgene <- as.data.frame(rawData[, 1:2])
    colnames(subgene)[1:2] <- c("Substrate_ACC", "Substrate_GN")
    subgene$Substrate_ACC <- gsub("-", ".", subgene$Substrate_ACC)
  } else {
    warning("rawData does not have at least 2 columns; skipping merge.")
    subgene <- data.frame(Substrate_ACC = character(0), Substrate_GN = character(0))
  }

  # 6. Merge with gene info
  merged_df <- suppressWarnings(left_join(aff_df, subgene, by = "Substrate_ACC"))

  # 7. Rename and reorder columns safely
  merged_df <- merged_df %>%
    rename(
      SUB_ACC_ID = Substrate_ACC,
      SUB_sites = sites,
      ENZYME_GN = ENZYME,
      SUBSTRATE_GN = Substrate_GN,
      AFFINITY = affinity
    )

  expected_cols <- c("ENZYME_GN", "SUB_ACC_ID", "SUBSTRATE_GN",
                     "SUB_sites", "AFFINITY", "z", "Pvalue")

  merged_df <- merged_df %>%
    select(any_of(expected_cols))

  merged_df <- merged_df %>%
    distinct(ENZYME_GN, SUB_ACC_ID, SUBSTRATE_GN, SUB_sites, AFFINITY, z, Pvalue, .keep_all = TRUE)

  return(merged_df)
}


catchEvaluations <- function(result_list) {
  # Load optional tidyverse tools if available
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install 'dplyr'")
  if (!requireNamespace("purrr", quietly = TRUE)) stop("Please install 'purrr'")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Please install 'tibble'")

  library(dplyr)
  library(purrr)
  library(tibble)

  eval_table <- map_dfr(names(result_list), function(enzyme) {
    metrics <- try(result_list[[enzyme]][["evaluation"]][[enzyme]][["metrics"]], silent = TRUE)
    if (!inherits(metrics, "try-error") && !is.null(metrics)) {
      as_tibble_row(as.list(metrics)) %>%
        mutate(ENZYME = enzyme, .before = 1)
    } else {
      NULL
    }
  })

  colnames(eval_table) <- c("ENZYME", "CFI", "TLI", "RMSEA", "SRMR")
  return(eval_table)
}


addEvalDescription <- function(eval_table) {
  description_text <- paste(
    "Metric\tMeaning\tInterpretation",
    "CFI\tComparative Fit Index\tMeasures model fit relative to a null model; >= 0.90 indicates acceptable fit",
    "TLI\tTucker-Lewis Index\tAdjusts fit for model complexity; >= 0.90 preferred",
    "RMSEA\tRoot Mean Square Error of Approximation\tPenalizes poor fit per degree of freedom; =< 0.08 acceptable, =< 0.05 good",
    "SRMR\tStandardized Root Mean Square Residual\tAverage difference between observed and predicted correlations; =< 0.08 good",
    "",   # ← one blank line
    "",   # ← another blank line for spacing
    sep = "\n"
  )

  attr(eval_table, "Description") <- description_text
  return(eval_table)
}



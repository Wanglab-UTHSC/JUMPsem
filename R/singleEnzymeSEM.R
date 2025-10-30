#' @description
#' Estimate the latent activity of a single enzyme using Structural Equation Modeling (SEM).
#' This function identifies substrates linked to an enzyme, filters based on data quality,
#' and infers latent enzyme activity using lavaan. It also evaluates the model fit
#' (CFI, TLI, RMSEA, SRMR) and provides a descriptive report.
#'
#' @title singleEnzymeSEM
#' @export
#' @import lavaan
#' @import EFAtools
#' @import Matrix
#'
#' @param enzyme Character. Name of the enzyme (must match column name in adjacency matrix).
#' @param input Data frame. Normalized and transformed phosphoproteomics data (rows = sites, cols = samples).
#' @param adj Matrix. Enzyme–substrate adjacency matrix (rows = substrates, columns = enzymes).
#' @param kmo.off Numeric. Cutoff for KMO (Kaiser-Meyer-Olkin) measure; default = 0.
#' @param mdsite Logical. If TRUE, match substrates using exact PTM site IDs (e.g., "CDK16_S12").
#'                If FALSE, match using gene-level mapping (e.g., "CDK16").
#'
#' @return A list with three components for the given enzyme:
#'   \describe{
#'     \item{activity}{Latent enzyme activity scores (sample × 1 matrix)}
#'     \item{affinity}{Table of substrate loadings, z-values, and p-values}
#'     \item{evaluation}{Model fit indices and descriptive interpretation}
#'   }


singleEnzymeSEM <- function(enzyme,
                            input,
                            adj,
                            kmo.off = 0,
                            mdsite = TRUE) {
  try({
    # Initialize output container ---------------------------------------------
    result_list <- list(
      activity = list(),   # latent scores per enzyme
      affinity = list(),   # substrate loadings (affinities)
      evaluation = list()  # model fit metrics
    )

    # Step 1. Identify enzyme substrates from adjacency matrix ----------------
    enzyme_substrates <- rownames(adj[which(adj[, enzyme] == 1), ])
    rownames(input) <- gsub("-", ".", rownames(input))  # sanitize names for modeling

    # Step 2. Find intersection between enzyme’s known substrates and available data -----
    if (mdsite) {
      # Match exact PTM identifiers (e.g., "Q9Y6R4_S123")
      intersect_substrates <- intersect(x = enzyme_substrates, y = rownames(input))
    } else {
      # Match by gene symbol (e.g., "CDK16" portion from "CDK16_S12")
      siteN <- rownames(input)
      first_parts <- sapply(strsplit(siteN, "_"), function(x) x[1])
      intersect_substrates <- intersect(x = enzyme_substrates, y = first_parts)
    }

    # Skip if no valid overlap between enzyme–substrate mappings ----------------
    if (length(intersect_substrates) > 1) {

      # Step 3. Prepare data subset for this enzyme ----------------------------
      data <- data.frame(input)

      if (mdsite) {
        # Directly subset by matching rownames (exact site match)
        enzyme_table_raw <- t(data[rownames(data) %in% intersect_substrates, ])
      } else {
        # Subset based on gene names only
        split_row_names <- sapply(rownames(data), function(x) unlist(strsplit(x, "_"))[1])
        enzyme_table_raw <- t(data[split_row_names %in% intersect_substrates, ])
        enzyme_table_raw <- as.data.frame(apply(enzyme_table_raw, 2, as.numeric))
      }

      enzyme_table_new <- enzyme_table_raw

      # Step 4. Compute correlation matrix for KMO filtering -------------------
      cor_new <- cor(enzyme_table_new, use = "pairwise.complete.obs")

      # Step 5. Apply nearPD and KMO filter to remove unstable variables -------
      if (length(colnames(enzyme_table_new)) > 2) {
        # Check matrix condition; repair if ill-conditioned
        recip_cond <- 1 / kappa(cor_new)
        if (!is.finite(recip_cond) || recip_cond < 1e-12) {
          cor_new <- as.matrix(nearPD(cor_new, corr = TRUE, maxit = 500, conv.tol = 1e-6)$mat)
        }
        # KMO filter (keep columns with KMO > cutoff)
        kmo_result <- KMO(cor_new)
        enzyme_table_md <- enzyme_table_new[, which(as.matrix(kmo_result[[2]]) > kmo.off), drop = FALSE]
      }else {
        enzyme_table_md <- enzyme_table_new
      }

      # Step 6. SEM modeling ---------------------------------------------------
      if (length(colnames(enzyme_table_md)) > 1) {

        # Transpose to get (sample × substrate) structure
        enzyme_new_t <- as.matrix(t(enzyme_table_md))
        enzyme_new_t_rowname <- rownames(enzyme_new_t)

        # Build model syntax: enzyme =~ substrate1 + substrate2 + ...
        substrate <- noquote(paste(enzyme_new_t_rowname, collapse = " + "))
        enzyme_model <- paste(enzyme, "=~", substrate)

        # Extract corresponding data for SEM
        data_s <- subset(data, rownames(data) %in% enzyme_new_t_rowname)
        data_sem <- t(data_s)

        # Get number of samples and variables
        n_samples <- nrow(data_sem)
        n_vars <- ncol(data_sem)

        # Step 6a. Dimensionality reduction (optional PCA if p > n) ------------
        if (n_vars >= n_samples) {
          pca <- prcomp(data_sem, scale. = TRUE)
          eigvals <- pca$sdev^2
          pca_scores <- as.data.frame(pca$x)
          pc1 <- pca_scores[, 1]

          # Select substrates most correlated with PC1
          cor_data <- cor(pc1, data_sem)
          cor_rank <- sort(cor_data[1, ], decreasing = TRUE)
          top_n <- min(n_samples - 1, length(cor_rank))
          selected_substrates <- names(cor_rank)[1:top_n]

          data_reduced <- as.data.frame(data_sem[, selected_substrates, drop = FALSE])
        } else {
          data_reduced <- as.data.frame(data_sem)
        }

        # Update model syntax based on reduced features
        substrate <- unique(colnames(data_reduced))
        enzyme_model <- paste(enzyme, "=~", paste(substrate, collapse = " + "))

        # Step 6b. Fit SEM using lavaan ----------------------------------------
        fit_sem <- suppressWarnings(sem(enzyme_model, data = data_sem))

        # Step 6c. Evaluate modification indices -------------------------------

        mi <- tryCatch(
          suppressWarnings(modindices(fit_sem)),
          error = function(e) NULL
        )

        if (!is.null(mi) && is.data.frame(mi) && "mi" %in% colnames(mi)) {
          mi_sorted <- mi[order(mi$mi, decreasing = TRUE), ]
          high_mi <- subset(mi_sorted, mi > 4 & op == "~~") # Only covariances with MI > 4
        } else {
          mi_sorted <- data.frame()
          high_mi <- data.frame()
        }

        # Refine model silently based on high MI pairs
        if (nrow(high_mi) > 0) {
          # Build additional covariance paths
          added_paths <- sprintf("%s ~~ %s", high_mi[["lhs"]], high_mi[["rhs"]])
          extra_model <- paste(added_paths, collapse = "\n")

          # Combine original model with new paths
          enzyme_model_refined <- paste(enzyme_model, extra_model, sep = "\n")

          # Refit model
          fit_sem <- suppressWarnings(sem(enzyme_model_refined, data = data_sem))
          refinement_status <- sprintf("Model refined by adding %d covariance paths.", nrow(high_mi))
        } else {
          refinement_status <- "No significant MI (>4); kept original model."
        }

        # Step 7. Extract latent scores ----------------------------------------
        scores <- as.matrix(lavPredict(fit_sem))
        rownames(scores) <- rownames(data_sem)

        # Step 8. Extract factor loadings (affinity table) ---------------------
        summary_PE <- tryCatch(
          quiet(summary(fit_sem)[[5]]),  # suppress verbose output
          error = function(e) NULL
        )

        if (!is.null(summary_PE)) {
          main_table <- summary_PE[summary_PE$op == "=~", c(3, 5, 7, 8)]
          colnames(main_table) <- c("substrate_sites", "affinity", "z", "Pvalue")
          main_table$affinity <- abs(main_table$affinity)  # take absolute loading
        } else {
          main_table <- data.frame()
        }

        # Step 9. Model fit evaluation metrics ---------------------------------
        eval_metrics <- tryCatch({
          fitMeasures(fit_sem, c("cfi", "tli", "rmsea", "srmr"))
        }, error = function(e) NULL)

        # Step 10. Interpret model fit quality ---------------------------------
        interpret_fit <- function(m) {
          if (is.null(m)) return("Model fit could not be evaluated.")

          # Helper function to classify each metric
          interpret_value <- function(value, good, ok, metric_name, direction = "higher") {
            if (is.na(value)) return(paste(metric_name, "NA"))
            if (direction == "higher") {
              if (value >= good) return("good")
              else if (value >= ok) return("acceptable")
              else return("poor")
            } else {
              if (value <= good) return("good")
              else if (value <= ok) return("acceptable")
              else return("poor")
            }
          }

          # Interpret individual indices
          cfi_desc  <- interpret_value(m["cfi"], 0.95, 0.90, "cfi", "higher")
          tli_desc  <- interpret_value(m["tli"], 0.95, 0.90, "tli", "higher")
          rmsea_desc <- interpret_value(m["rmsea"], 0.05, 0.08, "rmsea", "lower")
          srmr_desc <- interpret_value(m["srmr"], 0.05, 0.10, "srmr", "lower")

          # Combine into a descriptive summary message
          message <- paste0(
            sprintf(
              "[Model Fit for %s]\n  CFI = %.3f (%s fit)\n  TLI = %.3f (%s fit)\n  RMSEA = %.3f (%s fit)\n  SRMR = %.3f (%s fit)\n",
              enzyme, m["cfi"], cfi_desc, m["tli"], tli_desc,
              m["rmsea"], rmsea_desc, m["srmr"], srmr_desc
            )
          )

          # Overall interpretation summary
          overall_quality <- if (all(c(cfi_desc, tli_desc, rmsea_desc, srmr_desc) == "good")) {
            "Excellent model fit."
          } else if (any(cfi_desc == "poor", tli_desc == "poor", rmsea_desc == "poor", srmr_desc == "poor")) {
            "Model fit is poor; interpret latent activity cautiously."
          } else {
            "Model fit is acceptable."
          }

          message <- paste0(message, "Overall: ", overall_quality)
          return(message)
        }

        # Generate message and print to console
        eval_message <- interpret_fit(eval_metrics)
        cat(eval_message, "\n")

        # Step 11. Save all results --------------------------------------------
        result_list$activity[[enzyme]] <- scores
        result_list$affinity[[enzyme]] <- main_table
        result_list$evaluation[[enzyme]] <- list(
          metrics = eval_metrics,
          message = eval_message
        )

        print(paste0("Finished SEM for enzyme: ", enzyme))
        return(result_list)
      }
    }
  }, silent = TRUE)
}

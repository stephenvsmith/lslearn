#' SNL Algorithm
#'
#' `snl()` applies the SNL (baseline local structure-learning) algorithm to
#' a vector of target nodes in a network structure and returns an estimated
#' PDAG around their neighborhoods.
#'
#' Exactly one of `data` or `true_dag` (population/oracle version) must
#' ultimately be usable to identify each target's neighborhood: if
#' `true_dag` is `NULL`, the neighborhoods are first estimated from `data`
#' via `get_all_mbs()`; otherwise `true_dag` is used directly (with `data`,
#' if also supplied, used only for the conditional independence tests
#' themselves -- the "semi-sample" version).
#'
#' @inheritParams cml
#'
#' @returns A list with the estimated PDAG adjacency matrix (`amat`),
#'   separating sets (`S`), the number of conditional independence tests
#'   used (`NumTests`, `MBNumTests`), the nodes considered (`Nodes`, in
#'   1-based R numbering), timing information, the reference DAG used
#'   (`referenceDAG`), the estimated Markov Blanket list (`mbList`, if
#'   applicable), which orientation rules were used (`rules_used`), and,
#'   when `data` was supplied, the pre-scaling column means/covariance
#'   (`data_means`, `data_cov`).
#' @export
snl <- function(data = NULL, true_dag = NULL, targets,
                node_names = NULL, lmax = 3, tol = 0.01, mb_tol = 0.05,
                method = "MMPC", test = "testIndFisher", verbose = TRUE) {
  if (lmax < 0) {
    stop("Invalid lmax value")
  }

  p <- ifelse(is.null(data), ncol(true_dag), ncol(data))
  if (is.null(node_names)) {
    node_names <- paste0("V", 0:(p - 1))
  }

  data_means <- NA
  data_cov <- NA
  if (!is.null(data)) {
    if (is.data.frame(data)) {
      data <- as.matrix(data)
    }
    # Store data information
    data_means <- colMeans(data)
    data_cov <- stats::cov(data)

    # Scale the data
    data <- scale(data)
  }
  mb_num_tests <- 0
  if (is.null(true_dag)) {
    # Find Markov Blankets (mb-estimation.R)
    result <- get_all_mbs(targets, data, mb_tol, lmax, method, test, verbose)
    mb_list <- result$mb_list
    nodes_interest <- as.numeric(names(mb_list)) - 1
    mb_num_tests <- result$num_tests
    # Create adjacency matrix based on Markov Blankets (mb-estimation.R)
    true_dag <- get_est_initial_dag(mb_list, ncol(data), verbose)
    semi_sample_version <- FALSE
    est_dag <- TRUE
  } else {
    mb_list <- list()
    semi_sample_version <- TRUE
    nodes_interest <- seq(0, p - 1)
    est_dag <- FALSE
  }

  # Convert any data frame to a matrix
  if (is.data.frame(true_dag)) {
    true_dag <- as.matrix(true_dag)
  }

  # To account for zero-indexing in C++
  cpp_targets <- targets - 1
  if (verbose) {
    cat(
      "The node value for the C++ function is",
      paste(cpp_targets, collapse = ","),
      "\n"
    )
  }

  if (is.null(data)) {
    if (verbose) {
      cat("Population Version:\n")
    }
    results <- popSNL(
      true_dag, cpp_targets, nodes_interest,
      node_names, lmax, verbose
    )
  } else {
    if (semi_sample_version && verbose) {
      cat("Semi-Sample Version:\n")
    } else {
      if (verbose) cat("Sample Version:\n")
    }

    results <- sampleSNL(
      true_dag, data, cpp_targets, nodes_interest,
      node_names, lmax, tol, verbose, test, est_dag
    )
  }

  # We change the target to target - 1 in order to accommodate change to C++
  list(
    "amat" = results$G,
    "S" = results$S,
    "NumTests" = results$NumTests + mb_num_tests,
    "MBNumTests" = mb_num_tests,
    "Nodes" = results$allNodes + 1, # to convert to R numbering
    "targetSkeletonTimes" = results$targetSkeletonTimes,
    "totalTime" = results$totalTime,
    "referenceDAG" = true_dag,
    "mbList" = mb_list,
    "rules_used" = results$rulesUsed,
    "data_means" = data_means,
    "data_cov" = data_cov
  )
}

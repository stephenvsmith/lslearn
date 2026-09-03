#' CML Algorithm
#'
#' `cml()` applies the CML (Coordinated Multi-Neighborhood Learning)
#' algorithm to a dataset over certain target neighborhoods, which may be
#' provided by the user through the true DAG or may be estimated by a
#' Markov Blanket estimation algorithm.
#'
#' Exactly one of `data` or `true_dag` (population/oracle version) must
#' ultimately be usable to identify each target's neighborhood: if
#' `true_dag` is `NULL`, the neighborhoods are first estimated from `data`
#' via `get_all_mbs()`; otherwise `true_dag` is used directly (with `data`,
#' if also supplied, used only for the conditional independence tests
#' themselves -- the "semi-sample" version).
#'
#' @param data A data matrix or data frame, or `NULL` for the population
#'   (oracle) version.
#' @param true_dag A 0/1 adjacency matrix of the true DAG, or `NULL` to
#'   estimate the neighborhoods from `data` instead.
#' @param targets A vector of 1-based target node indices.
#' @param node_names Character vector of node names; defaults to `V0`,
#'   `V1`, ... if not supplied.
#' @param lmax Maximum size of the conditioning set considered during the
#'   skeleton search.
#' @param tol Significance level for the conditional independence tests
#'   used by the skeleton search.
#' @param mb_tol Significance level used for Markov Blanket estimation (see
#'   `get_all_mbs()`); only relevant when `true_dag` is `NULL`.
#' @param method Markov Blanket estimation algorithm (see `get_mb()`); only
#'   relevant when `true_dag` is `NULL`.
#' @param test The conditional independence test to use.
#' @param verbose Whether to provide detailed output.
#'
#' @returns A list with the estimated adjacency matrix (`amat`), separating
#'   sets (`S`), the number of conditional independence tests used
#'   (`NumTests`, `MBNumTests`), the nodes considered (`Nodes`, in 1-based R
#'   numbering), which FCI orientation rules were used (`RulesUsed`), timing
#'   information, the reference DAG used (`referenceDAG`), the estimated
#'   Markov Blanket list (`mbList`, if applicable), and, when `data` was
#'   supplied, the pre-scaling column means/covariance (`data_means`,
#'   `data_cov`).
#' @export
cml <- function(data = NULL, true_dag = NULL, targets,
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
    # Store data statistics
    data_means <- colMeans(data)
    data_cov <- stats::cov(data)

    # Standardize the data (for only continuous data)
    if (test == "testIndFisher") {
      data <- scale(data)
    } else if (test == "gSquare") {
      data <- data - min(data)
    }
  }
  mb_num_tests <- 0
  mb_time_track <- NA
  if (is.null(true_dag)) {
    # Find Markov Blankets (mb-estimation.R)
    mb_start <- Sys.time()
    result <- get_all_mbs(targets, data, mb_tol, lmax, method, test, verbose)
    mb_end <- Sys.time()
    mb_diff <- mb_end - mb_start
    units(mb_diff) <- "secs"

    # Track the time needed to estimate MB's
    mb_time_track <- as.numeric(mb_diff)
    mb_list <- result$mb_list

    # Store all of the nodes estimated to be needed (targets and
    # first-order neighbors). Subtract 1 for conversion from C++ numbering
    # scheme
    nodes_interest <- as.numeric(names(mb_list)) - 1
    mb_num_tests <- result$num_tests

    # Create adjacency matrix based on Markov Blankets (mb-estimation.R)
    true_dag <- get_est_initial_dag(mb_list, p, verbose)
    # We are using a DAG with estimated Markov Blankets encoded
    est_dag <- TRUE
  } else {
    mb_list <- list()
    nodes_interest <- seq(0, p - 1)
    est_dag <- FALSE
  }

  # Convert any data.frame to a matrix
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
    results <- popCML(
      true_dag, cpp_targets, nodes_interest,
      node_names, lmax, verbose
    )
  } else {
    results <- sampleCML(
      true_dag, data, cpp_targets, nodes_interest,
      node_names, lmax, tol, verbose, test, est_dag
    )
  }

  # We change the target to target - 1 in order to accommodate change to C++
  list(
    "amat" = results$G,
    "S" = results$S,
    "NumTests" = results$NumTests,
    "MBNumTests" = mb_num_tests,
    "RulesUsed" = results$RulesUsed,
    "Nodes" = results$allNodes + 1, # to convert to R numbering
    "mbEstTime" = mb_time_track,
    "totalSkeletonTime" = results$totalSkeletonTime,
    "targetSkeletonTimes" = results$targetSkeletonTimes,
    "totalTime" = results$totalTime,
    "referenceDAG" = true_dag,
    "mbList" = mb_list,
    "data_means" = data_means,
    "data_cov" = data_cov
  )
}

#' CML Algorithm (MAG-only output)
#'
#' `cml_mag()` applies the CML algorithm identically to [cml()], except it
#' skips the within-neighborhood mixed-graph conversion rules, leaving the
#' result as an ancestral (MAG-style) graph rather than converting it to
#' CML's own neighborhood notation.
#'
#' @inheritParams cml
#'
#' @returns Same shape as [cml()].
#' @export
cml_mag <- function(data = NULL, true_dag = NULL, targets,
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
    # Store data statistics
    data_means <- colMeans(data)
    data_cov <- stats::cov(data)

    # Standardize the data (for only continuous data)
    if (test == "testIndFisher") {
      data <- scale(data)
    } else if (test == "gSquare") {
      data <- data - min(data)
    }
  }
  mb_num_tests <- 0
  mb_time_track <- NA
  if (is.null(true_dag)) {
    # Find Markov Blankets (mb-estimation.R)
    mb_start <- Sys.time()
    result <- get_all_mbs(targets, data, mb_tol, lmax, method, test, verbose)
    mb_end <- Sys.time()
    mb_diff <- mb_end - mb_start
    units(mb_diff) <- "secs"

    # Track the time needed to estimate MB's
    mb_time_track <- as.numeric(mb_diff)
    mb_list <- result$mb_list

    # Store all of the nodes estimated to be needed (targets and
    # first-order neighbors). Subtract 1 for conversion from C++ numbering
    # scheme
    nodes_interest <- as.numeric(names(mb_list)) - 1
    mb_num_tests <- result$num_tests

    # Create adjacency matrix based on Markov Blankets (mb-estimation.R)
    true_dag <- get_est_initial_dag(mb_list, p, verbose)
    # We are using a DAG with estimated Markov Blankets encoded
    est_dag <- TRUE
  } else {
    mb_list <- list()
    nodes_interest <- seq(0, p - 1)
    est_dag <- FALSE
  }

  # Convert any data.frame to a matrix
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
    results <- popCML_mag(
      true_dag, cpp_targets, nodes_interest,
      node_names, lmax, verbose
    )
  } else {
    results <- sampleCML_mag(
      true_dag, data, cpp_targets, nodes_interest,
      node_names, lmax, tol, verbose, test, est_dag
    )
  }

  # We change the target to target - 1 in order to accommodate change to C++
  list(
    "amat" = results$G,
    "S" = results$S,
    "NumTests" = results$NumTests,
    "MBNumTests" = mb_num_tests,
    "RulesUsed" = results$RulesUsed,
    "Nodes" = results$allNodes + 1, # to convert to R numbering
    "mbEstTime" = mb_time_track,
    "totalSkeletonTime" = results$totalSkeletonTime,
    "targetSkeletonTimes" = results$targetSkeletonTimes,
    "totalTime" = results$totalTime,
    "referenceDAG" = true_dag,
    "mbList" = mb_list,
    "data_means" = data_means,
    "data_cov" = data_cov
  )
}

# Warning and Message Functions -------------------------------------------

#' Print a summary of the Markov blanket estimation settings
#'
#' @param method Character; the MB estimation algorithm name.
#' @param test Character; the conditional independence test used.
#' @param threshold Numeric; the significance threshold used.
#' @return Invisibly, `NULL`; called for its side effect of printing to the
#'   console.
#' @noRd
mb_est_message <- function(method, test, threshold) {
  cat(
    "Estimating Markov Blankets using\n",
    "Algorithm:", method, "\n",
    "Test:", test, "\n",
    "Tolerance:", threshold, "\n"
  )
}

#' Validate a Markov blanket estimation significance threshold
#'
#' @param threshold Numeric; must lie in (0, 1].
#' @return Invisibly, `NULL`. Raises an error if `threshold` is out of range.
#' @noRd
validate_threshold <- function(threshold) {
  if (threshold <= 0 || threshold > 1) {
    stop("MB Estimation threshold is invalid. Threshold must be in (0,1]")
  }
}

#' Validate a Markov blanket estimation algorithm name
#'
#' @param method Character; must be one of `"MMPC"`, `"SES"`, `"gOMP"`,
#'   `"pc.sel"`, or `"MMMB"`.
#' @return Invisibly, `NULL`. Raises an error if `method` is not supported.
#' @noRd
validate_method <- function(method) {
  if (!(method %in% c("MMPC", "SES", "gOMP", "pc.sel", "MMMB"))) {
    stop("Invalid MB estimation algorithm")
  }
}

#' Validate a target node index
#'
#' @param target Integer; the target node index.
#' @param p Integer; the total number of nodes.
#' @return Invisibly, `NULL`. Raises an error if `target` is not in `1:p`.
#' @noRd
validate_target <- function(target, p) {
  if (!(target %in% seq(p))) {
    stop(paste0("Invalid target index (t=", target, ")"))
  }
}

# d-separation conversion function ----------------------------------------

#' Population-level (oracle) d-separation query for a target DAG
#'
#' @param true_dag A 0/1 adjacency matrix of the true DAG.
#' @param x,y 0-based indices of the two nodes to test for d-separation.
#' @param z A matrix whose first column gives the 0-based indices of the
#'   conditioning set; may have zero rows for an empty conditioning set.
#' @return Numeric `0`/`1`: `1` if `x` and `y` are d-separated by `z` in
#'   `true_dag`, `0` otherwise.
#' @details Builds a `bnlearn` graph object from `true_dag` and delegates to
#'   `bnlearn::dsep()`. Called from C++ via `condIndTestPop()` (see
#'   `src/pCorTest.cpp`) as the population-level counterpart to the
#'   sample-based partial correlation test.
#' @noRd
my_dsep <- function(true_dag, x, y, z) {
  tmp <- bnlearn::empty.graph(nodes = as.character(seq_len(ncol(true_dag))))
  if (!is.null(colnames(true_dag))) {
    colnames(true_dag) <- NULL
  }
  if (!is.null(rownames(true_dag))) {
    rownames(true_dag) <- NULL
  }
  bnlearn::amat(tmp) <- true_dag

  if (length(z) == 0) {
    res <- as.numeric(
      bnlearn::dsep(tmp, as.character(x + 1), as.character(y + 1))
    )
  } else {
    if (!is.null(dim(z))) {
      z <- z[, 1]
    }
    res <- as.numeric(
      bnlearn::dsep(
        tmp, as.character(x + 1), as.character(y + 1),
        as.character(z + 1)
      )
    )
  }
  res
}

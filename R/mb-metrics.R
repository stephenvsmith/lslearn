# Measurement Functions -------------------------------------------------

#' Identify a Node's Connections in an Adjacency Matrix
#'
#' @param g A 0/1 adjacency matrix.
#' @param target Integer index of the node whose connections to find.
#' @returns The indices of nodes connected to `target` in either direction
#'   (parent, child, or undirected neighbor).
#' @noRd
get_connections <- function(g, target) {
  # Get all target children or parents or undirected
  which(g[target, ] == 1 | g[, target] == 1)
}

#' Compare Parent Recovery Between a Reference and Estimated Graph
#'
#' @param ref,est 0/1 adjacency matrices for the reference (true) and
#'   estimated graphs.
#' @param target Integer index of the node whose parents to compare.
#' @returns A named numeric vector with `tp` (true positives), `fn` (false
#'   negatives), and `fp` (false positives) for `target`'s parents.
#' @noRd
calc_parent_recovery <- function(ref, est, target) {
  # Get all the parent nodes in both graphs
  parent_ref <- which(ref[, target] == 1 & ref[target, ] != 1)
  parent_est <- which(est[, target] == 1 & est[target, ] != 1)

  # If a node is in both vectors, then it is a true positive
  # If it is in the reference but not in the estimated, then it is a false
  # negative
  # If it is in the estimated but not the reference, then it is a false
  # positive
  c(
    "tp" = length(intersect(parent_ref, parent_est)),
    "fn" = length(setdiff(parent_ref, parent_est)),
    "fp" = length(setdiff(parent_est, parent_ref))
  )
}

#' Compare Child Recovery Between a Reference and Estimated Graph
#'
#' @inheritParams calc_parent_recovery
#' @returns A named numeric vector with `tp`, `fn`, and `fp` for `target`'s
#'   children.
#' @noRd
calc_child_recovery <- function(ref, est, target) {
  # Get all the child nodes in both graphs
  child_ref <- which(ref[target, ] == 1)
  child_est <- which(est[target, ] == 1)
  # If a node is in both vectors, then it is a true positive
  # If it is in the reference but not in the estimated, then it is a false
  # negative
  # If it is in the estimated but not the reference, then it is a false
  # positive
  c(
    "tp" = length(intersect(child_ref, child_est)),
    "fn" = length(setdiff(child_ref, child_est)),
    "fp" = length(setdiff(child_est, child_ref))
  )
}

#' Identify the Spouses of a Target Node
#'
#' A spouse of `target` is a parent of one of `target`'s children that is
#' not itself connected to `target`.
#'
#' @param g A 0/1 adjacency matrix.
#' @param children Indices of `target`'s children in `g`.
#' @param target Integer index of the target node.
#' @returns A sorted, unique vector of `target`'s spouses.
#' @noRd
get_spouses <- function(g, children, target) {
  sort(
    unique(
      unlist(
        lapply(children, function(child) {
          # Find parents of the child that are not connected to target
          child_parents <- which(
            g[, child] == 1 & g[child, ] == 0 &
              g[target, ] == 0 & g[, target] == 0
          )
          setdiff(child_parents, target)
        })
      )
    )
  )
}

#' Compare Spouse Recovery Between a Reference and Estimated Graph
#'
#' @inheritParams calc_parent_recovery
#' @returns A named numeric vector with `tp`, `fn`, and `fp` for `target`'s
#'   spouses.
#' @noRd
calc_spouse_recovery <- function(ref, est, target) {
  # Get all the child nodes in both graphs
  child_ref <- which(ref[target, ] == 1 & ref[, target] == 0)
  child_est <- which(est[target, ] == 1 & est[, target] == 0)

  spouses_ref <- get_spouses(ref, child_ref, target)
  spouses_est <- get_spouses(est, child_est, target)

  # If a node is in both vectors, then it is a true positive
  # If it is in the reference but not in the estimated, then it is a false
  # negative
  # If it is in the estimated but not the reference, then it is a false
  # positive
  c(
    "tp" = length(intersect(spouses_ref, spouses_est)),
    "fn" = length(setdiff(spouses_ref, spouses_est)),
    "fp" = length(setdiff(spouses_est, spouses_ref))
  )
}

#' Per-Target Markov Blanket Recovery Counts, Split by Role
#'
#' Input the initial matrix from the MB estimation algorithm and determine
#' which MB nodes were correctly identified, broken down into children,
#' parents, and spouses.
#'
#' @param ref A 0/1 adjacency matrix for the reference (true) graph.
#' @param est A 0/1 adjacency matrix for the estimated graph.
#' @param targets A vector of target node indices.
#' @returns A list of one-row data frames (one per target) with columns
#'   `mb_children_fn`, `mb_children_tp`, `mb_parents_fn`, `mb_parents_tp`,
#'   `mb_spouses_fn`, `mb_spouses_tp`, and `mb_total_fp`.
#' @noRd
mb_recovery_metrics_list <- function(ref, est, targets) {
  lapply(targets, function(target) {
    # Get estimated MB nodes
    mb_nodes <- get_connections(est, target)
    children <- which(ref[target, ] == 1 & ref[, target] == 0)
    parents <- which(ref[, target] == 1 & ref[target, ] == 0)
    spouses <- get_spouses(ref, children, target)
    data.frame(
      "mb_children_fn" = sum(!(children %in% mb_nodes)),
      "mb_children_tp" = sum(children %in% mb_nodes),
      "mb_parents_fn" = sum(!(parents %in% mb_nodes)),
      "mb_parents_tp" = sum(parents %in% mb_nodes),
      "mb_spouses_fn" = sum(!(spouses %in% mb_nodes)),
      "mb_spouses_tp" = sum(spouses %in% mb_nodes),
      "mb_total_fp" = sum(!(mb_nodes %in% c(children, parents, spouses)))
    )
  })
}

#' Aggregate Markov Blanket Recovery Counts Across Targets
#'
#' @inheritParams mb_recovery_metrics_list
#' @returns A one-row data frame with the same columns as
#'   `mb_recovery_metrics_list()`, summed across `targets` when there is
#'   more than one.
#' @noRd
mb_recovery_metrics <- function(ref, est, targets) {
  metrics_list <- mb_recovery_metrics_list(ref, est, targets)
  df <- as.data.frame(do.call(rbind, metrics_list))
  if (nrow(df) == 1) {
    df
  } else {
    df <- apply(df, 2, unlist)
    data.frame(t(colSums(df)))
  }
}

#' Identify the Spouses of a Target Node in a Reference Graph
#'
#' @param g A 0/1 adjacency matrix.
#' @param target Integer index of the target node.
#' @returns A sorted, unique vector of `target`'s spouses in `g`.
#' @noRd
spouse_recovery <- function(g, target) {
  # First, obtain children
  children <- which(g[target, ] == 1 & g[, target] == 0)
  # obtain reference spouses
  get_spouses(g, children, target)
}

#' Compare Full Markov Blanket Recovery for a Single Target
#'
#' Unlike `mb_recovery_metrics_list()`, this pools children, parents, and
#' spouses into a single Markov Blanket comparison rather than reporting
#' them separately.
#'
#' @inheritParams calc_parent_recovery
#' @returns A named numeric vector with `mb_tp`, `mb_fn`, and `mb_fp` for
#'   `target`'s full Markov Blanket.
#' @noRd
mb_recovery_target <- function(ref, est, target) {
  # Obtain all children and parents from Reference and Target Graphs
  ref_nodes <- get_connections(ref, target)
  est_nodes <- get_connections(est, target)

  # Add spouse nodes from reference graph
  ref_nodes <- c(ref_nodes, spouse_recovery(ref, target))

  # Compare MB recovery
  c(
    "mb_tp" = length(intersect(ref_nodes, est_nodes)),
    "mb_fn" = length(setdiff(ref_nodes, est_nodes)),
    "mb_fp" = length(setdiff(est_nodes, ref_nodes))
  )
}

#' Aggregate Full Markov Blanket Recovery Across Targets
#'
#' @param ref A 0/1 adjacency matrix for the reference (true) graph.
#' @param est A 0/1 adjacency matrix for the estimated graph.
#' @param targets A vector of target node indices.
#' @returns A one-row data frame with columns `mb_tp`, `mb_fn`, and `mb_fp`,
#'   summed across `targets` via `mb_recovery_target()`.
#' @noRd
mb_recovery <- function(ref, est, targets) {
  recovery_list <- lapply(targets, function(t) mb_recovery_target(ref, est, t))
  recovery_vec <- Reduce("+", recovery_list)
  data.frame(t(recovery_vec))
}

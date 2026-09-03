#' Compare the Outputs from the Estimated and True Plots
#'
#' Plots the estimated graph (restricted to the nodes actually considered by
#' the algorithm) side-by-side with the corresponding subgraph of the true
#' DAG, for visual comparison.
#'
#' @param local_output A list of the output from one of the local learning
#'   algorithms (e.g. `cml()` or `snl()`), with `amat` and `Nodes` entries.
#' @param true_dag A matrix containing the adjacency matrix for the true DAG.
#' @export
plot_output <- function(local_output, true_dag) {
  # Setup
  nodes_used <- local_output$Nodes
  node_names <- colnames(true_dag)[nodes_used]

  # Estimated graph
  g_est <- bnlearn::empty.graph(node_names)
  bnlearn::amat(g_est) <- local_output$amat[nodes_used, nodes_used]

  # True DAG
  g_true <- bnlearn::empty.graph(node_names)
  bnlearn::amat(g_true) <- true_dag[nodes_used, nodes_used]

  # Plot them side-by-side (estimated | true)
  graphics::par(mfrow = c(1, 2))
  bnlearn::graphviz.plot(g_est)
  bnlearn::graphviz.plot(g_true)
}

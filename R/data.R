#' Structure of the Asia Bayesian Network
#'
#' The 0/1 adjacency matrix of the well-known "Asia" (lung cancer) Bayesian
#' network (Lauritzen and Spiegelhalter, 1988), an 8-node example network
#' used as ground truth in lslearn's graph and Markov Blanket estimation
#' tests.
#'
#' @format An 8x8 numeric matrix with row and column names `asia`, `tub`,
#'   `smoke`, `lung`, `bronc`, `either`, `xray`, `dysp`. `asiaDAG[i, j] == 1`
#'   means node `i` is a parent of node `j`.
#' @source Ported from the CML package
#'   (\url{https://github.com/stephenvsmith/CML}).
"asiaDAG"

#' Simulated Continuous Data from the Asia Network Structure
#'
#' 500 observations of 8 continuous (Gaussian) variables simulated from the
#' structure in `asiaDAG`, used as test data for Markov Blanket estimation.
#'
#' @format A data frame with 500 rows and 8 numeric columns (`V1`-`V8`,
#'   corresponding in order to the nodes of `asiaDAG`).
#' @source Ported from the CML package
#'   (\url{https://github.com/stephenvsmith/CML}).
"asiadf"

#' Structure of the Andes Bayesian Network
#'
#' The 0/1 adjacency matrix of the "Andes" tutoring-system Bayesian network
#' (Conati et al., 1997), a 223-node, 338-edge example network used to
#' exercise graph-metric functions on a network much larger than `asiaDAG`.
#'
#' @format A 223x223 numeric matrix with row and column names given by the
#'   network's node names. `andesDAG[i, j] == 1` means node `i` is a parent
#'   of node `j`.
#' @source Ported from the CML package
#'   (\url{https://github.com/stephenvsmith/CML}).
"andesDAG"

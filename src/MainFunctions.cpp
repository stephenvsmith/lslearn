// Public C++ entry points backing lslearn's user-facing R wrappers
// (currently just snl(); cml()'s sampleCML()/popCML() will be added here
// once the CML class is ported).
#include "SNL.h"
#include <chrono>

using namespace std::chrono;
using namespace Rcpp;

//' Run the sample (or semi-sample) version of SNL and collect its results
//'
//' Internal building block for `snl()`; end users should call `snl()`
//' rather than this function directly.
//'
//' @param true_dag 0/1 adjacency matrix: the true DAG structure (semi-sample
//'   version) or an estimated Markov Blanket inclusion matrix (`estDAG =
//'   TRUE`).
//' @param df Data matrix (0-indexed internally).
//' @param targets,nodes_interest 0-based target and neighborhood-of-interest
//'   node indices.
//' @param names Node names.
//' @param lmax Maximum conditioning set size.
//' @param signif_level Significance level for conditional independence
//'   tests.
//' @param verbose Whether to print detailed output.
//' @param test Which conditional independence test to use.
//' @param estDAG Whether `true_dag` is an estimated Markov Blanket matrix
//'   rather than the true DAG.
//' @returns A list with the estimated graph (`G`), separating sets (`S`),
//'   test counts, node ordering, orientation rule usage, and timing
//'   information; see `snl()` for the user-facing result shape.
//' @noRd
// [[Rcpp::export]]
List sampleSNL(NumericMatrix true_dag, arma::mat df, NumericVector targets,
               NumericVector nodes_interest, StringVector names, int lmax = 3,
               double signif_level = 0.01, bool verbose = true,
               std::string test = "testIndFisher", bool estDAG = false) {
  // Variable to keep track of timing
  auto start = high_resolution_clock::now();

  // Instantiate the Local PC object
  SNL snl(true_dag, df, targets, nodes_interest, names, lmax, signif_level,
          verbose, test, estDAG);

  snl.run();

  auto end = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(end - start);
  double total_time = duration.count() / 1e6;
  total_time /= 60; // Get time in minutes

  return List::create(_["G"] = snl.getAmat(), _["S"] = snl.getSepSetList(),
                      _["NumTests"] = snl.getNumTests(),
                      _["allNodes"] = snl.getNeighborhood(),
                      _["rulesUsed"] = snl.getRulesUsed(),
                      _["targetSkeletonTimes"] = snl.getTargetSkeletonTimes(),
                      _["algorithmTotalTime"] = snl.getTotalTime(),
                      _["totalTime"] = total_time);
}

//' Run the population (oracle) version of SNL and collect its results
//'
//' Internal building block for `snl()`; end users should call `snl()`
//' rather than this function directly.
//'
//' @inheritParams sampleSNL
//' @returns Same shape as `sampleSNL()`, minus the sample-only
//'   `algorithmTotalTime` entry.
//' @noRd
// [[Rcpp::export]]
List popSNL(NumericMatrix true_dag, NumericVector targets,
            NumericVector nodes_interest, StringVector names, int lmax = 3,
            bool verbose = true) {
  // Variable to keep track of timing
  auto start = high_resolution_clock::now();

  // Instantiate the Local FCI object
  SNL snl(true_dag, targets, nodes_interest, names, lmax, verbose);

  if (verbose) {
    Rcout << "Beginning algorithm over all neighborhoods.\n";
  }

  snl.run();

  auto end = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(end - start);
  double total_time = duration.count() / 1e6;
  total_time /= 60; // Get time in minutes

  return List::create(_["G"] = snl.getAmat(), _["S"] = snl.getSepSetList(),
                      _["NumTests"] = snl.getNumTests(),
                      _["allNodes"] = snl.getNeighborhood(),
                      _["rulesUsed"] = snl.getRulesUsed(),
                      _["targetSkeletonTimes"] = snl.getTargetSkeletonTimes(),
                      _["totalTime"] = total_time);
}

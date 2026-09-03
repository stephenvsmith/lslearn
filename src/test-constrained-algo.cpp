// Test-only helpers exported to R purely so
// tests/testthat/test-constrained-algo.R can exercise ConstrainedAlgo from
// testthat. ConstrainedAlgo is abstract (getSkeletonTarget() is pure
// virtual, implemented by each concrete algorithm, e.g. SNL/CML, ported in
// a later phase), so TestAlgo below is a minimal concrete subclass that
// exists only so these tests can construct a ConstrainedAlgo and exercise
// its own methods (checkSeparation(), getVStructures(), accessors) in
// isolation, without running a full skeleton search. Not part of lslearn's
// public API or algorithm implementation.
#include "ConstrainedAlgo.h"

class TestAlgo : public ConstrainedAlgo {
public:
  using ConstrainedAlgo::ConstrainedAlgo;
  void getSkeletonTarget(const size_t &t) override {}
};

// [[Rcpp::export]]
List testConstructSample(NumericMatrix true_dag, arma::mat df,
                         NumericVector targets, NumericVector nodes_interest,
                         StringVector names, int lmax, double signif_level,
                         bool verbose, std::string test, bool estDAG) {
  TestAlgo algo(true_dag, df, targets, nodes_interest, names, lmax,
                signif_level, verbose, test, estDAG);
  return List::create(
      _["p"] = algo.getSize(), _["N"] = algo.getNeighborhood().size(),
      _["neighborhood"] = algo.getNeighborhood(), _["amat"] = algo.getAmat(),
      _["sepset"] = algo.getSepSetList(), _["num_tests"] = algo.getNumTests());
}

// [[Rcpp::export]]
List testConstructPop(NumericMatrix true_dag, NumericVector targets,
                      NumericVector nodes_interest, StringVector names,
                      int lmax, bool verbose) {
  TestAlgo algo(true_dag, targets, nodes_interest, names, lmax, verbose);
  return List::create(
      _["p"] = algo.getSize(), _["N"] = algo.getNeighborhood().size(),
      _["neighborhood"] = algo.getNeighborhood(), _["amat"] = algo.getAmat(),
      _["sepset"] = algo.getSepSetList());
}

// [[Rcpp::export]]
List testCheckSeparationPop(NumericMatrix true_dag, NumericVector targets,
                            NumericVector nodes_interest, StringVector names,
                            int lmax, bool verbose, int l, size_t i, size_t j,
                            NumericMatrix kvals) {
  TestAlgo algo(true_dag, targets, nodes_interest, names, lmax, verbose);
  algo.checkSeparation(l, i, j, kvals);
  return List::create(
      _["neighborhood"] = algo.getNeighborhood(), _["amat"] = algo.getAmat(),
      _["sepset"] = algo.getSepSetList(), _["adjacent"] = algo.isAdjacent(i, j),
      _["pval"] = algo.getMostRecentPVal(),
      _["num_tests"] = algo.getNumTests());
}

// [[Rcpp::export]]
List testCheckSeparationSample(NumericMatrix true_dag, arma::mat df,
                               NumericVector targets,
                               NumericVector nodes_interest, StringVector names,
                               int lmax, double signif_level, bool verbose,
                               std::string test, bool estDAG, int l, size_t i,
                               size_t j, NumericMatrix kvals) {
  TestAlgo algo(true_dag, df, targets, nodes_interest, names, lmax,
                signif_level, verbose, test, estDAG);
  algo.checkSeparation(l, i, j, kvals);
  return List::create(
      _["neighborhood"] = algo.getNeighborhood(), _["amat"] = algo.getAmat(),
      _["sepset"] = algo.getSepSetList(), _["adjacent"] = algo.isAdjacent(i, j),
      _["pval"] = algo.getMostRecentPVal(),
      _["num_tests"] = algo.getNumTests());
}

// [[Rcpp::export]]
List testGetVStructuresManual(NumericMatrix true_dag, NumericVector targets,
                              NumericVector nodes_interest, StringVector names,
                              int lmax, bool verbose, NumericMatrix custom_amat,
                              size_t sep_i, size_t sep_j, NumericVector sep_k) {
  TestAlgo algo(true_dag, targets, nodes_interest, names, lmax, verbose);
  // Directly configure the skeleton and separating set with the setters
  // CML's own header describes as "useful for testing", rather than
  // driving checkSeparation() over every pair: with C_tilde still complete
  // outside the one pair being tested, checkSeparation() alone leaves every
  // other node in the neighborhood as a spurious "common neighbor",
  // producing extra v-structures that don't reflect the true graph.
  algo.setAmat(custom_amat);
  algo.setS(sep_i, sep_j, sep_k);
  algo.setS(sep_j, sep_i, sep_k);
  int times_used = algo.getVStructures();
  return List::create(_["neighborhood"] = algo.getNeighborhood(),
                      _["amat"] = algo.getAmat(), _["times_used"] = times_used);
}

// [[Rcpp::export]]
void testPrintElements(NumericMatrix true_dag, NumericVector targets,
                       NumericVector nodes_interest, StringVector names,
                       int lmax, bool verbose) {
  TestAlgo algo(true_dag, targets, nodes_interest, names, lmax, verbose);
  algo.print_elements();
}

// [[Rcpp::export]]
List testAlgoSettersAndAccessors(NumericMatrix true_dag, NumericVector targets,
                                 NumericVector nodes_interest,
                                 StringVector names, int lmax, bool verbose,
                                 NumericMatrix new_amat,
                                 NumericVector new_neighbors) {
  TestAlgo algo(true_dag, targets, nodes_interest, names, lmax, verbose);
  bool verbose_before = algo.getVerbose();
  algo.setVerboseTrue();
  bool verbose_after = algo.getVerbose();
  algo.setAmat(new_amat);
  NumericMatrix amat_after_set = algo.getAmat();
  algo.setNeighbors(new_neighbors);
  return List::create(_["verbose_before"] = verbose_before,
                      _["verbose_after"] = verbose_after,
                      _["amat_after_set"] = amat_after_set,
                      _["neighborhood_after_set"] = algo.getNeighborhood());
}

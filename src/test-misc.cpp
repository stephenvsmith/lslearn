// Test-only helpers exported to R purely so tests/testthat/test-misc.R can
// exercise specific Rcpp/Armadillo container semantics (copy-vs-reference
// behavior, STL container interop, etc.) from testthat. These are not part
// of lslearn's public API or algorithm implementation.
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Rcpp's set union, sugar-vectorized union_().
// [[Rcpp::export]]
NumericVector test_union(NumericVector x, NumericVector y) {
  x = union_(x, y);
  return x;
}

// clone() + std::sort(), confirming clone() produces an independent copy.
// [[Rcpp::export]]
NumericVector test_sort(NumericVector x) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
}

// Constructing and filling a NumericMatrix via std::fill().
// [[Rcpp::export]]
NumericVector test_fill(int nrow, int ncol, int value) {
  NumericMatrix x(nrow, ncol);
  std::fill(x.begin(), x.end(), value);
  return x;
}

// NumericMatrix's fill_diag() convenience method.
// [[Rcpp::export]]
NumericMatrix test_fill_diag(NumericMatrix x, int value) {
  x.fill_diag(value);
  return x;
}

// NumericVector::create(), Rcpp's vector literal constructor.
// [[Rcpp::export]]
NumericVector test_create(double x1, double x2) {
  return NumericVector::create(x1, x2);
}

// Rcpp sugar's setdiff().
// [[Rcpp::export]]
NumericVector test_setdiff(NumericVector v1, NumericVector v2) {
  return setdiff(v1, v2);
}

// Rcpp sugar's intersect().
// [[Rcpp::export]]
NumericVector test_intersect(NumericVector v1, NumericVector v2) {
  return intersect(v1, v2);
}

// Building a std::map<int, int> from two parallel IntegerVectors.
// [[Rcpp::export]]
std::map<int, int> test_map_insert(IntegerVector v1, IntegerVector v2) {
  int n = v1.length();
  std::map<int, int> my_map;
  for (int i = 0; i < n; ++i) {
    my_map.insert(std::pair<int, int>(v1[i], v2[i]));
  }
  return my_map;
}

// std::map::find() on a map built by test_map_insert(); result is only
// checked via its printed output (see expect_output() in test-misc.R).
// [[Rcpp::export]]
void test_map_find(IntegerVector v1, IntegerVector v2, int a) {
  std::map<int, int> m = test_map_insert(v1, v2);
  std::map<int, int>::iterator it = m.find(a);
  Rcout << "(" << it->first << "," << it->second << ")\n";
}

// Confirms a default-constructed arma::uvec has size 0 (used elsewhere to
// represent an empty conditioning set).
// [[Rcpp::export]]
arma::uvec test_sep_arma() {
  arma::uvec sep_arma;
  Rcout << "\nSize of arma vector when l=0: " << sep_arma.size() << std::endl;
  return sep_arma;
}

// Column subsetting and arma::cor(), as used by the partial correlation
// tests in pCorTest.cpp.
// [[Rcpp::export]]
arma::mat test_subset_mat(arma::mat m, NumericVector i) {
  arma::uvec ind = as<arma::uvec>(i);
  return arma::cor(m.cols(ind));
}

// Confirms that clone() produces an independent copy while the original
// NumericMatrix argument is still passed by reference (so mutating `G` here
// is visible to the caller, but mutating `Gc` is not).
// [[Rcpp::export]]
NumericMatrix test_NumMat_value(NumericMatrix G) {
  NumericMatrix Gc = clone(G);
  Gc(0, 0) = 10;
  G(0, 0) = 25;
  return Gc;
}

// Confirms that a NumericMatrix passed by reference lets an Rcpp function
// mutate the caller's matrix in place.
// [[Rcpp::export]]
void test_decrement_matrix(NumericMatrix &G) {
  --G(0, 0);
  ++G(1, 1);
}

// Test-only helper exported to R purely so tests can exercise Armadillo's
// correlation function directly, independent of lslearn's own CI-testing
// code in pCorTest.cpp/pDCorTest.cpp.
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat testArmaCor(arma::mat M) { return arma::cor(M); }

#include <testthat.h>
#include <Rcpp.h>

// Forward declaration of the C++ function under test
double euclidean_dist(Rcpp::NumericVector x, Rcpp::NumericVector y);

context("Euclidean Distance C++ Unit Tests") {

  test_that("euclidean_dist handles standard 3D points correctly") {
    Rcpp::NumericVector p1 = Rcpp::NumericVector::create(0.0, 0.0, 0.0);
    Rcpp::NumericVector p2 = Rcpp::NumericVector::create(3.0, 4.0, 0.0);

    double dist = euclidean_dist(p1, p2);
    expect_true(std::abs(dist - 5.0) < 1e-7);
  }

  test_that("euclidean_dist returns 0 for identical vectors") {
    Rcpp::NumericVector p1 = Rcpp::NumericVector::create(1.5, 2.5);
    Rcpp::NumericVector p2 = Rcpp::NumericVector::create(1.5, 2.5);

    expect_true(euclidean_dist(p1, p2) == 0.0);
  }

  test_that("euclidean_dist throws an error for mismatched vector lengths") {
    Rcpp::NumericVector p1 = Rcpp::NumericVector::create(1.0, 2.0);
    Rcpp::NumericVector p2 = Rcpp::NumericVector::create(1.0, 2.0, 3.0);

    expect_error(euclidean_dist(p1, p2));
  }

}

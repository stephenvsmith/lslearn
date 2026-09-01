#include <Rcpp.h>
using namespace Rcpp;

//' Calculate Euclidean Distance between two vectors
//' @param x NumericVector
//' @param y NumericVector
//' @return double Euclidean distance
//' @export
// [[Rcpp::export]]
double euclidean_dist(NumericVector x, NumericVector y) {
  if (x.size() != y.size()) {
    stop("Vectors x and y must be of the same length.");
  }

  double sum_sq = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    double diff = x[i] - y[i];
    sum_sq += diff * diff;
  }

  return std::sqrt(sum_sq);
}

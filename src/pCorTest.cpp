// Conditional independence testing for continuous (Gaussian) data via
// partial correlation and Fisher's z-transformation.
#include "pCorTest.h"

//' Partial correlation between two variables given a conditioning set
//'
//' @param C Correlation matrix of all variables under consideration.
//' @param i,j 0-based indices into `C` of the two variables being tested.
//' @param k 0-based indices into `C` of the conditioning set; may be empty.
//' @return The partial correlation between variables `i` and `j` given the
//'   variables in `k`. For an empty conditioning set this is just `C(i, j)`;
//'   for one conditioning variable it uses the closed-form recursion; for
//'   larger conditioning sets it inverts the relevant submatrix of `C`.
//' @noRd
// [[Rcpp::export]]
double getPartialCorrelation(arma::mat C, size_t i, size_t j, arma::uvec k) {
  if (i == j) {
    return (1.0);
  }

  double pc;
  size_t k_size = k.size();

  if (k_size == 0) {
    return C(i, j);
  } else if (k_size == 1) {
    pc = (C(i, j) - C(i, k(0)) * C(j, k(0))) /
         sqrt((1 - pow(C(i, k(0)), 2)) * (1 - pow(C(j, k(0)), 2)));
  } else {
    arma::uvec indices(k_size + 2);
    indices(0) = i;
    indices(1) = j;
    for (size_t l = 0; l < k_size; ++l) {
      indices(l + 2) = k(l);
    }
    arma::mat Cinv = arma::pinv(C(indices, indices));
    pc = -Cinv(0, 1) / sqrt(Cinv(0, 0) * Cinv(1, 1));
  }

  return pc;
}

// log_e(1 + 2r / (1 - r)), the building block of Fisher's z-transformation.
// r == 1 is capped at a large finite value rather than returning Inf, since
// this only arises from a (near-)deterministic relationship in finite
// samples, not a real infinite test statistic.
double logPart(double r) {
  double result;
  if (r == 1) { // Perfect correlation means this value will approach infinity
    result = 1000;
  } else {
    result = log1p((2 * r) / (1 - r));
  }
  return result; // returns log_e(1 + (2*r/(1-r)))
}

//' Fisher z-transformed test statistic for a partial correlation
//'
//' @param pc A partial correlation, e.g. from `getPartialCorrelation()`.
//' @param n Sample size.
//' @param k_size Size of the conditioning set used to compute `pc`.
//' @return The test statistic `sqrt(n - k_size - 3) * atanh(pc)`, which is
//'   asymptotically standard normal under the null of zero partial
//'   correlation (Kalisch and Buhlmann, 2007).
//' @noRd
// [[Rcpp::export]]
double fisherZ(double pc, size_t n, size_t k_size) {
  return sqrt(n - k_size - 3) * 0.5 * logPart(pc);
}

//' Conditional independence test for continuous data via partial correlation
//'
//' @param C Correlation matrix of all variables under consideration.
//' @param i,j 0-based indices into `C` of the two variables being tested.
//' @param k 0-based indices into `C` of the conditioning set; may be empty.
//' @param n Sample size.
//' @param signif_level Significance level of the test.
//' @return A list with `result` (`TRUE` if the null of conditional
//'   independence is not rejected, i.e. `Xi` and `Xj` are declared
//'   conditionally independent given `Xk`), `statistic` (the Fisher
//'   z-transformed test statistic), and `pval` (the corresponding two-sided
//'   p-value).
//' @noRd
// [[Rcpp::export]]
List condIndTest(arma::mat &C, const size_t &i, const size_t &j,
                 const arma::uvec &k, const size_t &n,
                 const double &signif_level) {
  double pc = getPartialCorrelation(C, i, j, k);
  double statistic = fisherZ(pc, n, k.size());

  bool lower = statistic < 0;

  double cutoff = R::qnorm(1 - signif_level / 2, 0.0, 1.0, true, false);

  bool accept_H0 = std::abs(statistic) <= cutoff;
  double pval = 2 * R::pnorm(statistic, 0.0, 1.0, lower, false);
  // The null hypothesis is accepted (p-value large) => H_0: r = 0
  // => Conditional independence
  return List::create(_["result"] = accept_H0, _["statistic"] = statistic,
                      _["pval"] = pval);
}

//' Population-level (oracle) conditional independence test via d-separation
//'
//' @param G 0/1 adjacency matrix of the true DAG.
//' @param i,j 0-based node indices to test for d-separation.
//' @param k 0-based indices of the conditioning set; may be empty.
//' @return A list with `result` (`TRUE` if `i` and `j` are d-separated by
//'   `k` in `G`, i.e. conditionally independent given a CI oracle) and
//'   `pval` (1 if d-separated, 0 otherwise), for interface consistency with
//'   `condIndTest()`.
//' @details Delegates to the R-level `my_dsep()` helper (see
//'   `R/misc-helper.R`), which in turn uses `bnlearn::dsep()`. `my_dsep()`
//'   is internal (not exported), so it must be looked up through the
//'   package namespace rather than via a plain `Function("my_dsep")`
//'   constructor: the latter only finds it when the whole package is
//'   attached to the search path unnamespaced (e.g. under
//'   `pkgload::load_all()`), not for a normally `library()`-loaded
//'   package, where non-exported functions aren't on the search path by
//'   name.
//' @noRd
// [[Rcpp::export]]
List condIndTestPop(NumericMatrix G, const size_t &i, const size_t &j,
                    const arma::uvec &k) {
  Environment lslearn_ns = Environment::namespace_env("lslearn");
  Function my_dsep = lslearn_ns["my_dsep"];
  NumericVector tmp = my_dsep(G, i, j, k);
  double pval = tmp[0];
  bool accept_H0 = static_cast<bool>(pval);

  // The null hypothesis is accepted (p-value large) => H_0: r = 0
  // => Conditional independence
  return List::create(_["result"] = accept_H0, _["pval"] = pval);
}

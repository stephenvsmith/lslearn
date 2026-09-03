// Small utility functions shared across lslearn's structure-learning
// algorithms: enumerating combinations, membership testing, and (for
// debugging) printing vectors/matrices and generating placeholder node
// names.
#include "SharedFunctions.h"
using namespace Rcpp;

//' All combinations of size l from a vector
//'
//' @param x Numeric vector of values to choose from.
//' @param l Size of each combination.
//' @return A numeric matrix with `l` rows, one column per combination
//'   (delegating to R's `combn()` when there is more than one element to
//'   choose from). When `l` is 0, returns a 1x1 matrix containing `NA`,
//'   representing the (single) empty combination.
//' @details
//' Tested: 12/16/20. Confirmed: 11/23/22.
//'
//' `l` is taken as a signed `int` (rather than `size_t`) specifically so a
//' negative value can be rejected explicitly here: converting a negative
//' double straight to `size_t` at the Rcpp boundary is implementation- and
//' version-defined (it has been observed to both wrap around to a huge
//' positive value and to clamp to 0 across different Rcpp/compiler
//' versions), so it isn't a reliable way to reject bad input.
//' @noRd
// [[Rcpp::export]]
NumericMatrix combn_cpp(NumericVector x, int l) {
  if (l < 0) {
    stop("Combination size l must be non-negative");
  }
  size_t l_size = static_cast<size_t>(l);
  size_t xlength = x.length();
  NumericMatrix result;

  if (l_size > xlength) {
    stop("There aren't enough neighbors for the current value of l");
  }

  if (l_size == 0) {
    result = NumericMatrix(1, 1);
    result(0, 0) = NA_REAL;
  } else if (xlength == 1 && l_size == 1) {
    result = NumericMatrix(1, 1);
    result(0, 0) = x(0);
  } else if (l_size >= 1 && xlength > 1) {
    Function f("combn");
    result = f(Named("x") = x, _["m"] = l_size);
  }

  return result;
}

//' Test whether a value is present in a numeric vector
//'
//' @param x Numeric vector to search.
//' @param i Value to search for.
//' @return `TRUE` if `i` is an element of `x`, `FALSE` otherwise.
//' @details Tested: 11/23/22.
//' @noRd
// [[Rcpp::export]]
bool isMember(NumericVector x, const size_t &i) {
  NumericVector::iterator it = x.begin();
  size_t j;

  for (; it != x.end(); ++it) {
    j = *it;
    if (i == j) {
      return true;
    }
  }

  return false;
}

// PRINTING FUNCTIONS (DEBUGGING)
// Not exported to R; for ad hoc debugging from C++ only.

// Prints the elements of `v`, using `v` as indices into `names` (i.e. `v` is
// treated as a vector of 0-based indices, and the corresponding node names
// are printed rather than the raw values), separated by spaces and
// surrounded by `opening`/`closing`.
void printVecElements(NumericVector v, StringVector names, String opening,
                      String closing) {
  size_t ln = v.length();
  Rcout << opening.get_cstring();
  for (size_t i = 0; i < ln; ++i) {
    Rcout << names(v(i));
    if (i < ln - 1) {
      Rcout << " ";
    }
  }
  Rcout << closing.get_cstring();
}

// Prints the raw elements of `v` (unlike `printVecElements()`, no name
// lookup), separated by `sep` and surrounded by `opening`/`closing`.
void printVecElementsNoNames(NumericVector v, String opening, String closing,
                             String sep) {
  size_t l = v.length();
  Rcout << opening.get_cstring();
  for (size_t i = 0; i < l; ++i) {
    Rcout << v(i);
    if (i < l - 1) {
      Rcout << sep.get_cstring();
    }
  }
  Rcout << closing.get_cstring();
}

// Prints `m` row by row, space-separated within a row and newline-separated
// between rows.
void printMatrix(NumericMatrix m) {
  size_t n = m.nrow();
  size_t p = m.ncol();
  String ending;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < p; ++j) {
      if (j == p - 1) {
        ending = "\n";
      } else {
        ending = " ";
      }
      Rcout << m(i, j) << ending.get_cstring();
    }
  }
}

// Appends p placeholder node names ("V0", "V1", ..., "V(p-1)") to
// `node_names`, for use when the caller hasn't supplied real variable names.
void makeNodeNames(int p, StringVector &node_names) {
  for (int i = 0; i < p; ++i) {
    String node("V");
    node += i;
    node_names.push_back(node);
  }
}

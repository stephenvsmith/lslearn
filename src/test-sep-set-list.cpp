// Test-only helpers exported to R purely so
// tests/testthat/test-sep-set-list.R can exercise the SepSetList class from
// testthat (there is no public R-level interface to a C++ class, so each
// helper constructs a SepSetList, calls one method, and returns the
// result). Not part of lslearn's public API or algorithm implementation.
#include "SepSetList.h"
using namespace Rcpp;

// [[Rcpp::export]]
void printS(NumericVector neighbors) {
  SepSetList *S = new SepSetList(neighbors);
  S->printSepSetList();
  delete S;
}

// [[Rcpp::export]]
NumericVector getInitialValues(NumericVector neighbors, int i, int j) {
  SepSetList *S = new SepSetList(neighbors);
  NumericVector results = S->getSepSet(i, j);
  delete S;
  return results;
}

// [[Rcpp::export]]
NumericVector setListEmptySet(NumericVector neighbors, int i, int j) {
  SepSetList *S = new SepSetList(neighbors);
  S->changeList(i, j);
  S->printSepSetList();
  NumericVector results = S->getSepSet(i, j);
  delete S;
  return results;
}

// [[Rcpp::export]]
NumericVector setListEfficient(NumericVector neighbors, int i, int j,
                               NumericVector kvals) {
  SepSetList *S = new SepSetList(neighbors);
  S->changeList(i, j, kvals);
  S->printSepSetList();
  NumericVector results = S->getSepSet(i, j);
  delete S;
  return results;
}

// [[Rcpp::export]]
bool checkSeparationFunc(NumericVector neighbors, int i, int j,
                         NumericVector sep, int val_to_check) {
  SepSetList *S = new SepSetList(neighbors);
  S->changeList(i, j, sep);
  bool results = S->isSepSetMember(i, j, val_to_check);
  delete S;
  return results;
}

// [[Rcpp::export]]
bool checkSeparationFuncCorrected(NumericVector neighbors, int i, int j,
                                  NumericVector sep, int val_to_check) {
  SepSetList *S = new SepSetList(neighbors);
  S->changeList(i, j, sep);
  S->changeList(j, i, sep); // correction from previous function
  bool results = S->isSepSetMember(i, j, val_to_check);
  delete S;
  return results;
}

// [[Rcpp::export]]
bool checkIsSepSetMember(NumericVector neighbors, size_t i, size_t j,
                         NumericVector sep1, NumericVector sep2,
                         size_t val_to_check) {
  SepSetList *S = new SepSetList(neighbors);
  S->changeList(i, j, sep1);
  S->changeList(j, i, sep2);
  bool results = S->isSepSetMember(i, j, val_to_check);
  delete S;
  return results;
}

// [[Rcpp::export]]
bool checkPotentialVStruct(NumericVector neighbors, size_t i, size_t j,
                           NumericVector sep, size_t k) {
  SepSetList *S = new SepSetList(neighbors);
  S->changeList(i, j, sep);
  S->changeList(j, i, sep);
  bool results = S->isPotentialVStruct(i, j, k);
  delete S;
  return results;
}

// [[Rcpp::export]]
List checkGetS(NumericVector neighbors, size_t i, size_t j, NumericVector sep) {
  SepSetList *S = new SepSetList(neighbors);
  S->changeList(i, j, sep);
  S->changeList(j, i, sep);
  List L = S->getS();
  delete S;
  return L;
}

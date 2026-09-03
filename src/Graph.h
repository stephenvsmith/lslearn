// A general graph over a fixed set of nodes, backed by a numeric adjacency
// matrix; see Graph.cpp for documentation of each function. This is the
// base class for DAG and PDAG.
#ifndef GRAPH_H
#define GRAPH_H

#include "SharedFunctions.h"

// Base class
class Graph {
public:
  // Constructors
  Graph(size_t nodes, bool verbose = false); // tested
  Graph(size_t nodes, StringVector node_names, NumericMatrix adj,
        bool verbose = false); // tested

  void validateIndex(const size_t &ind) const { // tested
    if (ind >= m_p) {
      stop("Invalid index: too large");
    }
  }

  // Validates that a double (e.g. an element of an R-facing NumericVector of
  // targets) is a finite, non-negative whole number and in range, then
  // returns it cast to size_t. Prevents silent truncation of fractional
  // values and mishandling of NA/NaN when indices arrive as doubles.
  size_t validateAndCast(double val) const {
    if (!R_finite(val) || val != std::floor(val) || val < 0) {
      stop("Invalid index: must be a non-negative whole number");
    }
    size_t ind = static_cast<size_t>(val);
    validateIndex(ind);
    return ind;
  }

  void validateAdjMatrix(const NumericMatrix adj) { // tested
    if (adj.nrow() != adj.ncol()) {
      stop("Dimensions of adjacency matrix do not match.");
    }
  }

  // Accessors
  size_t size() { return m_p; }                      // not tested
  StringVector getNodeNames() { return m_names; }    // tested
  NumericMatrix getAmat() { return amat; }           // tested
  int getAmatVal(const size_t &i, const size_t &j) { // tested
    validateIndex(i);
    validateIndex(j);
    return amat(i, j);
  }
  int operator()(const size_t &i, const size_t &j) const { // tested
    validateIndex(i);
    validateIndex(j);
    return amat(i, j);
  }

  NumericVector getAmatRow(const size_t &row) { // tested
    validateIndex(row);
    return amat(row, _);
  }

  NumericVector getAmatCol(const size_t &col) { // tested
    validateIndex(col);
    return amat(_, col);
  }

  int getNCol() { return amat.ncol(); }          // not tested
  int getNRow() { return amat.nrow(); }          // not tested
  NumericVector getAdjacent(const size_t &i);    // tested
  NumericVector getNonAdjacent(const size_t &i); // tested

  // Setters
  void setSize(const size_t &s) { m_p = s; }
  void setNames(StringVector n) { m_names = n; }
  void setAmat(NumericMatrix m) { // tested
    validateAdjMatrix(m);
    amat = clone(m);
    m_p = m.ncol();
  }
  double &operator()(const size_t &i, const size_t &j) { // tested
    validateIndex(i);
    validateIndex(j);
    return amat(i, j);
  }
  void setAmatVal(const size_t &i, const size_t &j,
                  const size_t &val) { // tested
    validateIndex(i);
    validateIndex(j);
    amat(i, j) = val;
  }

  void setVerboseTrue() { verbose = true; }

  // Fills the adjacency matrix with zeros
  void emptyGraph(); // tested
  void printAmat();  // tested
  // Determines if two nodes are adjacent
  bool areAdjacent(const size_t &i, const size_t &j) { // tested
    validateIndex(i);
    validateIndex(j);
    return amat(i, j) != 0 || amat(j, i) != 0;
  }

  // check if there is a directed edge from i to j (i -> j)
  bool isDirected(size_t i, size_t j) {
    validateIndex(i);
    validateIndex(j);
    return amat(i, j) == 1 && amat(j, i) == 0;
  }

  // check if there is an undirected edge from i to j (i - j)
  bool isUndirected(size_t i, size_t j) {
    validateIndex(i);
    validateIndex(j);
    return amat(i, j) == 1 && amat(j, i) == 1;
  }

  // Functions to identify graph features
  // 1. Minimum Discriminating Path
  // 2. Minimum Uncovered Partially Directed Path
  NumericVector get_d_vals(size_t a, const std::vector<bool> &visited);
  NumericVector get_r_vals(size_t d, const std::vector<bool> &visited);
  NumericVector minDiscPath(size_t a, size_t b, size_t c); // tested
  bool isPathUncovered(NumericVector p);
  void checkWronglyCovered(NumericVector p);
  bool areEdgesPotentiallyDirected(size_t alpha, size_t beta);
  NumericVector idThetaVals(size_t alpha, size_t beta,
                            const std::vector<bool> &visited);
  NumericVector idUncovPdPath(size_t alpha, size_t beta, size_t gamma, size_t d,
                              NumericVector mpath);
  NumericVector uncovPdPath(size_t alpha, size_t beta, size_t gamma);
  NumericVector minUncovPdPath(size_t alpha, size_t beta,
                               size_t gamma); // tested

protected:
  NumericMatrix amat;
  bool verbose;

private:
  size_t m_p;
  StringVector m_names;
};

#endif

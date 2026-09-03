// CML ("Coordinated Multi-Neighborhood Learning"), the paper's main
// algorithm: a per-target local FCI-style skeleton search over the union of
// all target neighborhoods, followed by v-structure identification and the
// FCI orientation rules (1-4, 8-10), producing a PAG-style mixed graph that
// is then converted into a final adjacency matrix; see CML.cpp for
// documentation of each function.
//
// NOTE: Code developed for the FCI rules adapted from the implementation
// provided in the source code of the pcalg package:
// https://rdrr.io/cran/pcalg/src/R/udag2pag.R
#ifndef CML_H
#define CML_H

#include "ConstrainedAlgo.h"

using namespace std::chrono;

class CML : public ConstrainedAlgo {
public:
  CML(NumericMatrix true_dag, arma::mat df, NumericVector targets,
      NumericVector nodes_interest, StringVector names, int lmax,
      double signif_level, bool verbose, std::string test = "testIndFisher",
      bool estDAG = false); // tested

  CML(NumericMatrix true_dag, // population version
      NumericVector targets, NumericVector nodes_interest, StringVector names,
      int lmax, bool verbose); // tested

  void getSkeletonTotal();                 // tested
  void getSkeletonTarget(const size_t &t); // tested
  int getVStructures();                    // tested

  // Accessors
  std::vector<double> getTargetSkeletonTimes() {
    return target_skeleton_times;
  };
  double getTotalSkeletonTime() { return total_skeleton_time; };
  NumericVector getRulesCount() { return rules_used; };

  // Orientation Rules (tested)
  // Rule 1
  void rule1search(size_t beta, size_t alpha, bool &track_changes);
  bool rule1(bool &track_changes);
  // Rule 2
  void rule2search(size_t beta, size_t alpha, bool condition1, bool condition2,
                   bool &track_changes);
  bool rule2(bool &track_changes);
  // Rule 3
  List rule3asearch(size_t beta, size_t alpha);
  void rule3bsearch(const size_t &alpha, const size_t &beta,
                    const size_t &gamma, bool &track_changes);
  bool rule3(bool &track_changes);
  // Rule 4
  bool check_sep_r4(size_t beta, NumericVector md_path);
  bool rule4(bool &track_changes);
  // Rule 8
  bool rule8(bool &track_changes);
  // Rule 9
  bool rule9(bool &track_changes);
  // Rule 10
  bool rule10simple(const size_t &alpha, const size_t &beta,
                    const size_t &gamma, const size_t &d);
  bool rule10(bool &track_changes);
  // All rules
  void allRules();

  // Function to run the algorithm
  void run();     // tested in wrapper test
  void run_mag(); // tested in wrapper test

  // Graph conversion
  void convertMixedGraph(); // tested
  void convertFinalGraph(); // tested
  // Ensures we are using proper notation for each pair of nodes
  void checkNotation(); // tested

private:
  std::map<int, int> node_numbering;
  NumericVector rules_used = NumericVector(11);
  std::vector<double> target_skeleton_times;
  double total_skeleton_time;
};

#endif

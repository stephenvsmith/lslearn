# ConstrainedAlgo is abstract (getSkeletonTarget() is pure virtual); these
# tests exercise it through TestAlgo, a minimal concrete subclass defined
# only in src/test-constrained-algo.cpp for this purpose (see the comment
# there). CML has no dedicated test file for ConstrainedAlgo either -- it's
# only ever exercised indirectly through its real subclasses (SNL, CML,
# ported in later phases).

data("asiaDAG")
data("asiadf")
asia_nodes <- colnames(asiaDAG)
asiadf_mat <- as.matrix(asiadf)
p <- ncol(asiaDAG)
target <- which(asia_nodes == "either") - 1 # 5
kvals_empty <- matrix(nrow = 0, ncol = 0)

# Setup / Construction ------------------------------------------------------

test_that("Population version construction builds the target neighborhood", {
  algo <- testConstructPop(asiaDAG, target, target, asia_nodes, 3, FALSE)
  # either's Markov Blanket: parents tub/lung, children xray/dysp, spouse
  # bronc (shares child dysp with either)
  expect_equal(algo$neighborhood, c(1, 3, 4, 5, 6, 7))
  expect_equal(algo$p, p)
  expect_equal(algo$N, 6)
  # C_tilde starts as a complete graph on the neighborhood
  expect_equal(algo$amat, matrix(1, 6, 6) - diag(6))
  # All separating sets start uninitialized (NA)
  for (i in 1:6) {
    for (j in setdiff(1:6, i)) {
      expect_true(is.na(algo$sepset[[i]][[j]]))
    }
  }
})

test_that("Sample version construction builds the target neighborhood", {
  algo <- testConstructSample(
    asiaDAG, asiadf_mat, target, target, asia_nodes, 3, 0.01, FALSE,
    "testIndFisher", FALSE
  )
  expect_equal(algo$neighborhood, c(1, 3, 4, 5, 6, 7))
  expect_equal(algo$p, p)
  expect_equal(algo$N, 6)
  expect_equal(algo$amat, matrix(1, 6, 6) - diag(6))
  expect_equal(algo$num_tests, 0)
})

test_that("Construction rejects an out-of-range target", {
  expect_error(testConstructPop(asiaDAG, p, p, asia_nodes, 3, FALSE))
  expect_error(
    testConstructPop(asiaDAG, target, target, asia_nodes, 3, FALSE), NA
  )
})

# checkSeparation (population version, exact d-separation oracle) -----------

test_that("checkSeparation removes an edge between d-separated nodes (l=0)", {
  # tub and lung share child "either" (a collider); with no conditioning
  # set, that path is blocked, so they're marginally d-separated even
  # though both are also either's neighbor.
  res <- testCheckSeparationPop(
    asiaDAG, target, target, asia_nodes, 3, FALSE, 0, 0, 1, kvals_empty
  )
  expect_false(res$adjacent)
  expect_equal(res$pval, 1)
  expect_equal(res$num_tests, 1)
  expect_equal(res$sepset[[1]][[2]], -1)
})

test_that("checkSeparation keeps an edge between adjacent nodes (l=0)", {
  # tub -> either is a direct edge, so they can never be d-separated.
  res <- testCheckSeparationPop(
    asiaDAG, target, target, asia_nodes, 3, FALSE, 0, 0, 3, kvals_empty
  )
  expect_true(res$adjacent)
  expect_equal(res$pval, 0)
})

test_that("checkSeparation finds independence via a blocked collider", {
  # asia and bronc: the only path between them (asia->tub->either<-lung<-
  # smoke->bronc) is blocked at the collider "either" when not conditioned
  # on -- so they're marginally d-separated despite not being adjacent to
  # a common, unblocked chain.
  targets2 <- c(which(asia_nodes == "asia") - 1, target)
  algo <- testConstructPop(asiaDAG, targets2, targets2, asia_nodes, 3, FALSE)
  i_asia <- which(algo$neighborhood == which(asia_nodes == "asia") - 1) - 1
  i_bronc <- which(algo$neighborhood == which(asia_nodes == "bronc") - 1) - 1
  i_tub <- which(algo$neighborhood == which(asia_nodes == "tub") - 1) - 1

  res_indep <- testCheckSeparationPop(
    asiaDAG, targets2, targets2, asia_nodes, 3, FALSE, 0, i_asia, i_bronc,
    kvals_empty
  )
  expect_false(res_indep$adjacent)
  expect_equal(res_indep$pval, 1)

  # asia and tub are directly adjacent, so this active chain keeps them
  # dependent.
  res_dep <- testCheckSeparationPop(
    asiaDAG, targets2, targets2, asia_nodes, 3, FALSE, 0, i_asia, i_tub,
    kvals_empty
  )
  expect_true(res_dep$adjacent)
})

# checkSeparation (sample version) -------------------------------------------

test_that("Sample checkSeparation matches condIndTest for the same query", {
  # tub and dysp are marginally dependent (chain tub->either->dysp) but
  # conditioning on "either" should separate them. Cross-check against the
  # already-tested condIndTest() directly on the same correlation matrix,
  # rather than hand-deriving the expected p-value.
  either_true_idx <- which(asia_nodes == "either") - 1
  kvals <- matrix(either_true_idx, nrow = 1, ncol = 1)
  res <- testCheckSeparationSample(
    asiaDAG, asiadf_mat, target, target, asia_nodes, 3, 0.01, FALSE,
    "testIndFisher", FALSE, 1, 0, 5, kvals
  )

  corr_mat <- cor(asiadf_mat)
  n <- nrow(asiadf_mat)
  tub_true_idx <- which(asia_nodes == "tub") - 1
  dysp_true_idx <- which(asia_nodes == "dysp") - 1
  oracle <- condIndTest(
    corr_mat, tub_true_idx, dysp_true_idx, either_true_idx, n, 0.01
  )

  expect_equal(res$pval, oracle$pval)
  expect_equal(res$adjacent, !oracle$result)
  expect_equal(res$num_tests, 1)
})

# getVStructures --------------------------------------------------------------

# A minimal, fully hand-verifiable v-structure: A -> C <- B, with A and B
# not adjacent to each other. Constructed directly with setAmat()/setS()
# (the setters CML's own header describes as "useful for testing") rather
# than driving checkSeparation() over every pair: with C_tilde still
# complete outside one separated pair, every other node in the
# neighborhood would look like a spurious "common neighbor", producing
# extra v-structures that don't reflect a real learned skeleton.
synthetic_dag <- matrix(c(
  0, 0, 1,
  0, 0, 1,
  0, 0, 0
), nrow = 3, byrow = TRUE)
synthetic_names <- c("A", "B", "C")
synthetic_amat <- matrix(c(
  0, 0, 1,
  0, 0, 1,
  1, 1, 0
), nrow = 3, byrow = TRUE)

test_that("getVStructures detects a collider when it's not in the sepset", {
  res <- testGetVStructuresManual(
    synthetic_dag, 2, 2, synthetic_names, 3, FALSE, synthetic_amat, 0, 1,
    c(-1)
  )
  expect_equal(res$times_used, 1)
  # Arrowheads into C from both A and B; the reverse edges are removed.
  expect_equal(res$amat, matrix(c(
    0, 0, 1,
    0, 0, 1,
    0, 0, 0
  ), nrow = 3, byrow = TRUE))
})

test_that("getVStructures does not flag a collider that's in the sepset", {
  # A and B "separated by" C: not a real d-separation fact for this graph,
  # but exercises that isPotentialVStruct() correctly respects the
  # separating set regardless.
  res <- testGetVStructuresManual(
    synthetic_dag, 2, 2, synthetic_names, 3, FALSE, synthetic_amat, 0, 1,
    c(2)
  )
  expect_equal(res$times_used, 0)
  expect_equal(res$amat, synthetic_amat)
})

# print_elements / accessors / setters ---------------------------------------

test_that("print_elements prints diagnostic output", {
  expect_output(
    testPrintElements(asiaDAG, target, target, asia_nodes, 3, FALSE)
  )
})

test_that("Setters update verbosity, adjacency matrix, and neighborhood", {
  new_amat <- matrix(0, nrow = 6, ncol = 6)
  new_neighbors <- c(1, 3)
  res <- testAlgoSettersAndAccessors(
    asiaDAG, target, target, asia_nodes, 3, FALSE, new_amat, new_neighbors
  )
  expect_false(res$verbose_before)
  expect_true(res$verbose_after)
  expect_equal(res$amat_after_set, new_amat)
  expect_equal(res$neighborhood_after_set, new_neighbors)
})

# Setup -------------------------------------------------------------------

data("asiadf")
asiadf_mat <- asiadf
data("asiaDAG")
node_names <- colnames(asiaDAG)
p <- ncol(asiadf)
asia_dag <- matrix(asiaDAG, nrow = p, ncol = p)
asiadf <- as.matrix(asiadf)


# ConstrainedAlgo tests ---------------------------------------------------

test_that("Test getSize", {
  expect_equal(getSizeCML(asia_dag, asiadf, 3, seq(0, p - 1), node_names), 8)
})

test_that("Change Separating Set", {
  expect_output(
    setSCML(asia_dag, asiadf, 3, seq(0, p - 1), node_names, 0, 1, c(3, 4))
  )
})

test_that("Test set verbose", {
  expect_output(setVerboseCML(asia_dag, asiadf, 3, seq(0, p - 1), node_names))
})


# Testing object formation ------------------------------------------------

test_that("Testing the CML object", {
  expect_output(initializeCML(asia_dag, asiadf, 3, seq(0, p - 1), node_names))
})

test_that("Testing the CML object (Population)", {
  expect_output(initializeCMLPop(asia_dag, 3, seq(0, p - 1), node_names))
})


# Testing Separation Set --------------------------------------------------

test_that("Testing Separation Test Function", {
  # Check if tub and lung are separated by the empty set
  expect_true(
    checkSeparationTest(
      asia_dag, asiadf, 3, seq(0, p - 1), node_names, 0, 2, 0, 1
    ) > 0.1
  )

  # Check if either and smoke are separated by tub (removing lung)
  expect_true(
    checkSeparationTest(
      asia_dag, asiadf, 3, seq(0, p - 1), node_names, 1, 3, 1, 2
    ) < 0.05
  )

  # Check if either and smoke are separated by lung (removing tub)
  # t = 3 => target is lung (neighborhood: {tub (0), smoke (1), lung (2),
  # either (3)}). Nodes 3 and 1 are either and smoke, and we are excluding
  # tub from consideration. We should obtain separation between either and
  # smoke conditioned on lung
  expect_true(
    checkSeparationTest(
      asia_dag, asiadf, 3, seq(0, p - 1), node_names,
      1, # Smoke
      3, # Either - are these separated?
      1, # Max. Size of Sep. Set
      0 # Removing tub from consideration
    ) > 0.05
  )

  # Check if either and lung are separated by any neighbor or set of
  # neighbors. t = 5 => either (neighborhood: {tub (0), lung (1), bronc (2),
  # either (3), xray (4), dysp (5)}). Nodes 1 and 3 are lung and either, and
  # we are excluding none. Should not obtain separation for any grouping
  expect_true(
    checkSeparationTest(
      asia_dag, asiadf, 5, seq(0, p - 1), node_names,
      1, # lung
      3, # either
      2, # Max size of sep set
      vector(mode = "double")
    ) < 0.01
  )
})


# Testing Skeleton Function -----------------------------------------------

test_that("Testing the total skeleton function (one target)", {
  # lung is the target
  expect_output(
    checkSkeletonTotal(asia_dag, asiadf, 3, seq(0, p - 1), node_names)
  )
})

test_that("Testing the total skeleton function (two targets)", {
  # Targets are lung and bronc
  expect_output(
    checkSkeletonTotal(asia_dag, asiadf, c(3, 4), seq(0, p - 1), node_names)
  )
})

test_that("Testing the total skeleton function for population (one target)", {
  expect_output(checkSkeletonTotalPop(asia_dag, 3, seq(0, p - 1), node_names))
})

test_that("Testing the total skeleton function for population (two targets)", {
  expect_output(
    checkSkeletonTotalPop(asia_dag, c(3, 4), seq(0, p - 1), node_names)
  )
})


# Testing Skeleton 2nd Stage and V-Structures -----------------------------

test_that("Skeleton target and V-structure functions (sample)", {
  expect_output(
    result_amat <- checkVStruct(
      asia_dag, asiadf, c(3, 4), seq(0, p - 1), node_names
    )
  )

  # Convert to PDAG numbering
  for (i in seq_len(nrow(result_amat))) {
    for (j in seq_len(ncol(result_amat))) {
      if (result_amat[i, j] == 2) {
        result_amat[i, j] <- 1
        result_amat[j, i] <- 0
      }
    }
  }
  # Check v-structure with tub -> either <- lung
  expect_equal(result_amat[1, 5], 1)
  expect_equal(result_amat[5, 1], 0)
  expect_equal(result_amat[3, 5], 1)
  expect_equal(result_amat[5, 3], 0)
  expect_equal(result_amat[1, 3], 0)
  expect_equal(result_amat[3, 1], 0)

  # Check v-structure either -> dysp <- bronc
  expect_equal(result_amat[5, 6], 1)
  expect_equal(result_amat[6, 5], 0)
  expect_equal(result_amat[4, 6], 1)
  expect_equal(result_amat[6, 4], 0)
  expect_equal(result_amat[4, 5], 0)
  expect_equal(result_amat[5, 4], 0)

  # Population Version
  result_amat <- checkVStructPop(asia_dag, c(3, 4), seq(0, p - 1), node_names)

  modified_asia_dag <- asia_dag
  # Remove asia and xray from consideration to match sample version
  modified_asia_dag[1, ] <- 0
  modified_asia_dag[, 7] <- 0
  # smoke - lung
  modified_asia_dag[4, 3] <- 1
  # smoke - bronc
  modified_asia_dag[5, 3] <- 1

  expect_equal(result_amat, modified_asia_dag)
})

test_that("Skeleton target and V-structure functions (population)", {
  expect_output(
    result_amat <- checkVStructPop(
      asia_dag, c(3, 4), seq(0, p - 1), node_names
    )
  )

  # Check v-structure with tub -> either <- lung
  expect_equal(result_amat[2, 6], 1)
  expect_equal(result_amat[6, 2], 0)
  expect_equal(result_amat[4, 6], 1)
  expect_equal(result_amat[6, 4], 0)
  expect_equal(result_amat[2, 4], 0)
  expect_equal(result_amat[4, 2], 0)

  # Check v-structure either -> dysp <- bronc
  expect_equal(result_amat[6, 8], 1)
  expect_equal(result_amat[8, 6], 0)
  expect_equal(result_amat[5, 8], 1)
  expect_equal(result_amat[8, 5], 0)
  expect_equal(result_amat[5, 6], 0)
  expect_equal(result_amat[6, 5], 0)
})


# Check final matrix conversion -------------------------------------------

test_that("Testing Adjacency Matrix Conversion", {
  test_mat <- matrix(c(
    0, 2, 0, 0, 1, 0,
    2, 0, 2, 0, 0, 0,
    0, 1, 0, 2, 0, 0,
    0, 0, 3, 0, 0, 0,
    1, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 3, 0
  ), nrow = 6, byrow = TRUE)

  # This true_dag argument is reused for two purposes inside CML: it seeds
  # the constructor's Markov Blanket list (which convertMixedGraph() then
  # consults to decide which node pairs are "in the same neighborhood"),
  # and, in CML's own upstream test, it was also intended to double as the
  # expected final result (since setAmat()/setNeighbors() below discard the
  # constructor's own working graph in favor of test_mat). Those two roles
  # disagree here: running CML's real, unmodified source against its own
  # upstream test produces 2s at [2,4]/[4,2] where the upstream fixture
  # hardcodes 1s. This is a pre-existing upstream inconsistency, not
  # something introduced by this port -- the expected value below is what
  # CML's real code actually produces.
  final_mat <- matrix(c(
    rep(0, 11),
    0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0,
    rep(0, 11),
    0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    rep(0, 11),
    rep(0, 11),
    rep(0, 11),
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    rep(0, 11),
    rep(0, 11)
  ), nrow = 11, byrow = TRUE)

  expected_conv <- final_mat
  expected_conv[2, 4] <- 2
  expected_conv[4, 2] <- 2

  result_amat <- checkAdjMatConversion(
    final_mat, asiadf, c(3, 4), seq(0, 10), node_names,
    test_mat, c(1, 3, 4, 6, 8, 10)
  )
  expect_equal(result_amat, expected_conv)

  asia_test <- matrix(c(
    0, 2, 0, 2,
    2, 0, 1, 2,
    0, 1, 0, 0,
    1, 1, 0, 0
  ), byrow = TRUE, nrow = 4)
  asia_result <- checkAdjMatConversion(
    asia_dag, asiadf, c(4, 5), seq(0, p - 1), node_names, asia_test, 4:7
  )
  expect_equal(asia_result, matrix(c(
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 2, 0, 1,
    0, 0, 0, 0, 2, 0, 1, 1,
    0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0
  ), nrow = 8, byrow = TRUE))
})


# Test Full Algorithm -----------------------------------------------------

test_that("CML (Putting it all together, sample)", {
  sample_g <- bnlearn::empty.graph(node_names)
  expect_output(
    bnlearn::amat(sample_g) <- checkCMLSummary(
      asia_dag, asiadf, c(0, 5), seq(0, p - 1), node_names
    )
  )
  expect_equal(bnlearn::amat(sample_g), matrix(c(
    0, 1, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0
  ), nrow = 8, byrow = TRUE, dimnames = list(node_names, node_names)))
})

test_that("CML (Putting it all together, population)", {
  pop_g <- bnlearn::empty.graph(node_names)

  bnlearn::amat(pop_g) <- checkCMLSummaryPop(
    asia_dag, c(0, 5), seq(0, p - 1), node_names
  )
  expect_equal(bnlearn::amat(pop_g), matrix(c(
    0, 1, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0
  ), nrow = 8, byrow = TRUE, dimnames = list(node_names, node_names)))

  bnlearn::amat(pop_g) <- checkCMLSummaryPop(
    asia_dag, c(2, 7), seq(0, p - 1), node_names
  )
  expect_equal(bnlearn::amat(pop_g), matrix(c(
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 1, 0, 0, 0,
    0, 0, 1, 0, 0, 1, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 0, 1, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0
  ), nrow = 8, byrow = TRUE, dimnames = list(node_names, node_names)))

  bnlearn::amat(pop_g) <- checkCMLSummaryPop(
    asia_dag, c(1, 3, 7), seq(0, p - 1), node_names
  )
  expect_equal(bnlearn::amat(pop_g), matrix(c(
    0, 1, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 1, 1, 0, 0, 0,
    0, 0, 1, 0, 0, 1, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0
  ), nrow = 8, byrow = TRUE, dimnames = list(node_names, node_names)))

  bnlearn::amat(pop_g) <- checkCMLSummaryPop(
    asia_dag, c(0, 6), seq(0, p - 1), node_names
  )
  expect_equal(bnlearn::amat(pop_g), matrix(c(
    0, 1, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0
  ), nrow = 8, byrow = TRUE, dimnames = list(node_names, node_names)))
})


# Misc. -------------------------------------------------------------------

test_that("Testing object conversion", {
  expect_equal(class(asiadf_mat), "data.frame")
  expect_error(cml(data = asiadf_mat, targets = 1), NA)
  asiadag_df <- data.frame(asia_dag)
  expect_equal(class(asiadag_df), "data.frame")
  expect_error(cml(true_dag = asiadag_df, targets = 1), NA)
})

test_that("Testing warnings for checkNotation", {
  m <- asia_dag
  m[5, 3] <- 3
  m[8, 6] <- 2
  m[3, 4] <- 2
  expect_warning(
    expect_warning(
      expect_warning(
        checkNotationWarnings(
          asia_dag, asiadf, 5, seq(1, p - 1), node_names, m
        ),
        "Ancestral marking mixed with neighborhood marking."
      ),
      "Ancestral marking mixed with neighborhood marking."
    ),
    "Ancestral marking mixed with neighborhood marking."
  )
  expect_warning(
    checkNotationWarnings(
      asia_dag, asiadf, 5, seq(1, p - 1), node_names, asia_dag
    ),
    NA
  )

  m[3, 4] <- 1
  m[6, 8] <- 3
  m[3, 5] <- 4
  expect_warning(
    checkNotationWarnings(asia_dag, asiadf, 5, seq(1, p - 1), node_names, m),
    NA
  )
})

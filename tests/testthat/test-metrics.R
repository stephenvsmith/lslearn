# Setup -------------------------------------------------------------------

data("asiaDAG")
data("asiadf")
data("andesDAG")
nodes <- colnames(asiaDAG)


# Check Supporting Functions ----------------------------------------------

test_that("Ensure the returned number of edges is correct", {
  # The asia graph has 8 edges/arcs
  expect_equal(getEdgeNumber(asiaDAG), 8)
  # The andes graph has 338 arcs
  expect_equal(getEdgeNumber(andesDAG), 338)
})

test_that("Check function to determine mutual neighborhoods", {
  # targets: tub, bronc | checking: either
  expect_true(inTargetNeighborhood(asiaDAG, c(2, 5) - 1, 5))
  # checking: dysp
  expect_true(inTargetNeighborhood(asiaDAG, c(2, 5) - 1, 7))
  # checking: lung
  expect_true(inTargetNeighborhood(asiaDAG, c(2, 5) - 1, 3))
  # checking: xray
  expect_false(inTargetNeighborhood(asiaDAG, c(2, 5) - 1, 6))
  # targets: asia, xray | checking: smoke
  expect_false(inTargetNeighborhood(asiaDAG, c(1, 7) - 1, 2))
  # target is lung, checking bronc
  expect_false(inTargetNeighborhood(asiaDAG, 3, 4))

  # targets: tub, bronc | checking: either, bronc
  expect_true(sharedNeighborhood(asiaDAG, c(2, 5) - 1, 5, 4))
  # checking: either, dysp
  expect_true(sharedNeighborhood(asiaDAG, c(2, 5) - 1, 5, 7))
  # targets: asia, dysp | checking: tub, either
  expect_false(sharedNeighborhood(asiaDAG, c(1, 8) - 1, 1, 5))
  # targets: asia, dysp | checking: tub, bronc
  expect_false(sharedNeighborhood(asiaDAG, c(1, 8) - 1, 1, 4))
  # targets: lung, either | checking: tub, either
  expect_true(sharedNeighborhood(asiaDAG, c(4, 6) - 1, 1, 5))
  # targets: either | checking: tub, either
  expect_true(sharedNeighborhood(asiaDAG, c(6) - 1, 1, 5))
  # targets: asia, either | checking: tub, dysp
  expect_true(sharedNeighborhood(asiaDAG, c(1, 6) - 1, 1, 7))
  # targets: asia, either | checking: asia, either
  expect_false(sharedNeighborhood(asiaDAG, c(1, 6) - 1, 0, 5))
  # targets: asia, either, tub | checking: asia, either
  expect_true(sharedNeighborhood(asiaDAG, c(1, 6, 2) - 1, 0, 5))

  test_amat <- matrix(c(
    0, 1, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 1, 0
  ), byrow = TRUE, ncol = 6)
  expect_false(sharedNeighborhood(test_amat, c(0, 2), 5, 4))
  expect_true(sharedNeighborhood(test_amat, c(0, 4), 5, 4))
})

test_that("Test inTargetNeighborhood with verbose", {
  expect_output(result <- inTargetNeighborhood(asiaDAG, c(2, 5) - 1, 5, TRUE))
  expect_true(result)
})

test_that("SharedNeighborhood with verbose (I)", {
  test_amat <- matrix(c(
    0, 1, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 1, 0
  ), byrow = TRUE, ncol = 6)
  expect_true(sharedNeighborhood(test_amat, c(0, 4), 5, 4))
  expect_output(result <- sharedNeighborhood(test_amat, c(0, 2), 5, 4, TRUE))
  expect_false(result)
})

test_that("SharedNeighborhood with verbose (II)", {
  test_amat <- matrix(c(
    0, 1, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 1, 0
  ), byrow = TRUE, ncol = 6)
  expect_output(result <- sharedNeighborhood(test_amat, c(0, 4), 5, 4, TRUE))
  expect_true(result)
})


# Check metric functions ----------------------------------------------------

### First set of checks

true_amat <- matrix(c(
  0, 1, 0, 0,
  1, 0, 1, 1,
  0, 0, 0, 1,
  0, 0, 0, 0
), byrow = TRUE, nrow = 4)

perfect_skel <- matrix(c(
  0, 1, 0, 0,
  0, 0, 0, 1,
  0, 1, 0, 1,
  0, 0, 0, 0
), byrow = TRUE, nrow = 4)

false_skel <- matrix(c(
  0, 0, 1, 1,
  0, 0, 1, 1,
  0, 0, 0, 0,
  0, 0, 0, 0
), byrow = TRUE, nrow = 4)

test_that("checking skeleton comparison function", {
  expect_equal(
    compareSkeletons(false_skel, true_amat),
    list("skel_fp" = 2, "skel_fn" = 2, "skel_correct" = 2)
  )
  expect_equal(
    compareSkeletons(perfect_skel, true_amat),
    list("skel_fp" = 0, "skel_fn" = 0, "skel_correct" = 4)
  )
})

test_that("checking v structure comparison functions", {
  expect_equal(
    compareVStructures(perfect_skel, true_amat, TRUE),
    list("missing" = 0, "added" = 1, "correct" = 0)
  )
  expect_equal(
    compareVStructures(false_skel, true_amat, TRUE),
    list("missing" = 0, "added" = 2, "correct" = 0)
  )
})

test_that("checking parent recovery accuracy functions", {
  expect_equal(
    parentRecoveryAccuracy(perfect_skel, true_amat, targets = 3),
    list("missing" = 0, "added" = 0, "correct" = 2, "potential" = 0)
  )
  expect_equal(
    parentRecoveryAccuracy(false_skel, true_amat, targets = 3),
    list("missing" = 1, "added" = 1, "correct" = 1, "potential" = 0)
  )
  expect_equal(
    parentRecoveryAccuracy(false_skel, true_amat, targets = c(0, 3)),
    list("missing" = 1, "added" = 1, "correct" = 1, "potential" = 0)
  )
  expect_equal(
    parentRecoveryAccuracy(false_skel, true_amat, targets = c(2, 3)),
    list("missing" = 1, "added" = 2, "correct" = 2, "potential" = 0)
  )
  # Adding a potential
  false_skel[3, 4] <- 1
  false_skel[4, 3] <- 1
  expect_equal(
    parentRecoveryAccuracy(false_skel, true_amat, targets = 3),
    list("missing" = 0, "added" = 1, "correct" = 1, "potential" = 1)
  )
  # Using 1-index numbering: Missing 3 -> 4, correctly has 2->3 and 2->4,
  # has 3-4 as undirected edge, added 1->4 and 1->3
  expect_equal(
    parentRecoveryAccuracy(false_skel, true_amat, targets = c(2, 3)),
    list("missing" = 0, "added" = 2, "correct" = 2, "potential" = 1)
  )
  expect_equal(
    parentRecoveryAccuracy(false_skel, true_amat, targets = c(1, 2, 3)),
    list("missing" = 0, "added" = 2, "correct" = 2, "potential" = 1)
  )
  false_skel[3, 4] <- 2
  false_skel[4, 3] <- 2
  expect_equal(
    parentRecoveryAccuracy(false_skel, true_amat, targets = 3),
    list("missing" = 1, "added" = 1, "correct" = 1, "potential" = 0)
  )
  expect_equal(
    parentRecoveryAccuracy(false_skel, true_amat, targets = c(2, 3)),
    list("missing" = 1, "added" = 2, "correct" = 2, "potential" = 0)
  )
  # This looks correct, but it is an ancestral edge
  false_skel[4, 3] <- 3
  expect_equal(
    parentRecoveryAccuracy(false_skel, true_amat, targets = 3),
    list("missing" = 1, "added" = 1, "correct" = 1, "potential" = 0)
  )
  expect_equal(
    parentRecoveryAccuracy(false_skel, true_amat, targets = c(2, 3)),
    list("missing" = 1, "added" = 2, "correct" = 2, "potential" = 0)
  )
  # This looks correct, but it is once again an ancestral edge
  false_skel[4, 3] <- 4
  expect_equal(
    parentRecoveryAccuracy(false_skel, true_amat, targets = 3),
    list("missing" = 1, "added" = 1, "correct" = 1, "potential" = 0)
  )
  expect_equal(
    parentRecoveryAccuracy(false_skel, true_amat, targets = c(2, 3)),
    list("missing" = 1, "added" = 2, "correct" = 2, "potential" = 0)
  )
})


# Advanced Skeleton Testing -------------------------------------------------

true_skeleton <- matrix(c(
  0, 1, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 1, 0, 1, 0, 0,
  0, 0, 0, 0, 0, 1, 1, 0, 0,
  0, 0, 0, 0, 0, 1, 0, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 1, 1, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 1, 0
), byrow = TRUE, nrow = 9, ncol = 9)

est_skeleton <- matrix(c(
  0, 4, 0, 0, 0, 0, 0, 0, 0,
  4, 0, 2, 0, 0, 0, 0, 0, 0,
  0, 2, 0, 0, 1, 0, 1, 0, 0,
  0, 0, 1, 0, 1, 2, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 0, 1, 1,
  0, 0, 0, 4, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 3,
  0, 0, 0, 0, 1, 0, 0, 2, 0
), byrow = TRUE, nrow = 9, ncol = 9)

test_that("Skeleton metrics with ancestral edges", {
  expect_equal(
    compareSkeletons(est_skeleton, true_skeleton, verbose = TRUE),
    list("skel_fp" = 2, "skel_fn" = 6, "skel_correct" = 4)
  )
})


# Advanced PRA Testing -------------------------------------------------------

true_amat <- matrix(0, ncol = 11, nrow = 11)
true_amat[3, 1] <- true_amat[3, 2] <- 1
true_amat[4, 3] <- true_amat[3, 4] <- 1
true_amat[5, 3] <- 1
true_amat[6, 5] <- 1
true_amat[7, 6] <- 1
true_amat[8, 11] <- 1
true_amat[9, 11] <- 1
true_amat[10, 7] <- 1
true_amat[11, 4] <- 1

test_amat <- matrix(0, ncol = 11, nrow = 11)
test_amat[1, 3] <- test_amat[3, 1] <- 1
test_amat[2, 3] <- test_amat[3, 2] <- 1
test_amat[3, 11] <- test_amat[11, 3] <- 1
test_amat[4, 3] <- 1
test_amat[5, 6] <- test_amat[6, 5] <- 1
test_amat[7, 3] <- 2
test_amat[3, 7] <- 4
test_amat[8, 11] <- test_amat[11, 8] <- 2
test_amat[9, 11] <- 1
test_that("Additional pra tests (1)", {
  # Mutations to test_amat inside a test_that() block don't persist to
  # later test_that() blocks (each one re-reads the file-scope matrix
  # above), so "(2)" and "(3)" below do not build cumulatively on this one.
  expect_output(
    result1 <- parentRecoveryAccuracy(
      test_amat, true_amat, c(2, 6, 10),
      verbose = TRUE
    )
  )
  expect_equal(
    result1, list(missing = 3, added = 1, correct = 1, potential = 0)
  )

  # make an undirected edge where the parent should be
  # adds a potential parent and removes a strictly missing parent
  test_amat[11, 8] <- test_amat[8, 11] <- 1
  expect_output(
    result2 <- parentRecoveryAccuracy(
      test_amat, true_amat, c(2, 6, 10),
      verbose = TRUE
    )
  )
  expect_equal(
    result2, list(missing = 2, added = 1, correct = 1, potential = 1)
  )
})

test_that("Additional pra tests (2)", {
  # test_amat here starts fresh from the file-scope matrix above (see the
  # note in "Additional pra tests (1)") -- these mutations do not build on
  # that earlier block's.
  # make a few adjustments
  # 11 ->(anc) 8 | missing parent (ancestral edge)
  test_amat[11, 8] <- 3
  test_amat[8, 11] <- 2
  # changes to 7 -> 3 (incorrect)
  test_amat[7, 3] <- 1
  test_amat[3, 7] <- 0
  # changes to undirected (correct)
  test_amat[3, 4] <- 1
  expect_output(
    result3 <- parentRecoveryAccuracy(
      test_amat, true_amat, c(2, 6, 10),
      verbose = TRUE
    )
  )
  expect_equal(
    result3, list(missing = 3, added = 1, correct = 1, potential = 0)
  )
})

test_that("Additional pra tests (3)", {
  # Make true parent an undirected edge (potential)
  test_amat[11, 8] <- test_amat[8, 11] <- 1
  expect_output(
    result4 <- parentRecoveryAccuracy(
      test_amat, true_amat, c(2, 6, 10),
      verbose = TRUE
    )
  )
  expect_equal(
    result4, list(missing = 2, added = 1, correct = 1, potential = 1)
  )
})


# Advanced V-Structure Testing -----------------------------------------------

test_that("additional v-structure comparison function tests", {
  true_amat <- matrix(c(
    0, 1, 1, 0, 0, 1,
    0, 0, 0, 0, 0, 1,
    0, 0, 0, 1, 0, 0,
    0, 1, 1, 0, 1, 1,
    0, 0, 0, 1, 0, 1,
    0, 0, 0, 1, 0, 0
  ), nrow = 6, byrow = TRUE)
  false_amat <- matrix(c(
    0, 0, 0, 0, 0, 1,
    1, 0, 0, 1, 0, 1,
    1, 0, 0, 0, 0, 0,
    0, 1, 1, 0, 0, 1,
    0, 0, 1, 0, 0, 1,
    0, 0, 0, 1, 0, 0
  ), nrow = 6, byrow = TRUE)

  expect_equal(
    compareVStructures(false_amat, true_amat, TRUE),
    list("missing" = 1, "added" = 2, "correct" = 2)
  )

  # change 3 - 4 - 5 in true graph to 3 -> 4 <- 5
  # (this combo is now incorrect and a false positive in estimated graph)
  true_amat[4, 3] <- true_amat[4, 5] <- 0
  expect_equal(
    compareVStructures(false_amat, true_amat, TRUE),
    list("missing" = 2, "added" = 2, "correct" = 2)
  )

  # Create additional v-structure by dropping 2 -> 1
  # also drops incorrect v-structure 2 -> 1 <- 3 from estimated graph
  true_amat[1, 2] <- 0
  expect_equal(
    compareVStructures(false_amat, true_amat, TRUE),
    list("missing" = 2, "added" = 2, "correct" = 2)
  )

  true_amat <- matrix(c(
    0, 0, 0, 0,
    1, 0, 0, 0,
    1, 0, 0, 0,
    0, 1, 0, 0
  ), byrow = TRUE, ncol = 4)
  expect_equal(
    compareVStructures(true_amat, true_amat, TRUE),
    list("missing" = 0, "added" = 0, "correct" = 1)
  )
})


# Ancestral Recovery Tests ----------------------------------------------------

test_that("Ancestral Relations", {
  # Basic, nothing informative (verbose defaults to FALSE, so no output)
  result <- interNeighborhoodEdgeMetrics(
    asiaDAG, asiaDAG, asiaDAG, seq(0, ncol(asiaDAG) - 1)
  )
  expect_equal(result, list(
    CorrectAncestors = 0, IncorrectAncestors = 0, TotalAncestralEdges = 0
  ))
})

test_that("Ancestral checks for cml on asia", {
  res <- cml(
    true_dag = asiaDAG, targets = c(1, 8), node_names = nodes,
    verbose = FALSE
  )
  # asia and dysp are targets: correct edge between tub and either, but it
  # is not oriented (verbose defaults to FALSE, so no output)
  result <- interNeighborhoodEdgeMetrics(
    res$amat, asiaDAG, asiaDAG, seq(0, ncol(asiaDAG) - 1)
  )
  expect_equal(result, list(
    CorrectAncestors = 0, IncorrectAncestors = 0, TotalAncestralEdges = 0
  ))
})

test_that("Ancestral checks for cml on asia (2)", {
  res <- cml(
    true_dag = asiaDAG, targets = c(3, 7), node_names = nodes,
    verbose = FALSE
  )
  # smoke and xray are the targets: lung and either are connected, but not
  # oriented (verbose defaults to FALSE, so no output)
  result <- interNeighborhoodEdgeMetrics(
    res$amat, asiaDAG, asiaDAG, seq(0, ncol(asiaDAG) - 1)
  )
  expect_equal(result, list(
    CorrectAncestors = 0, IncorrectAncestors = 0, TotalAncestralEdges = 0
  ))
})

true_amat <- matrix(0, ncol = 16, nrow = 16)
true_amat[1, 3] <- 1
true_amat[2, 4] <- true_amat[2, 11] <- 1
true_amat[1, 4] <- true_amat[4, 3] <- true_amat[4, 9] <- 1
true_amat[5, 4] <- 1
true_amat[6, 5] <- true_amat[6, 14] <- 1
true_amat[7, 4] <- true_amat[7, 8] <- 1
true_amat[9, 8] <- 1
true_amat[10, 5] <- 1
true_amat[11, 1] <- 1
true_amat[12, 16] <- 1
true_amat[13, 14] <- true_amat[13, 12] <- true_amat[13, 15] <- 1
true_amat[14, 10] <- 1
true_amat[15, 7] <- 1
true_amat[16, 11] <- 1

est_amat <- matrix(0, ncol = 16, nrow = 16)
est_amat[1, 12] <- 1
est_amat[2, 14] <- est_amat[14, 2] <- 1
est_amat[3, 1] <- 1
est_amat[4, 3] <- est_amat[4, 2] <- est_amat[4, 15] <- est_amat[15, 4] <- 1
est_amat[5, 4] <- est_amat[5, 15] <- est_amat[5, 14] <- est_amat[14, 5] <- 1
est_amat[7, 3] <- est_amat[7, 8] <- 1
est_amat[8, 9] <- est_amat[9, 8] <- 1
est_amat[12, 13] <- 1
est_amat[13, 14] <- est_amat[13, 9] <- 1
est_amat[15, 13] <- 1

test_that("Testing ancestral relations (1)", {
  # ancestral edge with no ancestry added
  est_amat[5, 14] <- est_amat[14, 5] <- 4
  # Incorrect ancestral edge
  est_amat[2, 14] <- 2
  est_amat[14, 2] <- 3
  # Correct ancestral edge
  est_amat[13, 9] <- 2
  est_amat[9, 13] <- 3
  num_nodes <- ncol(est_amat)
  expect_output(
    result <- interNeighborhoodEdgeMetrics(
      est_amat, true_amat, true_amat, seq(0, num_nodes - 1),
      verbose = TRUE
    )
  )
  expect_equal(result, list(
    CorrectAncestors = 1, IncorrectAncestors = 1, TotalAncestralEdges = 3
  ))
})

test_that("Testing ancestral relations (2)", {
  # est_amat here starts fresh from the file-scope matrix above -- mutations
  # made inside a test_that() block don't persist to later test_that()
  # blocks, so this does NOT build on the previous block's mutations.
  est_amat[12, 1] <- 2
  est_amat[1, 12] <- 4
  expect_output(
    result <- interNeighborhoodEdgeMetrics(
      est_amat, true_amat, true_amat, seq(0, ncol(est_amat) - 1),
      verbose = TRUE
    )
  )
  expect_equal(result, list(
    CorrectAncestors = 0, IncorrectAncestors = 0, TotalAncestralEdges = 1
  ))
})

test_that("Testing ancestral relations (3)", {
  # This should be overlooked because of the open mark
  est_amat[12, 1] <- 2

  # Get one true positive
  est_amat[1, 12] <- 3
  expect_output(
    result <- interNeighborhoodEdgeMetrics(
      est_amat, true_amat, true_amat, seq(0, ncol(est_amat) - 1),
      verbose = TRUE
    )
  )
  expect_equal(result, list(
    CorrectAncestors = 1, IncorrectAncestors = 0, TotalAncestralEdges = 1
  ))
})

true_amat <- matrix(c(
  0, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 1, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 0, 0, 0,
  0, 0, 0, 0, 0, 1, 0, 1,
  0, 0, 0, 0, 0, 0, 0, 0
), byrow = TRUE, ncol = 8)

est_amat <- matrix(c(
  0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 0, 0, 0, 0, 0, 0,
  0, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 1, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 4, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 2, 0, 0, 1,
  0, 0, 0, 0, 0, 0, 0, 0
), byrow = TRUE, ncol = 8)

est_amat2 <- matrix(c(
  0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 0, 0, 0, 0, 0, 0,
  0, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 1, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 3, 1, 0,
  0, 0, 0, 0, 2, 0, 0, 0,
  0, 0, 0, 0, 0, 1, 0, 1,
  0, 0, 0, 0, 0, 0, 0, 0
), byrow = TRUE, ncol = 8)

test_that("Testing ancestral relations (3b)", {
  result <- interNeighborhoodEdgeMetrics(
    est_amat, true_amat, true_amat, seq(0, ncol(true_amat) - 1),
    verbose = TRUE
  )
  expect_equal(result, list(
    CorrectAncestors = 0, IncorrectAncestors = 0, TotalAncestralEdges = 1
  ))
})

test_that("Testing ancestral relations (4)", {
  est_amat[5, 7] <- 3
  expect_output(
    result <- interNeighborhoodEdgeMetrics(
      est_amat, true_amat, true_amat, seq(0, ncol(true_amat) - 1),
      verbose = TRUE
    )
  )
  expect_equal(result, list(
    CorrectAncestors = 1, IncorrectAncestors = 0, TotalAncestralEdges = 1
  ))
})

test_that("Testing ancestral relations (5)", {
  expect_output(
    result <- interNeighborhoodEdgeMetrics(
      est_amat2, true_amat, true_amat, seq(0, ncol(true_amat) - 1),
      verbose = TRUE
    )
  )
  expect_equal(result, list(
    CorrectAncestors = 1, IncorrectAncestors = 0, TotalAncestralEdges = 1
  ))
})

test_that("Testing ancestral relations (6)", {
  est_amat2[5, 6] <- 2
  est_amat2[4, 1] <- 2
  est_amat2[1, 4] <- 3
  expect_output(
    result <- interNeighborhoodEdgeMetrics(
      est_amat2, true_amat, true_amat, seq(0, ncol(true_amat) - 1),
      verbose = TRUE
    )
  )
  expect_equal(result, list(
    CorrectAncestors = 0, IncorrectAncestors = 1, TotalAncestralEdges = 2
  ))
})

t_amat <- matrix(c(
  0, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 0, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 1,
  0, 0, 1, 0, 0, 0, 0, 0
), byrow = TRUE, nrow = 8)

test_that("Completing ancestral checks", {
  est_amat <- t_amat
  est_amat[4, 3] <- 2
  est_amat[3, 4] <- 3
  est_amat[, 8] <- est_amat[8, ] <- est_amat[7, ] <- est_amat[, 7] <- rep(0, 8)
  expect_output(
    result <- interNeighborhoodEdgeMetrics(
      est_amat, t_amat, t_amat, seq(0, 7),
      verbose = TRUE
    )
  )
  expect_equal(result, list(
    CorrectAncestors = 1, IncorrectAncestors = 0, TotalAncestralEdges = 1
  ))
})

t_amat <- matrix(c(
  0, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 1, 0, 0, 0, 0, 0, 1,
  0, 0, 0, 0, 1, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 1, 0
), byrow = TRUE, nrow = 8)

test_that("Completing ancestral checks (2)", {
  est_amat <- t_amat
  est_amat[3, 4] <- 2
  est_amat[4, 3] <- 3
  est_amat[, 8] <- est_amat[8, ] <- est_amat[7, ] <- est_amat[, 7] <- rep(0, 8)
  expect_output(
    result <- interNeighborhoodEdgeMetrics(
      est_amat, t_amat, t_amat, seq(0, 7),
      verbose = TRUE
    )
  )
  expect_equal(result, list(
    CorrectAncestors = 1, IncorrectAncestors = 0, TotalAncestralEdges = 1
  ))
})


# Testing Overall F1 Score ----------------------------------------------------

true_amat <- matrix(c(
  0, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 1, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 0, 0, 0,
  0, 0, 0, 0, 0, 1, 0, 1,
  0, 0, 0, 0, 0, 0, 0, 0
), byrow = TRUE, ncol = 8)

est_amat <- matrix(c(
  0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 0, 0, 0, 0, 0, 0,
  0, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 1, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 0, 0, 1,
  0, 0, 0, 0, 0, 0, 0, 0
), byrow = TRUE, ncol = 8)

est_amat2 <- matrix(c(
  0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 0, 0, 0, 0, 0, 0,
  0, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 1, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 0, 0, 0,
  0, 0, 0, 0, 0, 1, 0, 1,
  0, 0, 0, 0, 0, 0, 0, 0
), byrow = TRUE, ncol = 8)

test_that("Testing Overall F1 Score Function", {
  expect_output(
    result <- overallF1(est_amat, true_amat, c(0, 3, 6), verbose = TRUE)
  )
  expect_equal(result, 4 / 7)
})

test_that("Testing Overall F1 (2)", {
  expect_output(
    result <- overallF1(est_amat2, true_amat, c(0, 3, 6), verbose = TRUE)
  )
  expect_equal(result, 0.75)
})

test_that("Testing Overall F1 (3)", {
  est_amat3 <- est_amat2
  est_amat3[8, 6] <- est_amat3[6, 8] <- 1
  est_amat3[6, 7] <- 3
  est_amat3[7, 6] <- 2
  est_amat3[8, 7] <- 3
  est_amat3[7, 8] <- 2
  est_amat3[6, 8] <- est_amat3[8, 6] <- 3
  expect_output(
    result <- overallF1(est_amat3, true_amat, c(0, 3, 6), verbose = TRUE)
  )
  expect_equal(result, 2 / 7)
})


# All metric functions --------------------------------------------------------

test_that("checking metric functions for cml", {
  t <- c(1, 6, 7, 8)
  est <- cml(data = asiadf, true_dag = asiaDAG, targets = t, verbose = FALSE)

  asia_dag_sub <- asiaDAG
  asia_dag_sub[3, ] <- asia_dag_sub[, 3] <- rep(0, 8)

  # skeleton perfect, missing tub -> either <- lung and either -> dysp <-
  # bronc, have all parents except for tub and either, which are both
  # potential
  result <- allMetrics(
    est$amat, asia_dag_sub, t - 1, asiaDAG, seq(0, ncol(asiaDAG) - 1),
    algo = "cml"
  )
  expect_equal(result, data.frame(
    cml__skel_fp = 0, cml__skel_fn = 2, cml__skel_tp = 4,
    cml__v_fn = 2, cml__v_fp = 0, cml__v_tp = 0,
    cml_pra_fn = 2, cml_pra_fp = 0, cml_pra_tp = 3, cml_pra_potential = 0,
    cml_ancestors_correct = 0, cml_ancestors_incorrect = 0,
    cml_ancestors_total = 2, cml_overall_f1 = 0.8
  ))
})

test_that("checking metric functions for pc", {
  skip_if_not_installed("pcalg")
  t <- c(1, 6, 7, 8)
  est <- cml(data = asiadf, true_dag = asiaDAG, targets = t, verbose = FALSE)
  pc_fit <- as(pcalg::pc(
    suffStat = list(C = cor(asiadf), n = nrow(asiadf)),
    indepTest = pcalg::gaussCItest, ## indep.test: partial correlations
    alpha = 0.05, labels = colnames(asiaDAG),
    verbose = FALSE, m.max = 3
  ), "amat")
  pc_asia <- matrix(pc_fit, nrow = 8)
  # Remove smoke edges from consideration
  pc_asia[3, ] <- pc_asia[, 3] <- rep(0, 8)

  asia_dag_sub <- asiaDAG
  asia_dag_sub[3, ] <- asia_dag_sub[, 3] <- rep(0, 8)

  # skeleton perfect (smoke edges don't count), missing both v-structures
  # and added 1, missing all parents, either got two fp parents
  result <- allMetrics(
    pc_asia, asia_dag_sub, t - 1, asiaDAG, seq(0, ncol(asiaDAG) - 1),
    algo = "pc"
  )
  expect_equal(result, data.frame(
    pc__skel_fp = 0, pc__skel_fn = 0, pc__skel_tp = 6,
    pc__v_fn = 2, pc__v_fp = 1, pc__v_tp = 0,
    pc_pra_fn = 5, pc_pra_fp = 2, pc_pra_tp = 0, pc_pra_potential = 0,
    pc_ancestors_correct = 0, pc_ancestors_incorrect = 0,
    pc_ancestors_total = 0, pc_overall_f1 = 0
  ))
})

test_that("Check allMetrics (1)", {
  node_names_interest <- c("asia", "tub", "bronc", "xray", "dysp")
  nodes_int <- sapply(node_names_interest, function(x) which(x == nodes))
  names(nodes_int) <- NULL
  # targets are asia and bronc
  targets <- c(0, 2)
  est_mat <- matrix(c(
    0, 1, 0, 0, 0,
    0, 0, 0, 2, 2,
    0, 0, 0, 0, 0,
    0, 3, 1, 0, 1,
    0, 3, 1, 1, 0
  ), nrow = 5, byrow = TRUE)
  asia_g <- bnlearn::empty.graph(nodes)
  bnlearn::amat(asia_g) <- asiaDAG
  asia_cpdag <- bnlearn::cpdag(asia_g)
  asia_cpdag_mat <- bnlearn::amat(asia_cpdag)
  asia_cpdag_submat <- asia_cpdag_mat[nodes_int, nodes_int]
  result <- allMetrics(
    est_mat, asia_cpdag_submat, targets, asiaDAG, nodes_int - 1
  )
  expect_equal(result, data.frame(
    pc__skel_fp = 2, pc__skel_fn = 0, pc__skel_tp = 2,
    pc__v_fn = 0, pc__v_fp = 0, pc__v_tp = 0,
    pc_pra_fn = 0, pc_pra_fp = 2, pc_pra_tp = 0, pc_pra_potential = 0,
    pc_ancestors_correct = 2, pc_ancestors_incorrect = 0,
    pc_ancestors_total = 2, pc_overall_f1 = 0
  ))
})

test_that("Check allMetrics (2)", {
  node_names_interest <- c("asia", "tub", "bronc", "xray", "dysp")
  nodes_int <- sapply(node_names_interest, function(x) which(x == nodes))
  names(nodes_int) <- NULL
  # targets are asia and bronc
  targets <- c(0, 2)
  est_mat <- matrix(c(
    0, 1, 0, 0, 0,
    0, 0, 0, 2, 2,
    0, 0, 0, 0, 0,
    0, 3, 1, 0, 1,
    0, 3, 1, 1, 0
  ), nrow = 5, byrow = TRUE)
  asia_g <- bnlearn::empty.graph(nodes)
  bnlearn::amat(asia_g) <- asiaDAG
  asia_cpdag <- bnlearn::cpdag(asia_g)
  asia_cpdag_mat <- bnlearn::amat(asia_cpdag)
  asia_cpdag_submat <- asia_cpdag_mat[nodes_int, nodes_int]

  # Switch one ancestral edge
  est_mat[5, 2] <- 2
  est_mat[2, 5] <- 3
  result <- allMetrics(
    est_mat, asia_cpdag_submat, targets, asiaDAG, nodes_int - 1
  )
  expect_equal(result, data.frame(
    pc__skel_fp = 2, pc__skel_fn = 0, pc__skel_tp = 2,
    pc__v_fn = 0, pc__v_fp = 0, pc__v_tp = 0,
    pc_pra_fn = 0, pc_pra_fp = 2, pc_pra_tp = 0, pc_pra_potential = 0,
    pc_ancestors_correct = 1, pc_ancestors_incorrect = 1,
    pc_ancestors_total = 2, pc_overall_f1 = 0
  ))
})

test_that("Check allMetrics (3)", {
  node_names_interest <- c("asia", "tub", "bronc", "xray", "dysp")
  nodes_int <- sapply(node_names_interest, function(x) which(x == nodes))
  names(nodes_int) <- NULL
  # targets are asia and bronc
  targets <- c(0, 2)
  est_mat <- matrix(c(
    0, 1, 0, 0, 0,
    0, 0, 0, 2, 2,
    0, 0, 0, 0, 0,
    0, 3, 1, 0, 1,
    0, 3, 1, 1, 0
  ), nrow = 5, byrow = TRUE)
  asia_g <- bnlearn::empty.graph(nodes)
  bnlearn::amat(asia_g) <- asiaDAG
  asia_cpdag <- bnlearn::cpdag(asia_g)
  asia_cpdag_mat <- bnlearn::amat(asia_cpdag)
  asia_cpdag_submat <- asia_cpdag_mat[nodes_int, nodes_int]

  # Switch one ancestral edge
  est_mat[5, 2] <- 2
  est_mat[2, 5] <- 3

  # Ancestral edge without orientation
  est_mat[4, 2] <- est_mat[2, 4] <- 2
  result <- allMetrics(
    est_mat, asia_cpdag_submat, targets, asiaDAG, nodes_int - 1
  )
  expect_equal(result, data.frame(
    pc__skel_fp = 2, pc__skel_fn = 0, pc__skel_tp = 2,
    pc__v_fn = 0, pc__v_fp = 0, pc__v_tp = 0,
    pc_pra_fn = 0, pc_pra_fp = 2, pc_pra_tp = 0, pc_pra_potential = 0,
    pc_ancestors_correct = 0, pc_ancestors_incorrect = 1,
    pc_ancestors_total = 2, pc_overall_f1 = 0
  ))
})


# Miscellaneous tests -----------------------------------------------------

test_that("Miscellaneous function tests", {
  expect_equal(
    getNeighborhoodMetrics(asiaDAG), data.frame("size" = 8, "num_edges" = 8)
  )
  expect_equal(
    getNeighborhoodMetrics(andesDAG),
    data.frame("size" = 223, "num_edges" = 338)
  )
})

test_that("Testing warnings and stops", {
  m4 <- matrix(0, nrow = 4, ncol = 4)
  m4_3col <- matrix(0, nrow = 4, ncol = 3)
  m5 <- matrix(0, nrow = 5, ncol = 5)
  m6 <- matrix(0, nrow = 6, ncol = 6)

  # Invalid matrices
  expect_error(compareSkeletons(m4_3col, m4, 0))
  expect_error(compareSkeletons(m4, m4_3col, 0))
  expect_error(compareSkeletons(m5, m6, 1))

  # Invalid targets
  expect_error(parentRecoveryAccuracy(m5, m5, -1))
  expect_error(parentRecoveryAccuracy(m5, m5, 5))
  expect_error(parentRecoveryAccuracy(m5, m5, 6))
  expect_error(parentRecoveryAccuracy(m5, m5, c(2, 5, 4)))
  expect_error(parentRecoveryAccuracy(m5, m5, c(2, 4, -1)))

  # We have an undirected edge denoted by 3's for both entries
  test_amat <- asiaDAG
  test_amat[2, 8] <- test_amat[8, 2] <- 3
  expect_warning(
    interNeighborhoodEdgeMetrics(
      test_amat, asiaDAG, asiaDAG, seq(0, ncol(asiaDAG) - 1)
    )
  )
})

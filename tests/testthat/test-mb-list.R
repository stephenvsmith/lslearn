### Testing Markov Blanket List class

nodes <- seq(0, 6)
mb_mat <- matrix(c(
  0, 0, 0, 1, 0, 0, 0, # node 0's MB is 3
  0, 0, 1, 0, 0, 1, 0, # node 1's MB is 2, 5
  0, 1, 0, 0, 0, 1, 1, # node 2's MB is 1, 5, 6
  1, 0, 0, 0, 1, 0, 1, # node 3's MB is 0, 4, 6
  0, 0, 0, 1, 0, 0, 0, # node 4's MB is 3
  0, 1, 1, 0, 0, 0, 0, # node 5's MB is 1, 2
  0, 0, 1, 1, 0, 0, 0 # node 6's MB is 2, 3
), byrow = TRUE, nrow = 7)
data("asiaDAG")

test_that("Initialization of Markov Blanket Object works (Sample)", {
  expect_output(testInitializeMBList(nodes, mb_mat))
})

test_that("Access Markov Blanket for a node", {
  expect_equal(testAccessMB(nodes, mb_mat, 0), c(3))
  expect_equal(testAccessMB(nodes, mb_mat, 2), c(1, 5, 6))
  expect_equal(testAccessMB(nodes, mb_mat, 1), c(2, 5))
  expect_error(testAccessMB(nodes, mb_mat, 50))
})

test_that("Access Markov Blanket for multiple nodes", {
  expect_equal(
    testAccessMultipleMB(nodes, mb_mat, c(0, 1, 2)), c(1, 2, 3, 5, 6)
  )
  expect_equal(
    testAccessMultipleMB(nodes, mb_mat, c(0, 2, 1)), c(1, 2, 3, 5, 6)
  )
  expect_equal(testAccessMultipleMB(nodes, mb_mat, numeric()), numeric())

  expect_warning(testAccessMultipleMB(nodes, mb_mat, c(0, 2, 1), TRUE, TRUE))
  expect_equal(
    testAccessMultipleMB(nodes, mb_mat, c(0, 2, 1), TRUE, FALSE, TRUE),
    c(0, 1, 2, 3, 5, 6)
  )
  expect_equal(
    testAccessMultipleMB(nodes, mb_mat, c(0, 2, 1), FALSE, TRUE, TRUE),
    c(3, 5, 6)
  )
})

test_that("Membership in MB function is correct", {
  # We should find 5 in MB(1) and 6 in MB(3)
  expect_true(testIsMBMember(nodes, mb_mat, 1, 5))
  expect_true(testIsMBMember(nodes, mb_mat, 3, 6))
  # 0 not in MB(4) and 2 not in MB(3)
  expect_false(testIsMBMember(nodes, mb_mat, 4, 0))
  expect_false(testIsMBMember(nodes, mb_mat, 3, 2))

  expect_error(testIsMBMember(nodes, mb_mat, 3, 15))
})

nodes <- seq(0, 7)
node_names <- colnames(asiaDAG)
test_that("Initialization for Population Version works", {
  expect_output(testInitializeMBListPop(nodes, asiaDAG))
})

test_that("Initialization for Population Version with fewer nodes", {
  expect_output(testInitializeMBListPop(c(0, 3, 5), asiaDAG))
})

test_that("Test silencer", {
  # Round 1 and Round 3 should print the MB traversal output below their
  # header; Round 2, after calling silencer(), should print nothing between
  # its header and the next one.
  out <- capture.output(testSilencer(c(0, 1, 3, 5), asiaDAG, 0, 5))
  round1 <- which(out == "Round 1:")
  round2 <- which(out == "Round 2:")
  round3 <- which(out == "Round 3:")
  expect_length(round1, 1)
  expect_length(round2, 1)
  expect_length(round3, 1)
  expect_gt(round2 - round1, 1)
  expect_equal(round3 - round2, 1)
  expect_gt(length(out) - round3, 0)
})

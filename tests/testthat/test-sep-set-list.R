# Example setup -------------------------------------------------------------

data("asiaDAG")
nodes <- colnames(asiaDAG)
target <- 5 # "either" is the target

# The first-order Markov blanket neighborhood of "either" in asiaDAG,
# precomputed via CML's DAG-based `getNbhd()` helper (DAG has not been
# ported to lslearn yet; see Phase 2).
neighbors <- c(1, 3, 4, 5, 6, 7)
neighbor_names <- nodes[neighbors + 1]
size <- length(neighbors) - 1

# Initialization of Separating Sets ------------------------------------------

test_that("Checking initializations", {
  # Check S for neighborhood around either
  expect_equal(nodes[target + 1], "either")
  expect_equal(
    setdiff(
      neighbor_names, c("tub", "lung", "xray", "dysp", "bronc", "either")
    ),
    character(0)
  )
  expect_output(printS(neighbors))
})

test_that("Check initial values", {
  # Check all initial sets are NA
  for (n in 0:size) {
    for (n2 in setdiff(0:size, n)) {
      expect_true(is.na(getInitialValues(neighbors, n, n2)))
    }
  }
})

# Check member functions ------------------------------------------------------

test_that("Check member functions", {
  # Check that we can designate an empty set in S with -1
  set.seed(111)
  neighborhood_shuffled <- sample(1:size, size, replace = FALSE)
  expect_equal(
    setListEmptySet(
      neighbors, neighborhood_shuffled[1], neighborhood_shuffled[2]
    ),
    -1
  )
  expect_equal(
    setListEmptySet(
      neighbors, neighborhood_shuffled[2], neighborhood_shuffled[5]
    ),
    -1
  )

  # Check that we can designate another set for S[i,j]
  i <- 0
  j <- 1
  kvals <- 4:5
  expect_equal(setListEfficient(neighbors, i, j, kvals), kvals)
  expect_equal(setListEfficient(neighbors, j, i, kvals), kvals)

  # Checking if nodes 0 and 1 are separated by node 2 when S[0,1]={4,5}
  expect_false(checkSeparationFunc(neighbors, i, j, kvals, 2))
  # Checking if nodes 0 and 1 are separated by node 5 when S[0,1]={4,5}
  # Should produce warning and should be false because both separation sets
  # are not correctly defined
  expect_warning(tmp <- checkSeparationFunc(neighbors, i, j, kvals, 5))
  expect_false(tmp)
  # Corrected Version (i.e. S[0,1]=S[1,0])
  expect_true(checkSeparationFuncCorrected(neighbors, i, j, kvals, 5))
})

test_that("Checks for errors and warnings in values", {
  # Negative values
  expect_error(setListEfficient(neighbors, i = -1, j = 0, c(1, 2)))
  expect_error(setListEfficient(neighbors, i = 0, j = -1, c(1, 2)))
  expect_error(setListEfficient(neighbors, i = 5, j = 0, c(1, 2)), NA)
  expect_error(setListEfficient(neighbors, i = 6, j = 0, c(1, 2)))
  expect_error(setListEfficient(neighbors, i = 3, j = 6, c(1, 2)))

  # Separation sets disagree
  expect_warning(
    checkIsSepSetMember(neighbors, 0, 1, c(3, 4), c(4, 5), 4),
    NA
  )
  expect_warning(checkIsSepSetMember(neighbors, 0, 1, c(3, 4), c(4, 5), 3))
})

test_that("Check for potential v-structure", {
  neighbors <- c(0, 1, 3, 4, 5)
  expect_true(checkPotentialVStruct(neighbors, 0, 1, c(4, 5), 3))
  expect_true(checkPotentialVStruct(neighbors, 0, 1, c(4, 5), 2))
  expect_false(checkPotentialVStruct(neighbors, 0, 1, c(4, 5), 5))
})

test_that("Check that we can return S", {
  neighbors <- c(0, 1, 3, 4, 5)
  expected_s <- list(
    list(NA_real_, c(4, 5), NA_real_, NA_real_, NA_real_),
    list(c(4, 5), NA_real_, NA_real_, NA_real_, NA_real_),
    list(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
    list(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
    list(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)
  )
  expect_equal(checkGetS(neighbors, 0, 1, c(4, 5)), expected_s)
})

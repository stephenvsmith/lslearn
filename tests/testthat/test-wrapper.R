# Setup -------------------------------------------------------------------

data("asiadf")
data("asiaDAG")
node_names <- colnames(asiaDAG)
asia_dag <- matrix(asiaDAG, nrow = ncol(asiadf), ncol = ncol(asiadf))
p <- length(node_names)


test_that("Wrapper function works (Sample with true DAG)", {
  # CML
  results <- cml(
    true_dag = asia_dag, data = asiadf, targets = c(1, 6),
    node_names = node_names, verbose = FALSE
  )
  for (result_part in setdiff(
    names(results),
    c(
      "totalSkeletonTime", "targetSkeletonTimes", "mbEstTime", "totalTime",
      "MBNumTests"
    )
  )) {
    expect_output(print(results[[result_part]]))
  }
  expect_equal(results[["MBNumTests"]], 0)

  # SNL
  results <- snl(
    true_dag = asia_dag, data = asiadf, targets = c(1, 6),
    node_names = node_names, verbose = FALSE
  )
  for (result_part in setdiff(
    names(results),
    c("targetSkeletonTimes", "totalTime", "mbEstTime", "MBNumTests")
  )) {
    expect_output(print(results[[result_part]]))
  }
  expect_equal(results[["MBNumTests"]], 0)
})

test_that("Wrapper function works (Population)", {
  # CML
  results <- cml(
    true_dag = asia_dag, targets = c(1, 6), node_names = node_names,
    verbose = FALSE
  )
  for (result_part in setdiff(
    names(results),
    c(
      "totalSkeletonTime", "targetSkeletonTimes", "mbEstTime", "totalTime",
      "MBNumTests"
    )
  )) {
    expect_output(print(results[[result_part]]))
  }
  expect_equal(results[["MBNumTests"]], 0)

  # SNL
  results <- snl(
    true_dag = asia_dag, targets = c(1, 6), node_names = node_names,
    verbose = FALSE
  )
  for (result_part in setdiff(
    names(results),
    c("targetSkeletonTimes", "totalTime", "mbEstTime", "MBNumTests")
  )) {
    expect_output(print(results[[result_part]]))
  }
  expect_equal(results[["MBNumTests"]], 0)
})

test_that("Wrapper function works (Sample)", {
  # CML
  results <- cml(
    data = asiadf, targets = c(1, 6), node_names = node_names,
    verbose = FALSE
  )
  for (result_part in setdiff(
    names(results),
    c(
      "totalSkeletonTime", "targetSkeletonTimes", "totalTime", "mbList",
      "mbEstTime", "MBNumTests"
    )
  )) {
    expect_output(print(results[[result_part]]))
  }
  results1 <- results[["MBNumTests"]]

  # SNL
  results <- snl(
    data = asiadf, targets = c(1, 6), node_names = node_names,
    verbose = FALSE
  )
  for (result_part in setdiff(
    names(results),
    c("targetSkeletonTimes", "totalTime", "mbList", "MBNumTests", "mbEstTime")
  )) {
    expect_output(print(results[[result_part]]))
  }
  expect_output(print(results1))
  expect_output(print(results[["MBNumTests"]]))
})

test_that("Testing pre-checks", {
  expect_error(cml(
    data = asiadf, targets = c(1, 6), lmax = -1, node_names = node_names,
    verbose = FALSE
  ))
  expect_error(snl(
    data = asiadf, targets = c(1, 6), lmax = -1, node_names = node_names,
    verbose = FALSE
  ))
})

test_that("Discrete version", {
  data("asia", package = "bnlearn")
  asia_mat <- data.matrix(asia) - 1

  asia_d_cml <- cml(
    asia_mat,
    targets = c(2, 5, 6), test = "gSquare", verbose = FALSE
  )
  expect_output(print(asia_d_cml$S))

  # Ground truth pulled from a real (unmodified) CML build -- MMPC/gSquare
  # Markov Blanket estimation over discrete data is deterministic, so this
  # is reproducible.
  expect_equal(asia_d_cml$amat, matrix(c(
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0,
    0, 1, 0, 0, 0, 1, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0
  ), nrow = 8, byrow = TRUE))
  expect_equal(asia_d_cml$RulesUsed, c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
})

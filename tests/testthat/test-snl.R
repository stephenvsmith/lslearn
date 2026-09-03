# Setup -------------------------------------------------------------------

data("asiaDAG")
data("asiadf")
t <- 5
node_names <- colnames(asiaDAG)
df <- as.matrix(asiadf)
asia_dag <- as.matrix(asiaDAG)
p <- ncol(asia_dag)

test_that("Targets are properly validated", {
  targets <- c(1, 2, 4, 5, 8)
  expect_error(validateTargetSNL(targets, 10))
  expect_error(validateTargetSNL(targets, 5), NA)
  expect_error(validateTargetSNL(targets, -1))
})


# Test SNL ------------------------------------------------------------------

test_that("Sample SNL (given true Markov Blankets)", {
  expect_equal(
    checkInitializeSNL(asia_dag, df, t, seq(0, p - 1), node_names), 8
  )

  pc_dag <- bnlearn::empty.graph(node_names[c(2, 4, 5, 6, 7, 8)])
  expect_output(
    bnlearn::amat(pc_dag) <- checkGetTargetSkel(
      asia_dag, df, t, seq(0, p - 1), node_names
    )
  )
  # tub, lung, bronc, either, xray, dysp
  expect_equal(bnlearn::amat(pc_dag), matrix(c(
    0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1,
    0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0
  ), nrow = 6, dimnames = list(
    node_names[c(2, 4, 5, 6, 7, 8)], node_names[c(2, 4, 5, 6, 7, 8)]
  )))
})

test_that("Check V-Structures", {
  expect_equal(
    checkInitializeSNL(asia_dag, df, t, seq(0, p - 1), node_names), 8
  )

  pc_dag <- bnlearn::empty.graph(node_names[c(2, 4, 5, 6, 7, 8)])
  bnlearn::amat(pc_dag) <- checkGetTargetSkel(
    asia_dag, df, t, seq(0, p - 1), node_names
  )

  expect_output(
    bnlearn::amat(pc_dag) <- checkGetVStructures(
      asia_dag, df, t, seq(0, p - 1), node_names
    )
  )
  # tub, lung, bronc, either, xray, dysp: only the tub->either<-lung
  # v-structure is detected here (bronc->dysp<-either isn't, since bronc and
  # either remain adjacent at this point in the skeleton search)
  expect_equal(bnlearn::amat(pc_dag), matrix(c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0
  ), nrow = 6, dimnames = list(
    node_names[c(2, 4, 5, 6, 7, 8)], node_names[c(2, 4, 5, 6, 7, 8)]
  )))
})

test_that("Population SNL", {
  pc_dag_pop <- bnlearn::empty.graph(node_names[c(2, 4, 5, 6, 7, 8)])
  expect_equal(checkInitializeSNLPop(asia_dag, t, seq(0, p - 1), node_names), 8)
  bnlearn::amat(pc_dag_pop) <- checkGetTargetSkelPop(
    asia_dag, t, seq(0, p - 1), node_names
  )
  # tub, lung, bronc, either, xray, dysp
  expect_equal(
    bnlearn::amat(pc_dag_pop),
    matrix(c(
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1,
      1, 1, 0, 0, 1, 1,
      0, 0, 0, 1, 0, 0,
      0, 0, 1, 1, 0, 0
    ), nrow = 6, byrow = TRUE, dimnames = list(
      node_names[c(2, 4, 5, 6, 7, 8)], node_names[c(2, 4, 5, 6, 7, 8)]
    ))
  )

  expect_output(
    bnlearn::amat(pc_dag_pop) <- checkGetVStructuresPop(
      asia_dag, t, seq(0, p - 1), node_names
    )
  )
  expect_equal(bnlearn::amat(pc_dag_pop), matrix(c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0
  ), nrow = 6, dimnames = list(
    node_names[c(2, 4, 5, 6, 7, 8)], node_names[c(2, 4, 5, 6, 7, 8)]
  )))
})


# Test SNL Multiple Targets -------------------------------------------------

test_that("Initializing for multiple targets", {
  expect_output(
    checkInitializeSNLPop(asia_dag, c(0, 7), seq(0, p - 1), node_names)
  )
})

test_that("Multiple Targets", {
  amat_test <- checkSNL(asia_dag, df, c(3, 4), seq(0, p - 1), node_names)
  expect_equal(amat_test, matrix(c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0
  ), nrow = 8))
})

test_that("Multiple Targets Population Version", {
  amat_pop <- checkSNLPop(asia_dag, c(3, 4), seq(0, p - 1), node_names)
  expect_equal(amat_pop, matrix(c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0
  ), nrow = 8))
})

test_that("Additional Pop. Version test", {
  amat_pop2 <- checkSNLPop(asia_dag, c(0, 5), seq(0, p - 1), node_names)
  expect_equal(amat_pop2, matrix(c(
    0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0
  ), nrow = 8))
})


# Test Wrapper Functions -----------------------------------------------------

test_that("Wrapper function (Population Version)", {
  expect_output(
    wrapper_res <- snl(
      true_dag = asia_dag, targets = c(3, 4), node_names = node_names
    )
  )
  expect_named(wrapper_res, c(
    "amat", "S", "NumTests", "MBNumTests", "Nodes", "targetSkeletonTimes",
    "totalTime", "referenceDAG", "mbList", "rules_used", "data_means",
    "data_cov"
  ))
})

test_that("Wrapper function (Semi-Sample Version)", {
  expect_output(
    wrapper_res_semi <- snl(
      data = asiadf, true_dag = asia_dag, targets = c(3, 4),
      node_names = node_names
    )
  )
  expect_equal(wrapper_res_semi$MBNumTests, 0)
})

test_that("Wrapper function (Sample Version)", {
  expect_output(
    wrapper_res_sample <- snl(
      data = asiadf, targets = c(3, 4), node_names = node_names
    )
  )
  # Same as above, but without names
  expect_output(wrapper_res_sample <- snl(data = asiadf, targets = c(3, 4)))
  expect_true(wrapper_res_sample$MBNumTests > 0)
})

test_that("Run function (Sample Version)", {
  expect_output(checkSNLRun(asia_dag, df, c(3, 4), seq(0, p - 1), node_names))
})


# Testing Rules ---------------------------------------------------------------

test_that("Testing rule 1", {
  my_amat <- matrix(c(
    0, 1, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 1,
    0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 1, 0, 1, 0
  ), byrow = TRUE, nrow = 8)

  expect_output(
    amat_end <- checkRule1(asia_dag, my_amat, 0:7, 0:7, as.character(1:8))
  )
  expect_equal(amat_end, matrix(c(
    0, 1, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 1,
    0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 1, 0, 1, 0
  ), byrow = TRUE, nrow = 8))
})

test_that("Test rule 1 no change", {
  expect_output(
    res <- checkRule1(asia_dag, asia_dag, 0:7, 0:7, as.character(1:8))
  )
  expect_equal(res, asia_dag)
})

test_that("Testing rule 1 (2)", {
  my_amat <- matrix(c(
    0, 1, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1, 0, 1, 1,
    0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 1, 0, 1, 0
  ), byrow = TRUE, nrow = 8)

  expect_output(checkRule1(asia_dag, my_amat, 0:7, 0:7, as.character(1:8)))
})

test_that("Testing rule 2", {
  my_amat <- matrix(c(
    0, 1, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 1,
    0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 1, 0, 1, 0
  ), byrow = TRUE, nrow = 8)
  # Add an edge 1 - 8 and 1 - 6
  my_amat[1, 8] <- my_amat[8, 1] <- 1
  my_amat[1, 6] <- my_amat[6, 1] <- 1

  expect_output(
    amat_end <- checkRule2(asia_dag, my_amat, 0:7, 0:7, as.character(1:8))
  )
  expect_equal(amat_end, matrix(c(
    0, 1, 1, 0, 0, 1, 0, 1,
    0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 1,
    0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 1, 0, 1, 0
  ), byrow = TRUE, nrow = 8))
})

test_that("Test rule 2 no change", {
  expect_output(
    res <- checkRule2(asia_dag, asia_dag, 0:7, 0:7, as.character(1:8))
  )
  expect_equal(res, asia_dag)
})

test_that("Testing rule 3", {
  my_amat <- matrix(c(
    0, 1, 1, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 1, 1, 1, 0,
    1, 0, 0, 1, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 1,
    0, 1, 0, 0, 0, 0, 1, 0,
    0, 1, 0, 0, 0, 0, 1, 0,
    0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0
  ), byrow = TRUE, nrow = 8)

  expect_output(
    amat_end <- checkRule3(asia_dag, my_amat, 0:7, 0:7, as.character(1:8))
  )
  expect_equal(amat_end, matrix(c(
    0, 1, 1, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 1, 1, 1, 0,
    1, 0, 0, 1, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 1,
    0, 1, 0, 0, 0, 0, 1, 0,
    0, 1, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0
  ), byrow = TRUE, nrow = 8))
})

test_that("Test rule 3 no change", {
  expect_output(
    res <- checkRule3(asia_dag, asia_dag, 0:7, 0:7, as.character(1:8))
  )
  expect_equal(res, asia_dag)
})

test_that("Testing rule 4", {
  my_amat <- matrix(c(
    0, 1, 1, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 1, 1, 1, 0,
    1, 0, 0, 1, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 1,
    0, 1, 0, 0, 0, 0, 1, 0,
    0, 1, 0, 0, 0, 0, 1, 0,
    0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0
  ), byrow = TRUE, nrow = 8)
  my_amat[4, 3] <- 1
  my_amat[1, 4] <- my_amat[4, 1] <- 0
  my_amat[2, 3] <- 1

  expect_output(
    amat_end <- checkRule4(asia_dag, my_amat, 0:7, 0:7, as.character(1:8))
  )
  expect_equal(amat_end, matrix(c(
    0, 1, 1, 0, 0, 0, 0, 0,
    0, 0, 1, 1, 1, 1, 1, 0,
    1, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1,
    0, 1, 0, 0, 0, 0, 1, 0,
    0, 1, 0, 0, 0, 0, 1, 0,
    0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0
  ), byrow = TRUE, nrow = 8))
})

test_that("Test rule 4 no change", {
  expect_output(
    res <- checkRule4(asia_dag, asia_dag, 0:7, 0:7, as.character(1:8))
  )
  expect_equal(res, asia_dag)
})

test_that("Testing all rules together", {
  my_amat <- matrix(c(
    0, 1, 1, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 1, 1, 1, 0,
    1, 0, 0, 1, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 1,
    0, 1, 0, 0, 0, 0, 1, 0,
    0, 1, 0, 0, 0, 0, 1, 0,
    0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0
  ), byrow = TRUE, nrow = 8)
  my_amat[4, 3] <- 1
  my_amat[1, 4] <- my_amat[4, 1] <- 0
  my_amat[2, 3] <- my_amat[3, 2] <- 1
  my_amat[2, 8] <- my_amat[8, 2] <- 1

  expect_output(
    amat_end <- checkSNLRules(asia_dag, my_amat, 0:7, 0:7, as.character(1:8))
  )
  expect_equal(amat_end, matrix(c(
    0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0
  ), nrow = 8))
})

# Misc. -----------------------------------------------------------------------

test_that("Testing data structure correctness", {
  dag_df <- data.frame(asia_dag)
  expect_error(tmp <- snl(true_dag = dag_df, targets = t), NA)
})

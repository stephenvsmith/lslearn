# Setup -------------------------------------------------------------------

data("asiadf")
data("asiaDAG")
node_names <- colnames(asiaDAG)
asia_dag <- matrix(asiaDAG, nrow = ncol(asiadf), ncol = ncol(asiadf))
asiadf <- as.matrix(asiadf)

# Rcpp matrices alias the R-level SEXP: setAmat()/rule*() mutate their input
# matrix in place, which would silently corrupt a variable still holding the
# "before" version unless it's deep-copied first.
deep_copy <- function(m) matrix(as.vector(m), nrow = nrow(m), ncol = ncol(m))


test_that("Rule 1 is correct.", {
  nodes <- c("a", "b", "alpha", "beta", "gamma")
  adj_mat1 <- matrix(c(
    0, 1, 1, 0, 0,
    1, 0, 0, 2, 0,
    1, 0, 0, 2, 0,
    0, 3, 1, 0, 1,
    0, 0, 0, 1, 0
  ), nrow = 5, byrow = TRUE)
  adj_mat1orig <- deep_copy(adj_mat1)
  expect_output(
    adj_mat1 <- testRule1(asia_dag[1:5, 1:5], asiadf[, 1:5], 0, nodes, adj_mat1)
  )
  expect_equal(adj_mat1[4, 5], 2)
  expect_equal(adj_mat1[5, 4], 3)

  for (i in seq_len(nrow(adj_mat1))) {
    for (j in seq_len(ncol(adj_mat1))) {
      cond1 <- i == 4 & j == 5
      cond2 <- i == 5 & j == 4
      if (!(cond1 | cond2)) {
        expect_equal(adj_mat1orig[i, j], adj_mat1[i, j])
      }
    }
  }
})

test_that("Rule 1 Contradiction Gives Error", {
  # If we have beta o- gamma, this is an error and a contradiction of rule 1
  # Since this would mean we would have another v-structure
  amat_error <- matrix(c(
    0, 2, 0,
    1, 0, 3,
    0, 1, 0
  ), byrow = TRUE, ncol = 3)
  expect_warning(testRule1(asia_dag, asiadf, 0, node_names, amat_error))
})

test_that("Rule 2 is correct.", {
  adj_mat2 <- matrix(c(
    0, 2, 1, 1, 0, 0,
    1, 0, 2, 0, 0, 0,
    1, 3, 0, 0, 0, 0,
    1, 0, 0, 0, 2, 1,
    0, 0, 0, 3, 0, 2,
    0, 0, 0, 1, 1, 0
  ), nrow = 6, byrow = TRUE)
  p <- nrow(adj_mat2)
  adj_mat2orig <- deep_copy(adj_mat2)
  expect_output(
    adj_mat2 <- testRule2(
      asia_dag[1:p, 1:p], asiadf[, 1:p], 0, node_names[1:p], adj_mat2
    )
  )
  expect_equal(adj_mat2[1, 3], 2)
  expect_equal(adj_mat2[4, 6], 2)

  for (i in seq_len(nrow(adj_mat2))) {
    for (j in seq_len(ncol(adj_mat2))) {
      cond1 <- i == 1 & j == 3
      cond2 <- i == 4 & j == 6
      if (!(cond1 | cond2)) {
        expect_equal(adj_mat2[i, j], adj_mat2orig[i, j])
      }
    }
  }
})

test_that("Rule 3 is correct.", {
  adj_mat3 <- matrix(c(
    0, 2, 2, 1, 0, 0, 1,
    1, 0, 0, 0, 0, 1, 0,
    3, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 2, 0, 1,
    0, 0, 0, 1, 0, 1, 1,
    0, 1, 0, 0, 2, 0, 1,
    1, 0, 0, 1, 1, 1, 0
  ), nrow = 7, byrow = TRUE)
  adj_mat3orig <- deep_copy(adj_mat3)
  p <- ncol(adj_mat3)

  expect_output(
    adj_mat3 <- testRule3(
      asia_dag[1:p, 1:p], asiadf[, 1:p], 0, node_names[1:p], adj_mat3
    )
  )
  expect_equal(adj_mat3[7, 5], 2)
  for (i in seq_len(nrow(adj_mat3))) {
    for (j in seq_len(ncol(adj_mat3))) {
      cond1 <- i == 7 & j == 5
      if (!cond1) {
        expect_equal(adj_mat3[i, j], adj_mat3orig[i, j])
      }
    }
  }
})


adj_mat4 <- matrix(c(
  0, 2, 1, 0, 0, 2, 3,
  2, 0, 0, 2, 0, 2, 0,
  3, 0, 0, 0, 3, 0, 1,
  0, 2, 0, 0, 1, 2, 0,
  0, 0, 2, 2, 0, 1, 0,
  3, 3, 0, 3, 1, 0, 0,
  2, 0, 1, 0, 0, 0, 0
), nrow = 7, byrow = TRUE)
adj_mat4_v2 <- deep_copy(adj_mat4)
adj_mat4_v2_1 <- deep_copy(adj_mat4)
adj_mat4_v2_1_1 <- deep_copy(adj_mat4)
adj_mat4_v3 <- deep_copy(adj_mat4)
adj_mat4orig <- deep_copy(adj_mat4)
p <- ncol(adj_mat4)

nodes <- c("a", "b", "c", "alpha", "beta", "gamma", "theta")
sep <- c(4)
test_that("Rule 4 is correct.", {
  expect_output(
    adj_mat4 <- testRule4(
      asia_dag[1:p, 1:p], asiadf[, 1:p], seq(0, p - 1), nodes,
      adj_mat4, 5, 6, sep
    )
  )
  expect_equal(adj_mat4[5, 6], 2)
  expect_equal(adj_mat4[6, 5], 3)
  expect_equal(adj_mat4orig[5, 6], 1)
  expect_equal(adj_mat4orig[6, 5], 1)

  for (i in seq_len(nrow(adj_mat4))) {
    for (j in seq_len(ncol(adj_mat4))) {
      cond1 <- (i == 5 & j == 6) | (i == 6 & j == 5)
      if (!cond1) {
        expect_equal(adj_mat4[i, j], adj_mat4orig[i, j])
      }
    }
  }
})

test_that("Rule 4 is correct (beta not in separation set).", {
  # Same as above, but this time beta is NOT in the separating set
  sep <- c(1)
  expect_output(
    adj_mat4_v2 <- testRule4(
      asia_dag[1:p, 1:p], asiadf[, 1:p], seq(0, p - 1), nodes,
      adj_mat4_v2, 5, 6, sep
    )
  )
  expect_equal(adj_mat4_v2[5, 6], 2)
  expect_equal(adj_mat4_v2[6, 5], 2)
  expect_equal(adj_mat4_v2[4, 5], 2)
  expect_equal(adj_mat4orig[5, 6], 1)
  expect_equal(adj_mat4orig[6, 5], 1)
  expect_equal(adj_mat4orig[4, 5], 1)

  for (i in seq_len(nrow(adj_mat4_v2))) {
    for (j in seq_len(ncol(adj_mat4_v2))) {
      cond1 <- (i == 5 & j == 6) | (i == 6 & j == 5)
      cond2 <- i == 4 & j == 5
      if (!(cond1 | cond2)) {
        expect_equal(adj_mat4_v2[i, j], adj_mat4orig[i, j])
      }
    }
  }
})

test_that("Rule 4 Multiple Options:", {
  adj_mat4_v2_1[3, 4] <- 2
  adj_mat4_v2_1[4, 3] <- 2
  expect_output(
    res <- testRule4(
      asia_dag[1:p, 1:p], asiadf[, 1:p], seq(0, p - 1), nodes,
      adj_mat4_v2_1, 5, 6, sep
    )
  )
})

test_that("Rule 4 Multiple Options (2):", {
  adj_mat4_v2_1_1[3, 4] <- 2
  adj_mat4_v2_1_1[4, 3] <- 2
  adj_mat4_v2_1_1[3, 6] <- 1
  adj_mat4_v2_1_1[6, 3] <- 3
  expect_output(
    res <- testRule4(
      asia_dag[1:p, 1:p], asiadf[, 1:p], seq(0, p - 1), nodes,
      adj_mat4_v2_1_1, 5, 6, sep
    )
  )
})

test_that("Rule 4 testing conditions", {
  sep <- c(1)
  amat_noresult <- matrix(c(
    0, 1, 2, 1,
    2, 0, 1, 0,
    3, 1, 0, 0,
    1, 0, 0, 0
  ), byrow = TRUE, ncol = 4)
  expect_output(
    testRule4(
      asia_dag[1:p, 1:p], asiadf[, 1:p], seq(0, p - 1), nodes,
      amat_noresult, 5, 6, sep
    )
  )
})

test_that("Rule 4 contradiction testing", {
  sep <- c(1)
  a <- which(nodes == "alpha")
  b <- which(nodes == "beta")
  # This places a tail, but based on rule 4 it should be an arrowhead
  adj_mat4_v3[a, b] <- 3
  expect_warning(
    r4_contradiction <- testRule4(
      asia_dag[1:p, 1:p], asiadf[, 1:p], seq(0, p - 1), nodes,
      adj_mat4_v3, 5, 6, sep
    )
  )
  # Tail should have been converted to an arrowhead
  expect_equal(adj_mat4_v3[a, b], 2)
})

test_that("Rule 8 is correct.", {
  adj_mat8 <- matrix(c(
    0, 1, 0, 0, 0, 2, 0,
    3, 0, 2, 0, 0, 2, 0,
    0, 3, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 2, 2, 0,
    0, 0, 0, 3, 0, 2, 1,
    1, 3, 0, 1, 3, 0, 0,
    0, 1, 0, 0, 0, 0, 0
  ), nrow = 7, byrow = TRUE)
  adj_mat8orig <- deep_copy(adj_mat8)
  p <- ncol(adj_mat8)

  expect_output(
    adj_mat8 <- testRule8(
      asia_dag[1:p, 1:p], asiadf[, 1:p], seq(0, p - 1), node_names, adj_mat8
    )
  )

  expect_equal(adj_mat8[6, 1], 3)
  expect_equal(adj_mat8[1, 6], 2)
  expect_equal(adj_mat8[6, 4], 3)
  expect_equal(adj_mat8[4, 6], 2)
  for (i in seq_len(nrow(adj_mat8))) {
    for (j in seq_len(ncol(adj_mat8))) {
      cond1 <- (i == 6 & j == 1) | (i == 6 & j == 4)
      if (!cond1) {
        expect_equal(adj_mat8[i, j], adj_mat8orig[i, j])
      }
    }
  }
})

test_that("Rule 9 is correct.", {
  adj_mat9 <- matrix(c(
    0, 1, 0, 0, 0, 0, 0, 3,
    1, 0, 1, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 1, 0, 1, 2, 0,
    0, 0, 0, 0, 1, 0, 0, 2,
    0, 0, 1, 0, 1, 0, 0, 0,
    2, 0, 0, 0, 0, 1, 0, 0
  ), nrow = 8, byrow = TRUE)
  adj_mat9orig <- deep_copy(adj_mat9)
  p <- ncol(adj_mat9)

  expect_output(
    adj_mat9 <- testRule9(
      asia_dag[1:p, 1:p], asiadf[, 1:p], seq(0, p - 1), node_names, adj_mat9
    )
  )

  expect_equal(adj_mat9[5, 7], 2)
  expect_equal(adj_mat9[7, 5], 3)

  for (i in seq_len(nrow(adj_mat9))) {
    for (j in seq_len(ncol(adj_mat9))) {
      cond1 <- i == 7 & j == 5
      if (!cond1) {
        expect_equal(adj_mat9[i, j], adj_mat9orig[i, j])
      }
    }
  }
})

test_that("Rule 9 (completeness of upd function)", {
  adj_mat9_upd <- matrix(c(
    0, 0, 0, 3, 0, 2, 0, 0,
    0, 0, 0, 1, 0, 0, 2, 0,
    0, 0, 0, 1, 2, 1, 0, 0,
    1, 2, 1, 0, 0, 2, 0, 0,
    0, 0, 1, 0, 0, 3, 0, 1,
    3, 0, 1, 3, 2, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 2, 0, 1, 0
  ), byrow = TRUE, ncol = 8)
  p <- nrow(adj_mat9_upd)
  expect_output(
    testRule9(
      asia_dag[1:p, 1:p], asiadf[, 1:p], 0, node_names[1:p], adj_mat9_upd
    )
  )
})

test_that("Rule 10 is correct.", {
  adj_mat10 <- matrix(c(
    0, 0, 0, 2, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 1, 0, 1,
    0, 0, 0, 0, 2, 0, 2, 1,
    1, 0, 0, 0, 2, 0, 0, 0,
    0, 0, 1, 3, 0, 3, 0, 0,
    0, 1, 0, 0, 2, 0, 0, 0,
    1, 0, 3, 0, 0, 0, 0, 0,
    0, 2, 1, 0, 0, 0, 0, 0
  ), nrow = 8, byrow = TRUE)
  adj_mat10orig <- deep_copy(adj_mat10)
  p <- ncol(adj_mat10)
  expect_output(
    adj_mat10 <- testRule10(
      asia_dag[1:p, 1:p], asiadf[, 1:p], seq(0, p - 1), node_names, adj_mat10
    )
  )

  expect_equal(adj_mat10[3, 5], 2)
  expect_equal(adj_mat10[5, 3], 3)
  for (i in seq_len(nrow(adj_mat10))) {
    for (j in seq_len(ncol(adj_mat10))) {
      cond1 <- (i == 5 & j == 3)
      if (!cond1) {
        expect_equal(adj_mat10[i, j], adj_mat10orig[i, j])
      }
    }
  }
})

test_that("Rule 10 (simple)", {
  adj_mat10_simple <- matrix(c(
    0, 2, 2, 1,
    1, 0, 2, 0,
    1, 3, 0, 3,
    3, 0, 2, 0
  ), nrow = 4, byrow = TRUE)
  adj_mat10_s_orig <- deep_copy(adj_mat10_simple)
  p <- ncol(adj_mat10_simple)
  expect_output(
    res <- testRule10(
      asia_dag[1:p, 1:p], asiadf[, 1:p], 0, node_names[1:p], adj_mat10_simple
    )
  )
  expect_equal(res[3, 1], 3)
  expect_equal(adj_mat10_s_orig[3, 1], 1)

  for (i in seq_len(nrow(adj_mat10_simple))) {
    for (j in seq_len(ncol(adj_mat10_simple))) {
      cond1 <- (i == 3 & j == 1)
      if (!cond1) {
        expect_equal(adj_mat10_simple[i, j], adj_mat10_s_orig[i, j])
      }
    }
  }
})

test_that("Rule 10 (capture simple unprotected pd paths)", {
  adj_mat10_upd <- matrix(c(
    0, 0, 2, 0, 1, 2,
    0, 0, 2, 0, 0, 3,
    1, 3, 0, 3, 0, 0,
    0, 0, 2, 0, 1, 0,
    1, 0, 0, 2, 0, 0,
    3, 2, 0, 0, 0, 0
  ), byrow = TRUE, nrow = 6)
  p <- nrow(adj_mat10_upd)

  expect_output(
    testRule10(
      asia_dag[1:p, 1:p], asiadf[, 1:p], 0, node_names[1:p], adj_mat10_upd
    )
  )
})

test_that("All Rules is correct.", {
  adj_mat10 <- matrix(c(
    0, 0, 0, 2, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 1, 0, 1,
    0, 0, 0, 0, 2, 0, 2, 1,
    1, 0, 0, 0, 2, 0, 0, 0,
    0, 0, 1, 3, 0, 3, 0, 0,
    0, 1, 0, 0, 2, 0, 0, 0,
    1, 0, 3, 0, 0, 0, 0, 0,
    0, 2, 1, 0, 0, 0, 0, 0
  ), nrow = 8, byrow = TRUE)
  p <- ncol(adj_mat10)
  local_node_names <- paste0("V", seq(0, p - 1))
  expect_output(
    adj_mat10 <- testAllRules(
      asia_dag[1:p, 1:p], asiadf[, 1:p], 0, local_node_names, adj_mat10
    )
  )

  expect_equal(adj_mat10[7, 1], 2)
  expect_equal(adj_mat10[1, 7], 3)
  expect_equal(adj_mat10[1, 4], 2)
  expect_equal(adj_mat10[4, 1], 3)
  expect_equal(adj_mat10[2, 6], 2)
  expect_equal(adj_mat10[6, 2], 3)
  expect_equal(adj_mat10[3, 5], 2)
  expect_equal(adj_mat10[5, 3], 3)
})

test_that("Test conversion of Mixed Graph", {
  dag_amat <- matrix(0, nrow = 12, ncol = 12)
  dag_amat[1, 2] <- dag_amat[1, 5] <- 1
  dag_amat[2, 3] <- 1
  dag_amat[3, 6] <- 1
  dag_amat[4, 3] <- 1
  dag_amat[5, 7] <- 1
  dag_amat[6, 12] <- 1
  dag_amat[7, 8] <- 1
  dag_amat[10, 8] <- dag_amat[10, 9] <- 1
  dag_amat[11, 9] <- 1
  dag_amat[12, 11] <- dag_amat[12, 10] <- 1

  prelim_mixed_graph <- matrix(c(
    0, 0, 2, 2, 0, 0, 0, 0,
    0, 0, 2, 0, 0, 0, 0, 0,
    3, 1, 0, 2, 0, 0, 1, 2,
    2, 0, 1, 0, 1, 0, 0, 2,
    0, 0, 0, 1, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 1, 3,
    0, 0, 3, 0, 0, 2, 0, 0,
    0, 0, 3, 2, 2, 1, 0, 0
  ), byrow = TRUE, nrow = 8)

  result <- testConvertMixed(
    dag_amat, c(1, 9), paste0("V", 1:12), prelim_mixed_graph,
    c(2, 4, 3, 7, 8, 9, 11, 10) - 1
  )
  # This is a real ground-truth value pulled from CML's actual (unmodified)
  # source, not the upstream test's own hardcoded snapshot -- that snapshot
  # is stale relative to CML's current code (a pre-existing upstream
  # inconsistency, not something introduced by this port; see the similar
  # note in test-cml.R).
  expect_equal(result, matrix(c(
    0, 0, 1, 2, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 2, 0, 0, 4, 2,
    2, 0, 4, 0, 1, 0, 0, 2,
    0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 3, 0, 0, 1, 0, 0,
    0, 0, 3, 2, 1, 1, 0, 0
  ), nrow = 8, byrow = TRUE))
})

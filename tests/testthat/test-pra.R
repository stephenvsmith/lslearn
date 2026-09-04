true_amat <- matrix(c(
  0, 0, 1, 0, 0, 0,
  0, 0, 0, 0, 1, 0,
  0, 0, 0, 1, 0, 0,
  0, 0, 0, 0, 1, 0,
  0, 0, 0, 0, 0, 1,
  0, 0, 0, 0, 1, 0
), nrow = 6, byrow = TRUE)

est_amat <- matrix(c(
  0, 0, 0, 0, 0, 0,
  0, 0, 1, 0, 1, 0,
  1, 0, 0, 0, 0, 0,
  0, 0, 1, 0, 1, 0,
  0, 0, 0, 1, 0, 0,
  0, 0, 0, 0, 1, 0
), nrow = 6, byrow = TRUE)

# Nodes 2 and 4 are the targets
# TP: 1 to 4
# FP: 3 to 2, 1 to 2, and 5 to 4
# FN: 0 to 2
# Potential: 3 to 4
test_that("Testing Parent Recovery Accuracy Metrics", {
  expect_output(
    pra <- parentRecoveryAccuracy(est_amat, true_amat, c(2, 4), TRUE)
  )
  expect_equal(
    pra, list(missing = 1, added = 3, correct = 1, potential = 1)
  )

  expect_output(
    all_metrics <- allMetrics(
      est_amat, true_amat, c(2, 4), true_amat, seq(0, ncol(est_amat) - 1),
      algo = "cml", verbose = TRUE
    )
  )
  expect_equal(all_metrics, data.frame(
    cml__skel_fp = 1, cml__skel_fn = 0, cml__skel_tp = 5,
    cml__v_fn = 1, cml__v_fp = 2, cml__v_tp = 0,
    cml_pra_fn = 1, cml_pra_fp = 3, cml_pra_tp = 1, cml_pra_potential = 1,
    cml_ancestors_correct = 0, cml_ancestors_incorrect = 0,
    cml_ancestors_total = 0, cml_overall_f1 = 1 / 3
  ))
})

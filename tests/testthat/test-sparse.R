test_that("sparse_auc", {
  d_sparse <-
    data.frame(
      id = c(1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 4L, 5L, 6L, 4L, 5L, 6L, 7L, 8L, 9L, 7L, 8L, 9L),
      conc = c(0, 0, 0,  1.75, 2.2, 1.58, 4.63, 2.99, 1.52, 3.03, 1.98, 2.22, 3.34, 1.3, 1.22, 3.54, 2.84, 2.55, 0.3, 0.0421, 0.231),
      time = c(0, 0, 0, 1, 1, 1, 6, 6, 6, 2, 2, 2, 10, 10, 10, 4, 4, 4, 24, 24, 24),
      dose = c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100)
    )

  # calculated using the PK library with
  # v_batch <- PK::auc(data=d_sparse, method="t", design="batch")
  # v_serial <- PK::auc(data=d_sparse, method="t", design="ssd")
  auclast <- 39.4689 # using linear trapezoidal rule
  auclast_se_batch <- 7.30997787038754 # for a batch design (with multiple measures from the same animal taken into account)
  auclast_df_batch <- 2.74598236184576
  auclast_se_serial <- 6.86584835083522 # for a serial design (with multiple measures from the same animal not taken into account)
  auclast_df_serial <- 2.82631153092225

  expect_warning(
    sparse_batch <- pk.calc.sparse_auc(conc=d_sparse$conc, time=d_sparse$time, subject=d_sparse$id),
    regexp="Cannot yet calculate sparse degrees of freedom for multiple samples per subject"
  )
  expect_equal(sparse_batch$sparse_auc, auclast)
  expect_equal(sparse_batch$sparse_auc_se, auclast_se_batch)
  expect_equal(sparse_batch$sparse_auc_df, NA_real_)

  sparse_serial <- pk.calc.sparse_auc(conc=d_sparse$conc, time=d_sparse$time, subject=seq_len(nrow(d_sparse)))
  expect_equal(sparse_serial$sparse_auc, auclast)
  expect_equal(as.numeric(sparse_serial$sparse_auc_se), auclast_se_serial)
  expect_equal(sparse_serial$sparse_auc_df, auclast_df_serial)
})

test_that("sparse_auclast expected errors", {
  expect_error(
    pk.calc.sparse_auclast(auc.type = "foo"),
    class = "pknca_sparse_auclast_change_auclast"
  )
})

test_that("sparse_auc_df and sparse_auc_se are in the parameter list (#292)", {
  expect_true(
    all(c("sparse_auc_df", "sparse_auc_se") %in% names(get.interval.cols()))
  )
})

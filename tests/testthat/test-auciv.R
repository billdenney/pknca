test_that("pk.calc.auciv", {
  expect_equal(
    pk.calc.auciv(conc = 1:2, time = 1:2),
    structure(NA_real_, exclude="No time 0 in data")
  )
  expect_equal(
    # No check is done to confirm that the auc argument matches the data
    pk.calc.auciv(conc = 0:5, time = 0:5, c0 = 1, auc = 2.75),
    2.75 + 1 - 0.5
  )
  expect_equal(
    # No verifications are made on the data
    pk.calc.auciv(conc = 0:5, time = 0:5, c0 = 1, auc = 2.75, check=FALSE),
    2.75 + 1 - 0.5
  )
})

test_that("pk.calc.auciv_pbext", {
  expect_equal(
    pk.calc.auciv_pbext(auc = 1, auciv = 2.1),
    100 * (1 - 1/2.1)
  )
})

test_that("NA data are removed from concentrations for calculation of AUCiv (#353)", {
  d_iv_353alt <- data.frame(conc = c(NA, 4, 2, 1, 0.45), time = c(0, 5, 15, 30, 60))
  d_intervals <- data.frame(start = 0, end = Inf, aucivinf.obs = TRUE)
  o_conc_353alt <- PKNCAconc(data = d_iv_353alt, conc~time)
  o_dose <- PKNCAdose(data = data.frame(time = 0), ~time)
  o_data_353alt <- PKNCAdata(o_conc_353alt, o_dose, intervals = d_intervals)
  # The same warning is expected three times
  expect_warning(expect_warning(expect_warning(
    o_nca <- pk.nca(o_data_353alt),
    regexp = "Requesting an AUC range starting (0) before the first measurement (5) is not allowed",
    fixed = TRUE),
    regexp = "Requesting an AUC range starting (0) before the first measurement (5) is not allowed",
    fixed = TRUE),
    regexp = "Requesting an AUC range starting (0) before the first measurement (5) is not allowed",
    fixed = TRUE
  )
  expect_s3_class(o_nca, "PKNCAresults")
})

test_that("missing dose information does not cause NA time (#353)", {
  d_iv_nodose <- data.frame(conc = c(4, 2, 1, 0.45), time = c(5, 15, 30, 60))
  d_intervals <- data.frame(start = 0, end = Inf, aucivinf.obs = TRUE)
  o_conc_nodose <- PKNCAconc(data = d_iv_nodose, conc~time)
  o_data_nodose <- PKNCAdata(o_conc_nodose, intervals = d_intervals, impute = "start_conc0")
  expect_warning(
    o_nca <- pk.nca(o_data_nodose),
    regexp = "time.dose is NA"
  )
  expect_s3_class(o_nca, "PKNCAresults")
})

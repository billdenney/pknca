test_that("update.PKNCAresults", {
  d_conc <- generate.conc(2, 1, c(0, 2, 6, 12, 24))
  d_dose <- generate.dose(d_conc)
  o_conc <- PKNCAconc(d_conc, formula=conc~time|treatment+ID)
  o_dose <- PKNCAdose(d_dose, formula=dose~time|treatment+ID)
  o_data <- PKNCAdata(o_conc, o_dose, intervals = data.frame(start = 0, end = Inf, half.life = TRUE))
  o_result <- pk.nca(o_data)

  expect_message(
    update(o_result, o_data),
    regexp = "No changes detected in data"
  )
  o_data_changed <-
    PKNCAdata(
      o_conc, o_dose,
      intervals = data.frame(start = 0, end = Inf, half.life = TRUE),
      units = pknca_units_table(concu = "foo", doseu = "bar", timeu = "baz")
    )

  expect_warning(
    update(o_result, o_data_changed),
    regexp = "changes detected in data other than source concentration or dose data"
  )

  # Change concentration ----
  d_conc_changed <- d_conc
  d_conc_changed$conc[2] <- 1
  d_dose_changed <- d_dose
  d_dose_changed$dose[2] <- 2
  o_conc_changed <- PKNCAconc(d_conc_changed, formula=conc~time|treatment+ID)
  o_dose_changed <- PKNCAdose(d_dose_changed, formula=dose~time|treatment+ID)

  o_data_chconc <- PKNCAdata(o_conc_changed, o_dose, intervals = data.frame(start = 0, end = Inf, half.life = TRUE))
  o_data_chdose <- PKNCAdata(o_conc, o_dose_changed, intervals = data.frame(start = 0, end = Inf, half.life = TRUE))

  o_conc_changed_reordered <- PKNCAconc(d_conc_changed[order(-d_conc_changed$ID), ], formula=conc~time|treatment+ID)
  o_data_chconc_reordered <- PKNCAdata(o_conc_changed_reordered, o_dose, intervals = data.frame(start = 0, end = Inf, half.life = TRUE))

  o_nca_chconc_reordered <- pk.nca(o_data_chconc_reordered)
  expect_equal(
    update(o_result, o_data_chconc)$result,
    o_nca_chconc_reordered$result
  )
})

test_that("all units have been defined", {
  all_units <- pknca_units_table(concu="conc", doseu="dose", amountu="amount", timeu="time")
  expect_equal(
    setdiff(names(get.interval.cols()), c(all_units$PPTESTCD, "start", "end")),
    character()
  )
})

test_that("pknca_find_units_param expected errors", {
  expect_error(
    pknca_find_units_param(unit_type=1:2)
  )
  expect_error(
    pknca_find_units_param(unit_type=1)
  )
  expect_error(
    pknca_find_units_param(unit_type="foo"),
    regexp="No parameters found for unit_type"
  )
})

test_that("pknca_find_units_param", {
  expect_true(
    "cmax" %in% pknca_find_units_param(unit_type="conc")
  )
})

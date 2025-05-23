test_that("assert_intervaltime_single", {
  expect_equal(
    assert_intervaltime_single(interval = c(1, 2)),
    c(1, 2)
  )
  expect_equal(
    assert_intervaltime_single(start = 1, end = 2),
    c(1, 2)
  )
  expect_equal(
    assert_intervaltime_single(interval = c(1, 2), start = 1, end = 2),
    c(1, 2)
  )
  expect_error(
    assert_intervaltime_single(interval = c(1, 2), start = 1.1, end = 2),
    regexp = "`start` must be the same as the first value in the interval if both are given: 1.1!=1",
    fixed = TRUE
  )
  expect_error(
    assert_intervaltime_single(interval = c(1, 2), start = 1, end = 2.1),
    regexp = "`end` must be the same as the second value in the interval if both are given: 2.1!=2",
    fixed = TRUE
  )
})

test_that("assert_conc", {
  expect_warning(
    assert_conc(conc = -1),
    regexp="Negative concentrations found"
  )
  expect_warning(
    assert_conc(conc = c(NA, -1)),
    regexp="Negative concentrations found"
  )
  expect_warning(
    assert_conc(conc = c(NA, -1, 1)),
    regexp="Negative concentrations found"
  )
  expect_warning(
    assert_conc(conc = NA),
    regexp="All concentration data are missing"
  )
  expect_error(
    assert_conc(conc="A"),
    regexp="Assertion on 'conc' failed: Must be of type 'numeric', not 'character'."
  )
  expect_error(
    assert_conc(conc=factor("A")),
    regexp="Assertion on 'conc' failed: Must be of type 'numeric', not 'factor'."
  )
})

test_that("assert_time", {
  expect_error(
    assert_time(time=NA),
    regexp="Assertion on 'time' failed: Contains missing values (element 1).",
    fixed = TRUE
  )
  expect_error(
    assert_time(time = c(0, 0)),
    regexp="Assertion on 'time' failed: Contains duplicated values, position 2."
  )
  expect_error(
    assert_time(time = c(1, 0)),
    regexp="Assertion on 'time' failed: Must be sorted."
  )
  expect_error(
    assert_time(time="A"),
    regexp="Assertion on 'time' failed: Must be of type 'numeric', not 'character'."
  )
  expect_error(
    assert_time(time=factor("A")),
    regexp="Assertion on 'time' failed: Must be of type 'numeric', not 'factor'."
  )
})

test_that("assert_conc_time", {
  expect_error(
    assert_conc_time(conc = 1, time = 1:2),
    regexp="Assertion on 'conc' failed: Must have length 2, but has length 1."
  )
  expect_error(
    assert_conc_time(conc = 1:2, time = 2),
    regexp="Assertion on 'conc' failed: Must have length 1, but has length 2."
  )
})

test_that("assert_lambdaz", {
  expect_equal(assert_lambdaz(1), 1)
  expect_equal(assert_lambdaz(NA), NA) # NA is allowed by default
  expect_error(
    assert_lambdaz(NA, any.missing = FALSE),
    regexp = "Assertion on 'NA' failed: Contains missing values (element 1).",
    fixed = TRUE
  )
  expect_error(
    assert_lambdaz(-1),
    regexp = "Assertion on '-1' failed: Element 1 is not > 0"
  )
})

test_that("assert_PKNCAconc", {
  expect_error(
    assert_PKNCAconc("A"),
    regexp = "Must be a PKNCAconc object"
  )
})

test_that("assert_PKNCAdata", {
  expect_error(
    assert_PKNCAdata("A"),
    regexp = "Must be a PKNCAdata object"
  )
})

test_that("assert_PKNCAdose", {
  expect_error(
    assert_PKNCAdose("A"),
    regexp = "Must be a PKNCAdose object"
  )
})

test_that("assert_PKNCAresults", {
  expect_error(
    assert_PKNCAresults("A"),
    regexp = "Must be a PKNCAresults object"
  )
})

test_that("element_find", {
  values5<-c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)
  values10<-c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)
  expect_equal(
    element_find(values5),
    element_find(values10)
  )

  expect_equal(
    element_find(values5),
    "Elements 1, 2, 3, 4, 5"
  )

  expect_equal(
    element_find(values10),
    "Elements 1, 2, 3, 4, 5"
  )
})

test_that("assert_unit_col", {
  d <- data.frame(A = 1, B = "C")
  expect_equal(
    assert_unit_col(unit = "B", data = d),
    structure("B", unit_type = "column")
  )
  expect_error(
    assert_unit_col(unit = 1:2),
    regexp = "`unit` must be a single value"
  )
  expect_error(
    assert_unit_col(unit = 1),
    regexp = "`unit` must be a character string"
  )
  expect_error(
    assert_unit_col(unit = "D", data = "A"),
    regexp = "`data` must be a data.frame"
  )
  expect_error(
    assert_unit_col(unit = "D", data = d),
    regexp = "`unit` (D) must be a column name in the data",
    fixed = TRUE
  )
  expect_error(
    assert_unit_col(unit = "A", data = d),
    regexp = "`unit` (A) must contain character data",
    fixed = TRUE
  )
})

test_that("assert_unit_value", {
  expect_equal(assert_unit_value("A"), structure("A", unit_type = "value"))
  expect_null(assert_unit_value(NULL))
  expect_error(
    assert_unit_value(c("A", "B")),
    regexp = "`unit` must be a single value"
  )
  expect_error(
    assert_unit_value(1),
    regexp = "`unit` must be a character string"
  )
})

test_that("assert_unit", {
  d <- data.frame(A = 1, B = "C")

  expect_null(assert_unit(NULL, data = d))
  expect_equal(
    assert_unit(unit = "B", data = d),
    structure("B", unit_type = "column")
  )
  expect_equal(
    assert_unit(unit = "C", data = d),
    structure("C", unit_type = "value")
  )
  expect_error(
    assert_unit(unit = 1, data = d),
    regexp = "`unit` must be a character string"
  )
})

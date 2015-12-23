context("Check Conversion")

test_that("check.conversion", {
  good <- LETTERS
  expect_equal(check.conversion(good, as.character), good)
  expect_error(check.conversion(good, as.numeric),
               regexp="26 new NA value\\(s\\) created during conversion")
  good <- 1:5
  expect_equal(check.conversion(good, as.character),
               as.character(good))
})

context("Check concentration and time inputs")

test_that("check.conc.time", {
  ## Check all the invalid cases
  expect_warning(check.conc.time(conc=-1),
                 regexp="Negative concentrations found")
  expect_warning(check.conc.time(conc=c(NA, -1)),
                 regexp="Negative concentrations found")
  expect_warning(check.conc.time(conc=c(NA, -1, 1)),
                 regexp="Negative concentrations found")
  expect_error(check.conc.time(time=NA),
               regexp="Time may not be NA")
  expect_error(check.conc.time(time=c(0, 0)),
               regexp="Time must be monotonically increasing")
  expect_error(check.conc.time(time=c(1, 0)),
               regexp="Time must be monotonically increasing")
  expect_error(check.conc.time(conc=1, time=1:2),
               regexp="Conc and time must be the same length")
  expect_error(check.conc.time(conc=1:2, time=2),
               regexp="Conc and time must be the same length")  
})

context("Rounding to string values")

test_that("Rounding", {
  expect_error(roundString(1, c(2, 3)),
               regexp="digits must either be a scalar or the same length as x")
  expect_equal(roundString(11), "11")
  expect_equal(roundString(5), "5")
  expect_equal(roundString(0.05), "0")
  expect_equal(roundString(NA), "NA")
  expect_equal(roundString(NaN), "NaN")
  expect_equal(roundString(Inf), "Inf")
  expect_equal(roundString(-Inf), "-Inf")
  ## Respecting the digits
  expect_equal(roundString(0.05, 3), "0.050")
  expect_equal(roundString(123.05, 3), "123.050")
  ## Vectors work with different orders of magnitude work
  expect_equal(roundString(c(100, 0.1), 3), c("100.000", "0.100"))
  ## Vectors of digits work
  expect_equal(roundString(c(100, 0.1), c(0, 3)), c("100", "0.100"))
})

test_that("Significance", {
  expect_equal(signifString(11), "11.0000")
  expect_equal(signifString(5), "5.00000")
  expect_equal(signifString(0.05), "0.0500000")
  expect_equal(signifString(NA), "NA")
  expect_equal(signifString(NaN), "NaN")
  expect_equal(signifString(Inf), "Inf")
  expect_equal(signifString(-Inf), "-Inf")
  ## Respecting the digits
  expect_equal(signifString(0.05, 3), "0.0500")
  expect_equal(signifString(123.05, 3), "123")
  expect_equal(signifString(123456.05, 3), "123000")
  expect_equal(signifString(-123000.05, 3), "-123000")
  expect_equal(signifString(-.05, 3), "-0.0500")
  ## Exact orders of magnitude work on both sides of 0
  expect_equal(signifString(0.01, 3), "0.0100")
  expect_equal(signifString(1, 3), "1.00")
  expect_equal(signifString(100, 3), "100")
  ## Vectors work with different orders of magnitude work
  expect_equal(signifString(c(100, 0.1), 3), c("100", "0.100"))
  ## Rounding to a higher number of significant digits works correctly
  expect_equal(signifString(0.9999999, 3), "1.00")
})

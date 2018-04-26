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
  expect_warning(check.conc.time(conc=NA),
                 regexp="All concentration data is missing")
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
  expect_error(check.conc.time(conc="A"),
               regexp="Concentration data must be numeric and not a factor")
  expect_error(check.conc.time(conc=factor("A")),
               regexp="Concentration data must be numeric and not a factor")
  expect_error(check.conc.time(time="A"),
               regexp="Time data must be numeric and not a factor")
  expect_error(check.conc.time(time=factor("A")),
               regexp="Time data must be numeric and not a factor")
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
  expect_equal(roundString(c(100, 0.1), 3), c("100.000", "0.100"),
               info="Vectors work with different orders of magnitude work")
  expect_equal(roundString(c(100, 0.1), c(0, 3)), c("100", "0.100"),
               info="Vectors of digits work")
  expect_equal(roundString(c(0.1, NA), digits=3), c("0.100", "NA"),
               info="Mixed inputs (NA, NaN, Inf or numeric), NA")
  expect_equal(roundString(c(0.1, NA, NaN, Inf, -Inf), digits=3),
               c("0.100", "NA", "NaN", "Inf", "-Inf"),
               info="Mixed inputs (NA, NaN, Inf or numeric)")
  ## All zeros
  expect_equal(roundString(0, digits=3), "0.000")
  expect_equal(roundString(c(0, NA), digits=3), c("0.000", "NA"))
  # scientific notation
  expect_equal(roundString(1234567, digits=3, si_range=5), "1.234567000e6",
               info="si_range works with roundString (even if it looks odd)")
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
  expect_equal(signifString(123456.05, 3, si_range=6), "123000")
  expect_equal(signifString(123456.05, 3, si_range=5), "1.23e5")
  expect_equal(signifString(-123000.05, 3, si_range=5), "-1.23e5")
  expect_equal(signifString(999999, 3, si_range=6), "1.00e6",
               info="Rounding around the edge of the si_range works correctly (going up)")
  expect_equal(signifString(999999, 7, si_range=6), "999999.0",
               info="Rounding around the edge of the si_range works correctly (going staying the same)")
  expect_equal(signifString(-.05, 3), "-0.0500")
  ## Exact orders of magnitude work on both sides of 0
  expect_equal(signifString(0.01, 3), "0.0100")
  expect_equal(signifString(1, 3), "1.00")
  expect_equal(signifString(100, 3), "100")
  ## Vectors work with different orders of magnitude work
  expect_equal(signifString(c(100, 0.1), 3), c("100", "0.100"))
  ## Rounding to a higher number of significant digits works correctly
  expect_equal(signifString(0.9999999, 3), "1.00")
  ## Mixed inputs (NA, NaN, Inf or numeric)
  expect_equal(signifString(NA), "NA")
  expect_equal(signifString(c(0.1, NA), digits=3), c("0.100", "NA"))
  expect_equal(signifString(c(0.1, NA, NaN, Inf, -Inf), digits=3),
               c("0.100", "NA", "NaN", "Inf", "-Inf"))
  ## All zeros
  expect_equal(signifString(0, digits=3), "0.000")
  expect_equal(signifString(c(0, NA), digits=3), c("0.000", "NA"))
  
  # Data Frames
  expect_equal(signifString(data.frame(A=c(0, 1.111111),
                                       B=factor(LETTERS[1:2]),
                                       C=LETTERS[1:2],
                                       stringsAsFactors=FALSE),
                            digits=3),
               data.frame(A=c("0.000", "1.11"),
                          B=factor(LETTERS[1:2]),
                          C=LETTERS[1:2],
                          stringsAsFactors=FALSE),
               check.attributes=FALSE,
               info="Data frame significance is calculated correctly")
})

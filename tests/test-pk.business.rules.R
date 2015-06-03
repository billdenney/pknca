library(testthat)

test_that("geomean", {
  ## Test normal, nonzero scalars
  expect_equal(geomean(5), 5)
  expect_equal(geomean(c(5, 5)), 5)
  expect_equal(geomean(c(10, 1000)), 100)
  ## Test zeros
  expect_equal(geomean(0), 0)
  expect_equal(geomean(c(0, 1)), 0)
  ## Test NAs
  expect_equal(geomean(NA), as.numeric(NA))
  expect_equal(geomean(c(NA, 0)), as.numeric(NA))
  expect_equal(geomean(c(NA, 5), na.rm=TRUE), 5)
  expect_equal(geomean(c(NA, NA), na.rm=TRUE), NaN)
})

test_that("business.mean", {
  PKNCA.options(default=TRUE)
  PKNCA.options(max.missing=0.5)
  ## Test a normal mean of a scalar and vector
  expect_equal(business.mean(1), 1)
  expect_equal(business.mean(c(1, 2)), 1.5)
  ## Ensure that at the max.missing fraction a value is reported and
  ## above that, NA is returned.
  expect_equal(business.mean(c(NA, NA, 1, 2)), 1.5)
  expect_equal(business.mean(c(NA, NA, NA, 2)), NA)
  ## Ensure that it uses the current value of max.missing
  PKNCA.options(max.missing=0.3)
  expect_equal(business.mean(c(NA, 1, 2, 3)), 2)
  expect_equal(business.mean(c(NA, NA, 1, 2)), NA)
})

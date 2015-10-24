context("Business rule functions")

test_that("geomean", {
  ## Test normal, nonzero scalars
  expect_equal(geomean(5), 5)
  expect_equal(geomean(c(5, 5)), 5)
  expect_equal(geomean(c(10, 1000)), 100)
  ## Test zeros
  expect_equal(geomean(0), 0)
  expect_equal(geomean(c(0, 1)), 0)
  ## Test negative numbers
  expect_equal(geomean(c(-1, 2)), -geomean(c(1, 2)))
  ## Test NAs
  expect_equal(geomean(NA), as.numeric(NA))
  expect_equal(geomean(c(NA, 0)), as.numeric(NA))
  expect_equal(geomean(c(NA, 5), na.rm=TRUE), 5)
  expect_equal(geomean(c(NA, NA), na.rm=TRUE), NaN)
})

test_that("geosd", {
  expect_equal(geosd(c(1, 2)), exp(sd(log(c(1, 2)))))
  expect_equal(geosd(c(NA, 1, 2)), as.numeric(NA))
})

test_that("geocv", {
  expect_equal(geocv(c(1, 2)),
               sqrt(exp(sd(log(c(1, 2)))^2)-1)*100)
  expect_equal(geocv(c(NA, 1, 2)), as.numeric(NA))
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

test_that("pk.business", {
  PKNCA.options(default=TRUE)
  PKNCA.options(max.missing=0.5)
  b.mean <- pk.business(mean)
  expect_equal(b.mean(c(1, 2)), 1.5)
  ## Right at the border, it still reports
  expect_equal(b.mean(c(1, NA)), 1)
  ## When too much data is missing, NA is returned
  expect_equal(b.mean(c(1, NA, NA)), NA)
  ## It respects zero.missing
  b.mean.2 <- pk.business(mean, zero.missing=TRUE)
  expect_equal(b.mean(c(0, 1)), 0.5)
  expect_equal(b.mean.2(c(0, 1)), 1)
})

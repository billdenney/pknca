library(testthat)

test_that("check.conversion", {
  good <- LETTERS
  expect_equal(check.conversion(good, as.character), good)
  expect_error(check.conversion(good, as.numeric),
               regexp="26 new NA value\\(s\\) created during conversion")
  good <- 1:5
  expect_equal(check.conversion(good, as.character),
               as.character(good))
})

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

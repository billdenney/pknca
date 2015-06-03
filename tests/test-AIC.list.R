library(testthat)

test_that("AIC.list", {
  
})

test_that("get.first.model", {
  ## If only given NA, return NA.
  expect_equal(get.first.model(NA), NA)
  ## If given a list of NAs, return NA
  expect_equal(get.first.model(list(NA)), NA)
  expect_equal(get.first.model(list(NA, NA)), NA)
  ## Named lists are treated the same way
  expect_equal(get.first.model(list(A=NA, B=NA)), NA)
  ## If given a recursive list of all NAs, return NA
  expect_equal(get.first.model(list(list(NA), list(NA))), NA)
  expect_equal(get.first.model(list(list(NA, NA))), NA)
  expect_equal(get.first.model(list(list(list(NA, NA)))), NA)

  ## If given something usable, return that
  expect_equal(get.first.model(1), 1)
  expect_equal(get.first.model("B"), "B")

  ## Return the first usable object found
  expect_equal(get.first.model(list("A", "B")), "A")
  ## If the first usable object is deep within a list, return that
  expect_equal(get.first.model(list(list(list("A")), "B")), "A")
  ## If there are NAs before the first usable object, return the first
  ## thing that is not NA.
  expect_equal(get.first.model(list(NA, "A")), "A")
  expect_equal(get.first.model(list(NA, list("A"))), "A")
  expect_equal(get.first.model(list(list(NA, "A"), "B")), "A")
})

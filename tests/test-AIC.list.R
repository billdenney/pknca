library(testthat)

test_that("AIC.list", {
  tmpdat <- data.frame(x=c(1, 2, 3), y=c(1, 2.02, 3), z=c(0, 0.01, 0))
  mod1 <- glm(y~x, data=tmpdat)
  mod2 <- glm(y~x+z, data=tmpdat)
  AIC.list(list(mod1, mod2))
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

  ## Ensure that we do not recurse into items that are not of the list
  ## class (but are lists internally)
  tmpdat <- data.frame(x=c(1, 2, 3), y=c(1, 2.02, 3), z=c(0, 0.01, 0))
  mod1 <- glm(y~x, data=tmpdat)
  expect_equal(get.first.model(list(list(NA, mod1), "B")), mod1)
})

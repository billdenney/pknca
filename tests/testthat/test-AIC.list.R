context("AIC.list testing")
library(PKNCA)

test_that("AIC.list", {
  tmpdat <- data.frame(x=c(1, 2, 3), y=c(1, 2.02, 3), z=c(0, 0.01, 0))
  mod1 <- glm(y~x, data=tmpdat)
  mod2 <- glm(y~x+z, data=tmpdat)
  d1 <- list(list(mod1), mod2)
  result1 <- data.frame(AIC=c(as.numeric(AIC(mod1)), as.numeric(AIC(mod2))),
                        df=c(3, 4),
                        indentation=c(1, 0),
                        isBest=c("", "Best Model"),
                        stringsAsFactors=FALSE)
  expect_equal(AIC.list(d1), result1)

  ## Check that simple names apply correctly
  d2 <- list(list("A"=mod1), "B"=mod2)
  result2 <- data.frame(AIC=c(as.numeric(AIC(mod1)), as.numeric(AIC(mod2))),
                        df=c(3, 4),
                        indentation=c(1, 0),
                        isBest=c("", "Best Model"),
                        stringsAsFactors=FALSE,
                        row.names=c("A", "B"))
  expect_equal(AIC.list(d2), result2)

  d3 <- list("A"=mod1, "B"=mod2)
  result3 <- data.frame(AIC=c(as.numeric(AIC(mod1)), as.numeric(AIC(mod2))),
                        df=c(3, 4),
                        indentation=c(0, 0),
                        isBest=c("", "Best Model"),
                        stringsAsFactors=FALSE,
                        row.names=c("A", "B"))
  expect_equal(AIC.list(d3), result3)
  
  d4 <- list("C"=list("A"=mod1), "B"=mod2)
  result4 <- data.frame(AIC=c(as.numeric(AIC(mod1)), as.numeric(AIC(mod2))),
                        df=c(3, 4),
                        indentation=c(1, 0),
                        isBest=c("", "Best Model"),
                        stringsAsFactors=FALSE,
                        row.names=c("C A", "B"))
  expect_equal(AIC.list(d4), result4)

  d5 <- list("C"=list("A"=mod1), "B"=mod2, C=NA)
  result5 <- data.frame(AIC=c(as.numeric(AIC(mod1)), as.numeric(AIC(mod2)), NA),
                        df=c(3, 4, NA),
                        indentation=c(1, 0, 0),
                        isBest=c("", "Best Model", ""),
                        stringsAsFactors=FALSE,
                        row.names=c("C A", "B", "C"))
  expect_equal(AIC.list(d5), result5)

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

test_that("get.best.model", {
  tmpdat <- data.frame(x=c(1, 2, 3), y=c(1, 2.02, 3), z=c(0, 0.01, 0))
  mod1 <- glm(y~x, data=tmpdat)
  mod2 <- glm(y~x+z, data=tmpdat)
  d1 <- list(list(mod1), mod2)
  expect_equal(get.best.model(d1), mod2)
})

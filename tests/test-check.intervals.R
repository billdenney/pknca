library(testthat)

test_that("check.interval.specification", {
  ## Test that all valid combinations pass through without error or
  ## change
  good <-
    merge(merge(data.frame(start=0, end=1),
                data.frame(auc.type=c("AUClast", "AUCinf", "AUCall"))),
          data.frame(half.life=c(FALSE, TRUE)))
  expect_equal(check.interval.specification(good), good)

  ## A non-data.frame gets coerced but is otherwise unharmed
  good.matrix <- as.matrix(good)
  good.matrix[,'half.life'] <- rep(c("F", "T"), each=3)
  expect_warning(check.interval.specification(good.matrix),
                 regexp="AUC specification must be a data.frame")
  ## Missing the start column is an error
  expect_error(check.interval.specification(good[,setdiff(names(good), "start")]),
               regexp="AUC specification must have columns for start, end, auc.type, half.life Column\\(s\\) missing: start")

  ## half.life not as a logical is a warning.  Not coercable is an
  ## error
  for (n in c("half.life")) {
    almost <- good
    almost[,n] <- 0
    expect_warning(
      check.interval.specification(almost),
      regexp=sprintf("AUC specification for '%s' must be a logical vector, attempting conversion", n),
      info=n)
    almost[,n] <- "B"
    expect_error(check.interval.specification(almost),
                 regexp="6 new NA value\\(s\\) created during conversion",
                 info=n)
  }

  ## start and end not as a number is a warning.  Not coercable is an
  ## error
  for (n in c("start", "end")) {
    almost <- good
    almost[,n] <- "0.5"
    expect_warning(
      check.interval.specification(almost),
      regexp=sprintf("AUC specification for '%s' time must be numeric( or NA)?, attempting conversion", n),
      info=n)
    almost[,n] <- "B"
    expect_error(check.interval.specification(almost),
                 regexp="6 new NA value\\(s\\) created during conversion",
                 info=n)
  }

  ## Start can't be NA
  bad <- good
  bad$start <- NA
  expect_error(check.interval.specification(bad),
               regexp="AUC specification may not have NA for the starting time")

  ## end can't be NA when either last or all is specified
  for (n in c("AUClast", "AUCinf", "AUCall")) {
    bad <- data.frame(start=1, end=NA, auc.type=n, half.life=TRUE)
    expect_error(
      check.interval.specification(bad),
      regexp="AUC specification may not have NA for the end time and request an auc.type")
  }

  ## Can't have end=NA and half.life=FALSE
  bad <- data.frame(start=1, end=NA, auc.type=NA, half.life=FALSE)
  expect_error(
    check.interval.specification(bad),
    regexp="AUC specification may not have NA for the end time and not request half.life")

  ## Start after end is an error
  bad <- data.frame(start=1, end=0, auc.type="AUCinf", half.life=FALSE)
  expect_error(check.interval.specification(bad),
               regexp="AUC specification end must be after the start when end is given")
})

test_that("make.logical", {
  ## Simple identities
  expect_equal(make.logical(TRUE), TRUE)
  expect_equal(make.logical("TRUE"), TRUE)
  expect_equal(make.logical(factor("TRUE")), TRUE)
  expect_equal(make.logical("T"), TRUE)
  expect_equal(make.logical("YES"), TRUE)
  expect_equal(make.logical("Y"), TRUE)
  expect_equal(make.logical(1), TRUE)
  expect_equal(make.logical(Inf), TRUE)
  expect_equal(make.logical(FALSE), FALSE)
  expect_equal(make.logical("FALSE"), FALSE)
  expect_equal(make.logical(factor("FALSE")), FALSE)
  expect_equal(make.logical("F"), FALSE)
  expect_equal(make.logical("NO"), FALSE)
  expect_equal(make.logical("N"), FALSE)
  expect_equal(make.logical(0), FALSE)
  ## NA conversion works in all types
  expect_equal(make.logical(c(0, NA, 1), na.value=NA),
               c(FALSE, NA, TRUE))
  expect_equal(make.logical(c(0, NA, 1), na.value=FALSE),
               c(FALSE, FALSE, TRUE))
  expect_equal(make.logical(c(0, NA, 1), na.value=TRUE),
               c(FALSE, TRUE, TRUE))
  ## Default na.value is FALSE
  expect_equal(make.logical(c(0, NA, 1)),
               c(FALSE, FALSE, TRUE))

})

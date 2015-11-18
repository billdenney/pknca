context("Missing and BLQ data cleaners")

test_that("clean.conc.na", {
  ## Check that it gives errors for invalid conc.na values
  expect_error(clean.conc.na(conc=1, time=1, conc.na="foo"),
               regexp="conc.na must either be a finite number or the text 'drop'")

  ## It drops NA values if requested (even if they are the only value)
  expect_warning(v1 <-
    clean.conc.na(conc=as.numeric(NA), time=1, conc.na="drop"))
  expect_equal(v1,
               data.frame(conc=numeric(), time=numeric()))
  expect_warning(v2 <- 
    clean.conc.na(conc=as.numeric(c(NA, NA)), time=1:2,
                  conc.na="drop"))
  expect_equal(v2,
               data.frame(conc=numeric(), time=numeric()))
  expect_equal(clean.conc.na(conc=c(1, NA), time=1:2, conc.na="drop"),
               data.frame(conc=1, time=1))
  expect_equal(clean.conc.na(conc=c(1, NA, 2), time=1:3, conc.na="drop"),
               data.frame(conc=1:2, time=c(1, 3)),
               check.attributes=FALSE)

  ## It also works with a number as the conc.na value
  expect_warning(v3 <-
    clean.conc.na(conc=as.numeric(NA), time=1, conc.na=5))
  expect_equal(v3,
               data.frame(conc=5, time=1))
  expect_warning(v4 <-
    clean.conc.na(conc=c(NA, NA), time=1:2, conc.na=5))
  expect_equal(v4,
               data.frame(conc=c(5, 5), time=1:2))
  expect_equal(clean.conc.na(conc=c(1, NA), time=1:2, conc.na=5),
               data.frame(conc=c(1, 5), time=1:2))
  expect_equal(clean.conc.na(conc=c(1, NA, 2), time=1:3, conc.na=5),
               data.frame(conc=c(1, 5, 2), time=1:3))

  ## It doesn't touch BLQ
  expect_equal(clean.conc.na(conc=0, time=1, conc.na="drop"),
               data.frame(conc=0, time=1))
  expect_equal(clean.conc.na(conc=c(0, 0), time=1:2, conc.na="drop"),
               data.frame(conc=c(0, 0), time=1:2))
  expect_equal(clean.conc.na(conc=c(1, 0), time=1:2, conc.na="drop"),
               data.frame(conc=c(1, 0), time=1:2))
  expect_equal(clean.conc.na(conc=c(1, 0, 2), time=1:3, conc.na="drop"),
               data.frame(conc=c(1, 0, 2), time=1:3))

  ## Additional arguments are returned into the data frame and
  ## filtered if NA conc goes with them.
  expect_equal(clean.conc.na(conc=1:3,
                             time=1:3,
                             extra=c("a", "b", "c"),
                             conc.na="drop"),
               data.frame(conc=1:3,
                          time=1:3,
                          extra=c("a", "b", "c")))
  expect_equal(clean.conc.na(conc=c(1, NA, 2),
                             time=1:3,
                             extra=c("a", "b", "c"),
                             conc.na="drop"),
               data.frame(conc=1:2,
                          time=c(1, 3),
                          extra=c("a", "c")),
               check.attributes=FALSE)
})

test_that("clean.conc.blq", {
  ## If there are no BLQ values or NA values, it does nothing
  d.test <- data.frame(conc=1, time=1)
  expect_equal(clean.conc.blq(d.test$conc, d.test$time,
                                conc.blq="drop", conc.na="drop"),
               d.test)

  ## If there are NA values, it runs conc.na on those and returns that
  d.test <- data.frame(conc=c(1, NA), time=1:2)
  expect_equal(clean.conc.blq(d.test$conc, d.test$time,
                                conc.blq="drop", conc.na="drop"),
               clean.conc.na(d.test$conc, d.test$time, conc.na="drop"))

  ## If there are BLQ values at the beginning, drops those if given a
  ## simple drop
  d.test <- data.frame(conc=c(0, 1), time=1:2)
  expect_equal(clean.conc.blq(d.test$conc, d.test$time,
                                conc.blq="drop", conc.na="drop"),
               d.test[2,])
  ## It leaves beginning BLQs if instructed with a list.
  d.test <- data.frame(conc=c(0, 1), time=1:2)
  expect_equal(clean.conc.blq(d.test$conc, d.test$time,
                              conc.blq=list(
                                first="keep",
                                middle="drop",
                                last="drop"),
                              conc.na="drop"),
               d.test)

  ## Errors in how to handle first/middle/last rules are caught.
  d.test <- data.frame(conc=c(0, 1), time=1:2)
  expect_error(clean.conc.blq(d.test$conc, d.test$time,
                              conc.blq=list(
                                first="foo",
                                middle="drop",
                                last="drop"),
                              conc.na="drop"),
               regexp="conc.blq must either be a finite number or the text 'drop' or 'keep'")
  expect_error(clean.conc.blq(d.test$conc, d.test$time,
                              conc.blq=list(
                                first="keep",
                                middle="foo",
                                last="drop"),
                              conc.na="drop"),
               regexp="conc.blq must either be a finite number or the text 'drop' or 'keep'")
  expect_error(clean.conc.blq(d.test$conc, d.test$time,
                              conc.blq=list(
                                first="keep",
                                middle="drop",
                                last="foo"),
                              conc.na="drop"),
               regexp="conc.blq must either be a finite number or the text 'drop' or 'keep'")

  ## If there are BLQ values at the end, it drops them if the
  ## instructions are generic drop.
  d.test <- data.frame(conc=c(1, 0), time=1:2)
  expect_equal(clean.conc.blq(d.test$conc, d.test$time,
                                conc.blq="drop", conc.na="drop"),
               d.test[1,])
  ## And it keeps them if instructed with a list.
  d.test <- data.frame(conc=c(1, 0), time=1:2)
  expect_equal(clean.conc.blq(d.test$conc, d.test$time,
                              conc.blq=list(
                                first="drop",
                                middle="drop",
                                last="keep"),
                              conc.na="drop"),
               d.test)

  ## If there are BLQ values at the beginning and end, it drops those
  ## if given a single instruction.
  d.test <- data.frame(conc=c(0, 1, 0), time=1:3)
  expect_equal(clean.conc.blq(d.test$conc, d.test$time,
                                conc.blq="drop", conc.na="drop"),
               d.test[2,])

  ## If all values are BLQ, drops all rows
  d.test <- data.frame(conc=0, time=1:3)
  expect_equal(clean.conc.blq(d.test$conc, d.test$time,
                                conc.blq="drop", conc.na="drop"),
               d.test[c(),])

  ## If there are BLQ values in the middle, it drops or keeps those or
  ## sets them to a number
  d.test <- data.frame(conc=c(1, 0, 2), time=1:3)
  expect_equal(clean.conc.blq(d.test$conc, d.test$time,
                                conc.blq="drop", conc.na="drop"),
               d.test[-2,])

  d.test <- data.frame(conc=c(1, 0, 2), time=1:3)
  expect_equal(clean.conc.blq(d.test$conc, d.test$time,
                                conc.blq="keep", conc.na="drop"),
               d.test)

  d.test <- data.frame(conc=c(1, 0, 2), time=1:3)
  d.result <- data.frame(conc=c(1, 0.5, 2), time=1:3)
  expect_equal(clean.conc.blq(d.test$conc, d.test$time,
                                conc.blq=0.5, conc.na="drop"),
               d.result)

  ## If there are BLQ values at the beginning, middle, and end, it
  ## only drops all of them or drops them selectively as instructed.
  d.test <- data.frame(conc=c(0, 1, 0, 2, 0), time=1:5)
  expect_equal(clean.conc.blq(d.test$conc, d.test$time,
                                conc.blq="drop", conc.na="drop"),
               d.test[c(2, 4),])
  for (first in c("drop", "keep")) {
    for (middle in c("drop", "keep")) {
      for (last in c("drop", "keep")) {
        expect_equal(clean.conc.blq(d.test$conc, d.test$time,
                                    conc.blq=list(
                                      first=first,
                                      middle=middle,
                                      last=last),
                                    conc.na=0),
                     d.test[c(first %in% "keep",
                              TRUE,
                              middle %in% "keep",
                              TRUE,
                              last %in% "keep"),],
                     info=paste(first, middle, last))
      }
    }
  }
  
  ## When conc.na is 0, it drops those.
  d.test <- data.frame(conc=c(0, 1, NA, 2, 0), time=1:5)
  expect_equal(clean.conc.blq(d.test$conc, d.test$time,
                                conc.blq="drop", conc.na=0),
               d.test[c(2, 4),])

  ## When conc.na is a number, it keeps those
  d.test <- data.frame(conc=c(0, 1, NA, 2, 0), time=1:5)
  d.result <- data.frame(conc=c(0, 1, 0.5, 2, 0), time=1:5)
  expect_equal(clean.conc.blq(d.test$conc, d.test$time,
                                conc.blq="drop", conc.na=0.5),
               d.result[2:4,])

  ## It passes additional to be part of the output data frame
  d.test <- data.frame(conc=c(0, 1, NA, 2, 0), time=1:5, more=6:10)
  d.result <- data.frame(conc=c(0, 1, 0.5, 2, 0), time=1:5, more=6:10)
  expect_equal(clean.conc.blq(d.test$conc, d.test$time,
                              more=d.test$more,
                              conc.blq="drop", conc.na=0.5),
               d.result[2:4,])
  d.test <- data.frame(conc=c(0, 1, NA, 2, 0), time=1:5, more=6:10)
  expect_equal(clean.conc.blq(d.test$conc, d.test$time,
                              more=d.test$more,
                              conc.blq="drop", conc.na="drop"),
               d.test[c(2,4),])
  
})

context("AUC and AUMC")

test_that("pk.calc.auxc", {
  ## #####
  ## Verify input checks

  ## Interval is decreasing
  expect_error(pk.calc.auxc(conc=1:2, time=0:1, interval=2:1, method="linear"),
               regexp="The AUC interval must be increasing")
  ## AUC should start at or after the first measurement and should be
  ## before the last measurement
  expect_warning(pk.calc.auxc(conc=1:2, time=2:3, interval=c(1, 3),
                              method="linear"),
                 regexp="Requesting an AUC range starting \\(1\\) before the first measurement \\(2\\) is not allowed")
  ## But starting before the beginning does return NA (not an error)
  expect_warning(v1 <- pk.calc.auxc(conc=1:2, time=2:3, interval=c(1, 3),
                                    method="linear"))
  expect_equal(v1, NA)

  ## All concentrations are NA, return NA
  expect_warning(v2 <- pk.calc.auxc(conc=c(NA, NA), time=2:3, interval=c(1, 3),
                                    method="linear"))
  expect_equal(v2, NA)
  ## All concentrations are 0, return 0
  expect_equal(pk.calc.auxc(conc=c(0, 0), time=2:3, interval=c(2, 3),
                            method="linear"),
               0)
  ## Concentrations mix 0 and NA, return 0
  expect_equal(pk.calc.auxc(conc=c(NA, 0, NA), time=2:4, interval=c(1, 3),
                            method="linear"),
               0)
})

test_that("pk.calc.auc", {
  nullfun <- function(x) x
  ## Linear AUC when the conc at the end of the interval is above LOQ;
  ## lambda.z is unused
  tests <- list(AUCinf=as.numeric(NA),
                AUClast=1.5,
                AUCall=1.5)
  for (t in names(tests)) {
    ## Note: using this structure ensures that there will not be
    ## excessive warnings during testing.
    if (t == "AUCinf") {
      runner <- expect_warning
    } else {
      runner <- nullfun
    }
    runner(v1 <- pk.calc.auc(conc=c(0, 1, 1),
                             time=0:2,
                             interval=c(0, 2),
                             lambda.z=NA,
                             auc.type=t,
                             method="linear"))
    expect_equal(v1, tests[[t]], info=t)
  }

  
  ## Linear AUC when the conc at the end of the interval is BLQ;
  ## lambda.z is used to extrapolate to the end of the interval.
  ## Since lambda.z is NA, the result is NA.
  tests <- list(AUCinf=as.numeric(NA),
                AUClast=0.5,
                AUCall=1)
  for (t in names(tests)) {
    ## Note: using this structure ensures that there will not be
    ## excessive warnings during testing.
    if (t == "AUCinf") {
      runner <- expect_warning
    } else {
      runner <- nullfun
    }
    runner(v1 <- pk.calc.auc(conc=c(0, 1, 0),
                             time=0:2,
                             interval=c(0, 2),
                             lambda.z=NA,
                             auc.type=t,
                             method="linear"))
    expect_equal(v1, tests[[t]], info=t)
  }

  ## The same when lambda.z is given
  tests <- list(AUCinf=1.5,
                AUClast=0.5,
                AUCall=1)
  for (t in names(tests)) {
    ## Note: using this structure ensures that there will not be
    ## excessive warnings during testing.
    if (t == "AUCinf") {
      runner <- expect_warning
    } else {
      runner <- nullfun
    }
    runner(v1 <- pk.calc.auc(conc=c(0, 1, 0),
                             time=0:2,
                             interval=c(0, 2),
                             lambda.z=1,
                             auc.type=t,
                             method="linear"))
    expect_equal(v1, tests[[t]], info=t)
  }

  ## And when there are multiple BLQ values at the end
  tests <- list(AUCinf=1.5,
                AUClast=0.5,
                AUCall=1)
  for (t in names(tests)) {
    ## Note: using this structure ensures that there will not be
    ## excessive warnings during testing.
    if (t == "AUCinf") {
      runner <- expect_warning
    } else {
      runner <- nullfun
    }
    runner(v1 <- pk.calc.auc(conc=c(0, 1, rep(0, 5)),
                             time=0:6,
                             interval=c(0, 6),
                             lambda.z=1,
                             auc.type=t,
                             method="linear"))
    expect_equal(v1, tests[[t]], info=t)
  }

  ## Confirm that center BLQ points are dropped, kept, or imputed
  ## correctly.  Do this with both "linear" and "lin up/log down"
  tests <- list(
    "linear"=list(
      AUCinf=1+1+0.5+1.5+1.5+1,
      AUClast=1+1+0.5+1.5+1.5,
      AUCall=1+1+0.5+1.5+1.5),
    "lin up/log down"=list(
      AUCinf=1+1+0.5+1.5+1/log(2)+1,
      AUClast=1+1+0.5+1.5+1/log(2),
      AUCall=1+1+0.5+1.5+1/log(2)))
  for (t in names(tests))
    for (n in names(tests[[t]])) {
      ## Note: using this structure ensures that there will not be
      ## excessive warnings during testing.
      if (n == "AUCinf") {
        runner <- expect_warning
      } else {
        runner <- nullfun
      }
      runner(v1 <- pk.calc.auc(conc=c(0, 2, 0, 1, 2, 1),
                               time=0:5,
                               interval=c(0, 5),
                               lambda.z=1,
                               auc.type=n,
                               conc.blq="keep",
                               method=t))
      expect_equal(v1, tests[[t]][[n]], info=paste(t, n))
    }

  ## AUCall looks different when there are BLQs at the end
  tests <- list(
    "linear"=list(
      AUCinf=1+1+0.5+1.5+1.5+1,
      AUClast=1+1+0.5+1.5+1.5,
      AUCall=1+1+0.5+1.5+1.5+0.5),
    "lin up/log down"=list(
      AUCinf=1+1+0.5+1.5+1/log(2)+1,
      AUClast=1+1+0.5+1.5+1/log(2),
      AUCall=1+1+0.5+1.5+1/log(2)+0.5))
  for (t in names(tests))
    for (n in names(tests[[t]])) {
      ## Note: using this structure ensures that there will not be
      ## excessive warnings during testing.
      if (n == "AUCinf") {
        runner <- expect_warning
      } else {
        runner <- nullfun
      }
      runner(v1 <- pk.calc.auc(conc=c(0, 2, 0, 1, 2, 1, 0, 0),
                               time=0:7,
                               interval=c(0, 7),
                               lambda.z=1,
                               auc.type=n,
                               conc.blq="keep",
                               method=t))
      expect_equal(v1, tests[[t]][[n]], info=paste(t, n))
    }

  ## When BLQ in the middle and end are dropped, you get different
  ## answers.  (Note that first BLQ dropping would cause errors due to
  ## starting times differing, so not tested here.)
  tests <- list(
    "linear"=list(
      AUCinf=1+3+1.5+1.5+1,
      AUClast=1+3+1.5+1.5,
      AUCall=1+3+1.5+1.5+0.5),
    "lin up/log down"=list(
      AUCinf=1+2/log(2)+1.5+1/log(2)+1,
      AUClast=1+2/log(2)+1.5+1/log(2),
      AUCall=1+2/log(2)+1.5+1/log(2)+0.5))
  for (t in names(tests))
    for (n in names(tests[[t]])) {
      ## Note: using this structure ensures that there will not be
      ## excessive warnings during testing.
      if (n == "AUCinf") {
        runner <- expect_warning
      } else {
        runner <- nullfun
      }
      runner(v1 <- pk.calc.auc(conc=c(0, 2, 0, 1, 2, 1, 0, 0),
                               time=0:7,
                               interval=c(0, 7),
                               lambda.z=1,
                               auc.type=n,
                               conc.blq=list(
                                 first="keep",
                                 middle="drop",
                                 last="keep"),
                               method=t))
      expect_equal(v1,
                   tests[[t]][[n]],
                   info=paste(t, n))
    }
  
  ## When AUCinf is requested with NA for lambda.z, the result is NA.
  tests <- list(
    "linear"=list(
      AUCinf=as.numeric(NA),
      AUClast=1+3+1.5+1.5,
      AUCall=1+3+1.5+1.5+0.5),
    "lin up/log down"=list(
      AUCinf=as.numeric(NA),
      AUClast=1+2/log(2)+1.5+1/log(2),
      AUCall=1+2/log(2)+1.5+1/log(2)+0.5))
  for (t in names(tests))
    for (n in names(tests[[t]])) {
      ## Note: using this structure ensures that there will not be
      ## excessive warnings during testing.
      if (n == "AUCinf") {
        runner <- expect_warning
      } else {
        runner <- nullfun
      }
      runner(v1 <- pk.calc.auc(conc=c(0, 2, 0, 1, 2, 1, 0, 0),
                               time=0:7,
                               interval=c(0, 7),
                               lambda.z=NA,
                               auc.type=n,
                               conc.blq=list(
                                 first="keep",
                                 middle="drop",
                                 last="keep"),
                               method=t))
      expect_equal(v1,
                   tests[[t]][[n]],
                   info=paste(t, n))
    }
  ## Test NA at the beginning

  ## Test NA at the end
  tests <- list(
    "linear"=list(
      AUCinf=as.numeric(NA),
      AUClast=1+3+1.5+1.5,
      AUCall=1+3+1.5+1.5+1),
    "lin up/log down"=list(
      AUCinf=as.numeric(NA),
      AUClast=1+2/log(2)+1.5+1/log(2),
      AUCall=1+2/log(2)+1.5+1/log(2)+1))
  for (t in names(tests))
    for (n in names(tests[[t]])) {
      ## Note: using this structure ensures that there will not be
      ## excessive warnings during testing.
      if (n == "AUCinf") {
        runner <- expect_warning
      } else {
        runner <- nullfun
      }
      runner(v1 <- pk.calc.auc(conc=c(0, 2, 0, 1, 2, 1, NA, 0),
                               time=0:7,
                               interval=c(0, 7),
                               lambda.z=NA,
                               auc.type=n,
                               conc.blq=list(
                                 first="keep",
                                 middle="drop",
                                 last="keep"),
                               method=t))
      expect_equal(v1,
                   tests[[t]][[n]],
                   info=paste(t, n))
    }
  
  ## Test interpolation of times within the time interval
  tests <- list(
    "linear"=list(
      AUCinf=1+3+1.5+1.5+1,
      AUClast=1+3+1.5+1.5,
      AUCall=1+3+1.5+1.5+0.75),
    "lin up/log down"=list(
      AUCinf=1+2/log(2)+1.5+1/log(2)+1,
      AUClast=1+2/log(2)+1.5+1/log(2),
      AUCall=1+2/log(2)+1.5+1/log(2)+0.5/log(2)))
  for (t in names(tests))
    for (n in names(tests[[t]])) {
      ## Note: using this structure ensures that there will not be
      ## excessive warnings during testing.
      if (n == "AUCinf") {
        runner <- expect_warning
      } else {
        runner <- nullfun
      }
      runner(v1 <- pk.calc.auc(conc=c(0, 2, 0, 1, 2, 1, 0),
                               time=c(0:5, 7),
                               interval=c(0, 6),
                               lambda.z=1,
                               auc.type=n,
                               conc.blq=list(
                                 first="keep",
                                 middle="drop",
                                 last="keep"),
                               method=t))
      expect_equal(v1,
                   tests[[t]][[n]],
                   info=paste(t, n))
    }
  
  ## Confirm warning with beginning of interval before the beginning of
  ## time
  tests <- list(
    "linear"=list(
      AUCinf=NA,
      AUClast=NA,
      AUCall=NA),
    "lin up/log down"=list(
      AUCinf=NA,
      AUClast=NA,
      AUCall=NA))
  for (t in names(tests))
    for (n in names(tests[[t]])) {
      ## Note: using this structure ensures that there will not be
      ## excessive warnings during testing.
      expect_warning(
        v1 <- pk.calc.auc(conc=c(0, 2, 0, 1, 2, 1, 0),
                          time=c(0:5, 7),
                          interval=c(-1, 6),
                          lambda.z=1,
                          auc.type=n,
                          conc.blq=list(
                            first="keep",
                            middle="drop",
                            last="keep"),
                          method=t))
      expect_equal(v1,
                   tests[[t]][[n]],
                   info=paste(t, n))
    }
  
  ## Requesting an AUC interval that starts after the last measurement
  ## results in a warning, but will provide an answer (to be tested
  ## elsewhere).
  expect_warning(pk.calc.auc(conc=1:2, time=0:1, interval=2:3,
                             method="linear"),
                 regexp="AUC start time \\(2\\) is after the maximum observed time \\(1\\)")

  ## Confirm error with NA at the same time as the beginning of the
  ## interval
  tests <- list(
    "linear"=list(
      AUCinf=NA,
      AUClast=NA,
      AUCall=NA),
    "lin up/log down"=list(
      AUCinf=NA,
      AUClast=NA,
      AUCall=NA))
  for (t in names(tests))
    for (n in names(tests[[t]])) {
      ## Note: using this structure ensures that there will not be
      ## excessive warnings during testing.
      if (n == "AUCinf") {
        runner <- expect_warning
      } else {
        runner <- nullfun
      }
      expect_warning(
        v1 <- pk.calc.auc(conc=c(NA, 2, 0, 1, 2, 1, 0),
                          time=c(0:5, 7),
                          interval=c(0, 6),
                          lambda.z=1,
                          auc.type=n,
                          conc.blq=list(
                            first="keep",
                            middle="drop",
                            last="keep"),
                          method=t))
      expect_equal(v1,
                   tests[[t]][[n]],
                   info=paste(t, n))
    }

  ## Confirm error with concentration and time not of equal lengths
  expect_error(pk.calc.auc(conc=c(1, 2, 3), time=c(1, 2)),
               regexp="Conc and time must be the same length")

  ## Confirm error with time not monotonically increasing (less than)
  expect_error(pk.calc.auc(conc=c(1, 2, 3), time=c(1, 2, 1)),
               regexp="Time must be monotonically increasing")

  ## Confirm error with time not monotonically increasing (equal)
  expect_error(pk.calc.auc(conc=c(1, 2, 3), time=c(1, 2, 2)),
               regexp="Time must be monotonically increasing")

  ## Confirm that AUC method checking works
  expect_error(pk.calc.auc(conc=c(1, 2, 3), time=c(1, 2, 3), method="wrong"),
               regex='should be one of',
               info="Method names are confirmed for pk.calc.auc")

  ## Confirm that everything works even when check is FALSE
  expect_equal(
    pk.calc.auc.inf(conc=c(0, 1, 1, 0),
                    time=0:3,
                    interval=c(0, Inf),
                    lambda.z=1,
                    check=FALSE,
                    method="linear"),
    pk.calc.auc.inf(conc=c(0, 1, 1, 0),
                    time=0:3,
                    interval=c(0, Inf),
                    lambda.z=1,
                    method="linear"))
})

test_that("pk.calc.auc.last", {
  expect_equal(
    pk.calc.auc.last(conc=c(0, 1, 1, 0),
                     time=0:3,
                     interval=c(0, 3),
                     method="linear"),
    pk.calc.auc(conc=c(0, 1, 1, 0),
                time=0:3,
                interval=c(0, 3),
                auc.type="AUClast",
                method="linear"))
  expect_error(
    pk.calc.auc.last(conc=c(0, 1, 1, 0),
                     time=0:3,
                     interval=c(0, 3),
                     auc.type="foo",
                     method="linear"),
    regexp="auc.type cannot be changed when calling pk.calc.auc.last, please use pk.calc.auc")
})

test_that("pk.calc.auc.inf", {
  ## Just ensuring that it is a simple wrapper.  Computation testing
  ## is done elsewhere.
  expect_equal(
    pk.calc.auc.inf(conc=c(0, 1, 1, 0),
                    time=0:3,
                    interval=c(0, Inf),
                    lambda.z=1,
                    method="linear"),
    pk.calc.auc(conc=c(0, 1, 1, 0),
                time=0:3,
                interval=c(0, Inf),
                lambda.z=1,
                auc.type="AUCinf",
                method="linear"))
  expect_error(
    pk.calc.auc.inf(conc=c(0, 1, 1, 0),
                    time=0:3,
                    interval=c(0, 3),
                    auc.type="foo",
                    method="linear"),
    regexp="auc.type cannot be changed when calling pk.calc.auc.inf, please use pk.calc.auc")
})

test_that("pk.calc.auc.all", {
  ## Just ensuring that it is a simple wrapper.  Computation testing
  ## is done elsewhere.
  expect_equal(
    pk.calc.auc.all(conc=c(0, 1, 1, 0),
                    time=0:3,
                    interval=c(0, 3),
                    method="linear"),
    pk.calc.auc(conc=c(0, 1, 1, 0),
                time=0:3,
                interval=c(0, 3),
                auc.type="AUCall",
                method="linear"))
  expect_error(
    pk.calc.auc.all(conc=c(0, 1, 1, 0),
                    time=0:3,
                    interval=c(0, 3),
                    auc.type="foo",
                    method="linear"),
    regexp="auc.type cannot be changed when calling pk.calc.auc.all, please use pk.calc.auc")
})

test_that("pk.calc.aumc", {
  expect_equal(
    pk.calc.aumc(
      conc=c(0, 1, 1, 0.5),
      time=0:3,
      interval=c(0, 3),
      method="linear"),
    3.75)
  expect_equal(
    pk.calc.aumc(
      conc=c(0, 1, 1, 0.5),
      time=0:3,
      interval=c(0, 3),
      method="lin up/log down"),
    2-0.5/log(0.5)+0.5/(log(0.5)^2))
  expect_equal(
    pk.calc.aumc(
      conc=c(0, 1, 1, 0.5),
      time=0:3,
      interval=c(0, Inf),
      auc.type="AUCinf",
      lambda.z=1,
      method="lin up/log down"),
    2-0.5/log(0.5)+0.5/(log(0.5)^2)+1.5+0.5)
})


test_that("pk.calc.aumc.inf", {
  ## Just ensuring that it is a simple wrapper.  Computation testing
  ## is done elsewhere.
  expect_equal(
    pk.calc.aumc.inf(conc=c(0, 1, 1, 0),
                     time=0:3,
                     interval=c(0, Inf),
                     lambda.z=1,
                     method="linear"),
    pk.calc.aumc(conc=c(0, 1, 1, 0),
                 time=0:3,
                 interval=c(0, Inf),
                 lambda.z=1,
                 auc.type="AUCinf",
                 method="linear"))
  expect_error(
    pk.calc.aumc.inf(conc=c(0, 1, 1, 0),
                    time=0:3,
                    interval=c(0, 3),
                    auc.type="foo",
                    method="linear"),
    regexp="auc.type cannot be changed when calling pk.calc.aumc.inf, please use pk.calc.aumc")
})

test_that("pk.calc.aumc.all", {
  ## Just ensuring that it is a simple wrapper.  Computation testing
  ## is done elsewhere.
  expect_equal(
    pk.calc.aumc.all(conc=c(0, 1, 1, 0),
                    time=0:3,
                    interval=c(0, 3),
                    method="linear"),
    pk.calc.aumc(conc=c(0, 1, 1, 0),
                 time=0:3,
                 interval=c(0, 3),
                 auc.type="AUCall",
                 method="linear"))
  expect_error(
    pk.calc.aumc.all(conc=c(0, 1, 1, 0),
                    time=0:3,
                    interval=c(0, 3),
                    auc.type="foo",
                    method="linear"),
    regexp="auc.type cannot be changed when calling pk.calc.aumc.all, please use pk.calc.aumc")
})

test_that("pk.calc.aumc.last", {
  ## Just ensuring that it is a simple wrapper.  Computation testing
  ## is done elsewhere.
  expect_equal(
    pk.calc.aumc.last(conc=c(0, 1, 1, 0),
                    time=0:3,
                    interval=c(0, 3),
                    method="linear"),
    pk.calc.aumc(conc=c(0, 1, 1, 0),
                 time=0:3,
                 interval=c(0, 3),
                 auc.type="AUClast",
                 method="linear"))
  expect_error(
    pk.calc.aumc.last(conc=c(0, 1, 1, 0),
                    time=0:3,
                    interval=c(0, 3),
                    auc.type="foo",
                    method="linear"),
    regexp="auc.type cannot be changed when calling pk.calc.aumc.last, please use pk.calc.aumc")
})

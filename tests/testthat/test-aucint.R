context("aucint")

test_that("AUCint gives errors appropriately", {
  expect_error(pk.calc.aucint(conc=1, time=1),
               regexp="If interval is not given, start and end must be given.",
               fixed=TRUE,
               info="Neither start nor end is given.")
  expect_error(pk.calc.aucint(conc=1, time=1, end=1),
               regexp="If interval is not given, start and end must be given.",
               fixed=TRUE,
               info="start is not given.")
  expect_error(pk.calc.aucint(conc=1, time=1, start=1),
               regexp="If interval is not given, start and end must be given.",
               fixed=TRUE,
               info="end is not given.")
  expect_error(pk.calc.aucint(conc=1, time=1, start=1:2, end=1),
               regexp="start must be a scalar",
               fixed=TRUE)
  expect_error(pk.calc.aucint(conc=1, time=1, start=1, end=1:2),
               regexp="end must be a scalar",
               fixed=TRUE)

  expect_error(pk.calc.aucint(conc=1, time=1, start="A", end=1),
               regexp="start must be a number",
               fixed=TRUE,
               info="start not a number")
  expect_error(pk.calc.aucint(conc=1, time=1, start=factor("A"), end=1),
               regexp="start must be a number",
               fixed=TRUE,
               info="start not a number (factor)")
  expect_error(pk.calc.aucint(conc=1, time=1, start=Inf, end=1),
               regexp="interval beginning (or start) must be finite",
               fixed=TRUE,
               info="start not finite")

  expect_error(pk.calc.aucint(conc=1, time=1, end="A", start=1),
               regexp="end must be a number",
               fixed=TRUE,
               info="end not a number")
  expect_error(pk.calc.aucint(conc=1, time=1, end=factor("A"), start=1),
               regexp="end must be a number",
               fixed=TRUE,
               info="end not a number (factor)")

  expect_error(pk.calc.aucint(conc=1, time=1, start=1, end=2, interval=c(1, 2)),
               regexp="start and end cannot be given if interval is given",
               fixed=TRUE)
  expect_error(pk.calc.aucint(conc=1, time=1, interval=1:3),
               regexp="interval must be a vector with 2 elements",
               fixed=TRUE)
  expect_error(pk.calc.aucint(conc=1, time=1, interval=c("A", "B")),
               regexp="interval must be numeric",
               fixed=TRUE,
               info="interval not a number")
  expect_error(pk.calc.aucint(conc=1, time=1, interval=factor(c("A", "B"))),
               regexp="interval must be numeric",
               fixed=TRUE,
               info="interval not a number (factor)")
  expect_error(pk.calc.aucint(conc=1, time=1, interval=c(Inf, 1)),
               regexp="interval beginning (or start) must be finite",
               fixed=TRUE,
               info="interval not finite (start)")
  expect_error(pk.calc.aucint(conc=1, time=1, interval=c(-Inf, 1)),
               regexp="interval beginning (or start) must be finite",
               fixed=TRUE,
               info="interval not finite (start)")
  expect_error(pk.calc.aucint(conc=1, time=1, interval=c(1, 1)),
               regexp="interval start must be before interval end.",
               fixed=TRUE,
               info="interval must be increasing (equal)")
  expect_error(pk.calc.aucint(conc=1, time=1, interval=c(1, 0)),
               regexp="interval start must be before interval end.",
               fixed=TRUE,
               info="interval must be increasing (decreasing)")
})

test_that("AUCint gives the same value when no interpolation/extrapolation is required", {
  tmpdata <- data.frame(conc=c(8, 4, 2, 1),
                        time=0:3)
  expect_equal(pk.calc.aucint(conc=tmpdata$conc, time=tmpdata$time,
                              interval=c(0, 3)),
               pk.calc.auc(conc=tmpdata$conc, time=tmpdata$time,
                           interval=c(0, 3)),
               info="No interpolation/extrapolation is equivalent to normal AUC")
  expect_equal(pk.calc.aucint(conc=tmpdata$conc, time=tmpdata$time,
                              start=0, end=3),
               pk.calc.auc(conc=tmpdata$conc, time=tmpdata$time,
                           interval=c(0, 3)),
               info="Giving interval and start+end are the same, no interp/extrap (test 1)")
  expect_equal(pk.calc.aucint(conc=tmpdata$conc, time=tmpdata$time,
                              start=0, end=2),
               pk.calc.auc(conc=tmpdata$conc, time=tmpdata$time,
                           interval=c(0, 2)),
               info="Giving interval and start+end are the same, no interp/extrap (test 2)")
})

test_that("AUCint gives a warning and NA when it cannot interpolate or extrapolate a value", {
  tmpdata <- data.frame(conc=c(8, 4, 2, 1),
                        time=0:3)
  expect_warning(
    over_dose <-
      pk.calc.aucint(conc=tmpdata$conc,
                     time=tmpdata$time,
                     interval=c(0, 4),
                     time.dose=1.5,
                     lambda.z=NA,
                     auc.type="AUCinf"),
    regexp="Some interpolated/extrapolated concentration values are missing (may be due to interpolating or extrapolating over a dose with lambda.z=NA). Time points with missing data are:  1.5, 4",
    fixed=TRUE,
    info="warned when integrating over a dose with lambda.z=NA")
  expect_equal(over_dose, NA_real_,
               info="When you cannot integrate over a dose, you get NA")

  expect_warning(
    before_time <-
      pk.calc.aucint(conc=tmpdata$conc,
                     time=tmpdata$time,
                     interval=c(-1, 4),
                     time.dose=c(-1.5, -0.5),
                     lambda.z=1,
                     auc.type="AUCinf"),
    regexp="Some interpolated/extrapolated concentration values are missing Time points with missing data are:  -1, -0.5",
    fixed=TRUE,
    info="warned when integrating over a dose with lambda.z=NA")
  expect_equal(before_time, NA_real_,
               info="When you cannot interpolate a point, you get NA")
})

test_that("AUCint respects auc.type and does the correct calculations for each AUC type", {
  tmpdata <- data.frame(conc=c(8, 4, 2, 1),
                        time=0:3)
  tmpdata_blq <- data.frame(conc=c(8, 4, 2, 1, 0),
                            time=0:4)
  expect_equal(pk.calc.aucint(conc=tmpdata$conc, time=tmpdata$time,
                              interval=c(0, 3), auc.type="AUClast"),
               pk.calc.aucint(conc=tmpdata$conc, time=tmpdata$time,
                              interval=c(0, 3)),
               info="Default AUC type is AUClast")

  expect_equal(pk.calc.aucint(conc=tmpdata$conc, time=tmpdata$time,
                              interval=c(0, 4), auc.type="AUCinf", clast=1, lambda.z=log(2)),
               pk.calc.auc(conc=c(tmpdata$conc, 0.5), time=c(tmpdata$time, 4),
                           interval=c(0, 4)),
               info="AUCinf is traced")
  expect_equal(pk.calc.aucint(conc=tmpdata$conc, time=tmpdata$time,
                              interval=c(0, 4), auc.type="AUCinf", clast=2, lambda.z=log(2)),
               pk.calc.auc(conc=tmpdata$conc, time=tmpdata$time,
                           interval=c(0, 3)) +
                 pk.calc.auc(conc=c(2, 1), time=c(3, 4),
                             interval=c(3, 4)),
               info="AUCinf is traced with clast respected")
  expect_equal(pk.calc.aucint(conc=tmpdata$conc, time=tmpdata$time,
                              interval=c(0, 4), auc.type="AUCinf", clast=1, lambda.z=log(2)*2),
               pk.calc.auc(conc=c(tmpdata$conc, 0.25), time=c(tmpdata$time, 4),
                           interval=c(0, 4)),
               info="AUCinf is traced with lambda.z respected")
  expect_equal(pk.calc.aucint(conc=tmpdata$conc, time=tmpdata$time,
                              interval=c(0, 4), auc.type="AUCinf", clast=2, lambda.z=2*log(2)),
               pk.calc.auc(conc=tmpdata$conc, time=tmpdata$time,
                           interval=c(0, 3)) +
                 pk.calc.auc(conc=c(2, 0.5), time=c(3, 4),
                             interval=c(3, 4)),
               info="AUCinf is traced with clast and lambda.z respected")

  expect_equal(pk.calc.aucint(conc=tmpdata$conc, time=tmpdata$time,
                              interval=c(0, 3), auc.type="AUCall"),
               pk.calc.auc(conc=tmpdata$conc, time=tmpdata$time,
                           interval=c(0, 3), auc.type="AUCall"),
               info="AUCall is the same as AUClast when no BLQ follow tlast (both AUCall)")
  expect_equal(pk.calc.aucint(conc=tmpdata$conc, time=tmpdata$time,
                              interval=c(0, 3), auc.type="AUCall"),
               pk.calc.auc(conc=tmpdata$conc, time=tmpdata$time,
                           interval=c(0, 3), auc.type="AUClast"),
               info="AUCall is the same as AUClast when no BLQ follow tlast (test AUClast)")
  expect_equal(pk.calc.aucint(conc=tmpdata_blq$conc, time=tmpdata_blq$time,
                              interval=c(0, 4), auc.type="AUCall"),
               pk.calc.auc(conc=tmpdata_blq$conc, time=tmpdata_blq$time,
                           interval=c(0, 4), auc.type="AUCall"),
               info="AUCall is the same the normal calculation when no interpolation/extrapolation happens")
  expect_equal(pk.calc.aucint(conc=tmpdata_blq$conc, time=tmpdata_blq$time,
                              interval=c(0, 3.5), auc.type="AUCall"),
               pk.calc.auc(conc=c(tmpdata$conc, 0.5), time=c(tmpdata$time, 3.5),
                           interval=c(0, 4), auc.type="AUClast"),
               info="AUCall traces correctly")
})

test_that("aucint respects doses", {
  tmpdata <- data.frame(conc=c(8, 4, 2, 1),
                        time=0:3)
  tmpdata_blq <- data.frame(conc=c(8, 4, 2, 1, 0),
                            time=0:4)
  time.dose_at <- 1
  time.dose_between <- 1.5
  time.dose_after_tlast <- 3.5
  time.dose_after_all <- 4.5
  expect_equal(pk.calc.aucint(conc=tmpdata_blq$conc, time=tmpdata_blq$time,
                              interval=c(0, 4), auc.type="AUCall", time.dose=time.dose_at),
               pk.calc.aucint(conc=tmpdata_blq$conc, time=tmpdata_blq$time,
                              interval=c(0, 4), auc.type="AUCall"),
               info="Calculation with dosing at the same time as an observation causes no change.")
  expect_equal(pk.calc.aucint(conc=tmpdata_blq$conc, time=tmpdata_blq$time,
                              interval=c(0, 4), auc.type="AUCall", time.dose=time.dose_after_all),
               pk.calc.aucint(conc=tmpdata_blq$conc, time=tmpdata_blq$time,
                              interval=c(0, 4), auc.type="AUCall"),
               info="Calculation with dosing at a time after all observations causes no change.")
  ## FINDME
  # expect_warning(no_lambda.z_dose <-
  #                  pk.calc.aucint(conc=tmpdata_blq$conc, time=tmpdata_blq$time,
  #                                 interval=c(0, 4), auc.type="AUCall", time.dose=time.dose_between),
  #                regexp="Cannot interpolate/extrapolate to dosing times in pk.calc.aucint without lambda.z",
  #                fixed=TRUE)
  # expect_equal(no_lambda.z_dose, NA_real_,
  #              info="Calculating with dosing times, without observations, and without lambda.z")
  expect_equal(pk.calc.aucint(conc=tmpdata_blq$conc, time=tmpdata_blq$time,
                              time.dose=5,
                              interval=c(0, 4), auc.type="AUCall"),
               pk.calc.aucint(conc=tmpdata_blq$conc, time=tmpdata_blq$time,
                              interval=c(0, 4), auc.type="AUCall"),
               info="Calculation with dosing at a time after all observations causes no change.")
})

test_that("aucint works with infinite intervals", {
  tmpdata <- data.frame(conc=c(8, 4, 2, 1),
                        time=0:3)
  expect_equal(pk.calc.aucint.last(conc=tmpdata$conc, time=tmpdata$time, start=0, end=Inf),
               pk.calc.auc.last(conc=tmpdata$conc, time=tmpdata$time),
               info="Simple AUClast = aucint.last")
  expect_equal(pk.calc.aucint.all(conc=tmpdata$conc, time=tmpdata$time, start=0, end=Inf),
               pk.calc.auc.all(conc=tmpdata$conc, time=tmpdata$time),
               info="Simple AUCall = aucint.all")
  expect_equal(pk.calc.aucint.inf.obs(conc=tmpdata$conc, time=tmpdata$time,
                                      start=0, end=Inf,
                                      clast.obs=1, lambda.z=log(2)),
               pk.calc.auc.inf.obs(conc=tmpdata$conc, time=tmpdata$time,
                                   clast.obs=1, lambda.z=log(2)),
               info="Simple AUCinf.obs = aucint.inf.obs")
  expect_equal(pk.calc.aucint.inf.pred(conc=tmpdata$conc, time=tmpdata$time,
                                       start=0, end=Inf,
                                       clast.pred=2, lambda.z=log(2)),
               pk.calc.auc.inf.pred(conc=tmpdata$conc, time=tmpdata$time,
                                    clast.pred=2, lambda.z=log(2)),
               info="Simple AUCinf.pred = aucint.inf.pred")
})
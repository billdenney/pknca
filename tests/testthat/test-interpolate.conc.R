context("Interpolation and extrapolation of concentration")

test_that("extrapolate.conc", {
  ## Confirm that the ecxtrap.method must be AUCinf, AUClast, AUCall
  expect_error(extrapolate.conc(conc=1,
                                time=1,
                                time.out=2,
                                extrap.method="wrong"),
               regexp="extrap.method must be one of 'AUCinf', 'AUClast', or 'AUCall'")

  ## Confirm that time.out may only be a scalar
  expect_error(extrapolate.conc(conc=1,
                                time=1,
                                time.out=c(2, 3),
                                extrap.method="AUCinf"),
               regexp="Only one time.out value may be estimated at once.")

  ## Ensure that if Clast is NA that the extrapolated concentration is
  ## NA.
  expect_warning(v1 <-
    extrapolate.conc(conc=NA,
                     time=1,
                     time.out=2,
                     extrap.method="AUCinf"))
  expect_equal(v1, NA)

  ## Confirm that it is an error to extrapolate at or before Tlast
  expect_error(extrapolate.conc(conc=1,
                                time=1,
                                time.out=0.5,
                                extrap.method="AUCinf"),
               regexp="extrapolate.conc can only work beyond Tlast, please use interp.extrap.conc to combine both interpolation and extrapolation.")
  expect_error(extrapolate.conc(conc=1,
                                time=1,
                                time.out=1,
                                extrap.method="AUCinf"),
               regexp="extrapolate.conc can only work beyond Tlast, please use interp.extrap.conc to combine both interpolation and extrapolation.")

  ## Extrapolating AUClast is always 0
  expect_equal(extrapolate.conc(conc=c(0, 1, 0),
                                time=1:3,
                                time.out=4,
                                extrap.method="AUClast"),
               0)
  expect_equal(extrapolate.conc(conc=c(0, 1, 0),
                                time=1:3,
                                time.out=2.5,
                                extrap.method="AUClast"),
               0)
  
  ## Extrapolating AUCall after the last value and it is 0 is 0
  expect_equal(extrapolate.conc(conc=c(0, 1, 0),
                                time=1:3,
                                time.out=4,
                                extrap.method="AUCall"),
               0)

  ## Extrapolating AUCall after the last value and it is nonzero is 0
  expect_equal(extrapolate.conc(conc=c(0, 1, 1),
                                time=1:3,
                                time.out=4,
                                lambda.z=2,
                                extrap.method="AUCall"),
               0)

  ## Extrapolating AUCall between Tlast and the first BLQ is a linear
  ## interpolation.
  expect_equal(extrapolate.conc(conc=c(0, 1, 0),
                                time=1:3,
                                time.out=2.5,
                                extrap.method="AUCall"),
               0.5)

  ## Extrapolating AUCinf between the Clast value and a 0 ignores the
  ## 0.
  expect_equal(extrapolate.conc(conc=c(0, 1, 0),
                                time=1:3,
                                time.out=2.5,
                                lambda.z=1,
                                extrap.method="AUCinf"),
               1*exp(-1*0.5))

  ## A slightly less trivial example
  expect_equal(extrapolate.conc(conc=c(0, 5, 0),
                                time=1:3,
                                time.out=2.5,
                                lambda.z=3,
                                extrap.method="AUCinf"),
               5*exp(-3*0.5))

  ## Extrapolating AUCinf with lambda.z=NA is NA.
  expect_equal(extrapolate.conc(conc=c(0, 5, 0),
                                time=1:3,
                                time.out=2.5,
                                lambda.z=NA,
                                extrap.method="AUCinf"),
               as.numeric(NA))

  ## Extrapolating with all NA is NA.
  expect_warning(v1 <-
    extrapolate.conc(conc=rep(NA, 3),
                     time=1:3,
                     time.out=2.5,
                     lambda.z=NA,
                     extrap.method="AUCinf"))
  expect_equal(v1, NA)

  ## Ensure that extrapolation beyond the last point works if the last
  ## point is 0
  extrapolations <- list(AUCinf=exp(-2),
                         AUCall=0,
                         AUClast=0)
  for (n in names(extrapolations))
    expect_equal(extrapolate.conc(conc=c(0, 1, 0.5, 1, 0),
                                  time=0:4, 5,
                                  lambda.z=1,
                                  extrap.method=n),
                 extrapolations[[n]],
                 info=n)

  ## Ensure that extrapolation to the last point of 0 works
  extrapolations <- list(AUCinf=exp(-1),
                         AUCall=0,
                         AUClast=0)
  for (n in names(extrapolations))
    expect_equal(extrapolate.conc(conc=c(0, 1, 0.5, 1, 0),
                                  time=0:4, 4,
                                  lambda.z=1,
                                  extrap.method=n),
                 extrapolations[[n]],
                 info=n)

  ## Ensure that extrapolation between the last points with the last point 0 works
  extrapolations <- list(AUCinf=exp(-0.5),
                         AUCall=0.5,
                         AUClast=0)
  for (n in names(extrapolations))
    expect_equal(extrapolate.conc(conc=c(0, 1, 0.5, 1, 0),
                                  time=0:4, 3.5,
                                  lambda.z=1,
                                  extrap.method=n),
                 extrapolations[[n]],
                 info=n)
})

test_that("interpolate.conc", {
  ## Confirm that interpolating to a given point results in the point
  ## itself
  interpolations <- list(linear=0,
                         "lin up/log down"=0)
  for (n in names(interpolations))
    expect_equal(interpolate.conc(conc=c(0, 1),
                                  time=0:1,
                                  time.out=0,
                                  interp.method=n),
                 interpolations[[n]],
                 info=n)

  interpolations <- list(linear=1,
                         "lin up/log down"=1)
  for (n in names(interpolations))
    expect_equal(interpolate.conc(conc=c(0, 1),
                                  time=0:1,
                                  time.out=1,
                                  interp.method=n),
                 interpolations[[n]],
                 info=n)


  ## Confirm that interpolating to a given point when that point is NA
  ## gives a interpolates as expected.
  interpolations <- list(linear=1,
                         "lin up/log down"=1)
  for (n in names(interpolations))
    expect_equal(interpolate.conc(conc=c(0, 1, NA, 1),
                                  time=0:3,
                                  time.out=2,
                                  interp.method=n),
                 interpolations[[n]],
                 info=n)

  ## Slightly less trivial tests
  interpolations <- list(linear=0.5,
                         "lin up/log down"=0.5)
  for (n in names(interpolations))
    expect_equal(interpolate.conc(conc=c(0, 1, 0),
                                  time=0:2,
                                  time.out=0.5,
                                  interp.method=n),
                 interpolations[[n]],
                 info=n)

  interpolations <- list(linear=1.75,
                         "lin up/log down"=2^0.75)
  for (n in names(interpolations))
    expect_equal(interpolate.conc(conc=c(0, 2, 1),
                                  time=0:2,
                                  time.out=1.25,
                                  interp.method=n),
                 interpolations[[n]],
                 info=n)

  ## Relative time is the only thing that matters.  Negative time
  ## values are fine.
  interpolations <- list(linear=1.75,
                         "lin up/log down"=2^0.75)
  for (n in names(interpolations))
    expect_equal(interpolate.conc(conc=c(0, 2, 1),
                                  time=seq(-10, -8, by=1),
                                  time.out=-8.75,
                                  interp.method=n),
                 interpolations[[n]],
                 info=n)

  ## Ensure that linear-up/log down interpolation works linearly for
  ## the "up" even if it is after a down.
  interpolations <- list(linear=0.25,
                         "lin up/log down"=0.25)
  for (n in names(interpolations))
    expect_equal(interpolate.conc(conc=c(0, 1, 0, 1, 0),
                                  time=0:4,
                                  time.out=2.25,
                                  conc.blq="keep",
                                  interp.method=n),
                 interpolations[[n]],
                 info=n)

  ## Ensure that interpolation with linear works with a BLQ in the
  ## middle: going down and up from that BLQ
  interpolations <- list(linear=0.5,
                         "lin up/log down"=0.5)
  for (n in names(interpolations))
    expect_equal(interpolate.conc(conc=c(0, 1, 0, 1, 0),
                                  time=0:4,
                                  time.out=1.5,
                                  conc.blq="keep",
                                  interp.method=n),
                 interpolations[[n]],
                 info=n)

  interpolations <- list(linear=0.5,
                         "lin up/log down"=0.5)
  for (n in names(interpolations))
    expect_equal(interpolate.conc(conc=c(0, 1, 0, 1, 0),
                                  time=0:4,
                                  time.out=2.5,
                                  conc.blq="keep",
                                  interp.method=n),
                 interpolations[[n]],
                 info=n)

  ## Ensure that conc.blq is appropriately observed
  interpolations <- list(linear=1,
                         "lin up/log down"=1)
  for (n in names(interpolations))
    expect_equal(interpolate.conc(conc=c(0, 1, 0, 1, 0),
                                  time=0:4,
                                  time.out=2.25,
                                  conc.blq="drop",
                                  interp.method=n),
                 interpolations[[n]],
                 info=n)

  ## Ensure that the "log down" part works, too.
  expect_equal(interpolate.conc(conc=c(0, 1, 0.5, 1, 0),
                                time=0:4,
                                time.out=1.5,
                                interp.method="lin up/log down"),
               exp(mean(log(c(1, 0.5)))))

  ## Ensure that it is an error to extrapolate past the last point.
  expect_error(interpolate.conc(conc=c(0, 1, 0.5, 1, 0),
                                time=0:4,
                                time.out=5,
                                interp.method="lin up/log down"),
               regexp="interpolate.conc can only works through Tlast, please use interp.extrap.conc to combine both interpolation and extrapolation.")

  ## ##############################
  ## Confirm errors that should happen

  ## Confirm that extrapolating before the first time is an error
  expect_error(interpolate.conc(conc=0:1,
                                time=0:1,
                                time.out=-1),
               regexp="Cannot interpolate backward in time")

  ## Confirm that more than one time.out requested is an error (change
  ## documentation when this is not true)
  expect_error(interpolate.conc(conc=0:1,
                                time=0:1,
                                time.out=0:1),
               regexp="Can only interpolate for one time point per function call")

  ## Confirm that invalid interpolation methods are an error.
  expect_error(interpolate.conc(conc=0:1,
                                time=0:1,
                                time.out=0.5,
                                interp.method="this doesn't work"),
               regexp="interp.method must be one of 'linear' or 'lin up/log down'")

  ## FIXME: Future feature: Interpolation backward in time from a 0
  ## starting value can be 0.
})

test_that("interp.extrap.conc", {
  ## Ensure that data checking works correctly
  expect_equal(
    interp.extrap.conc(conc=c(0, 1, 0.5, 1, 0),
                       time=0:4,
                       time.out=1.5,
                       interp.method="lin up/log down"),
    exp(mean(log(c(1, 0.5)))))

  ## Ensure that NA for time.out results in NA output
  expect_warning(v1 <-
    interp.extrap.conc(conc=c(0, 1, 0.5, 1, 0),
                       time=0:4,
                       time.out=c(1.5, NA),
                       interp.method="lin up/log down"))
  expect_equal(v1, c(exp(mean(log(c(1, 0.5)))), NA))

  ## Ensure a warning with NA for time.out
  expect_warning(
    interp.extrap.conc(conc=c(0, 1, 0.5, 1, 0),
                       time=0:4,
                       time.out=c(1.5, NA),
                       interp.method="lin up/log down"),
    regexp="An interpolation/extrapolation time is NA")

})

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

  ## Confirm that extrapolating before the first time uses conc.origin
  expect_equal(interpolate.conc(conc=0:1,
                                time=0:1,
                                time.out=-1),
               0,
               info="conc.origin defaults to zero")
  expect_equal(interpolate.conc(conc=0:1,
                                time=0:1,
                                time.out=-1,
                                conc.origin=NA),
               NA,
               info="conc.origin is honored as NA")
  expect_equal(interpolate.conc(conc=0:1,
                                time=0:1,
                                time.out=-1,
                                conc.origin=5),
               5,
               info="conc.origin is honored as a number")
  
  expect_equal(interpolate.conc(conc=c(NA, 1),
                                time=c(0, 1),
                                time.out=0.5,
                                check=FALSE),
               NA_real_,
               info="Skipping the checks with an NA bounding the interpolation gives NA")
  
  expect_equal(interpolate.conc(conc=c(NA, 1, 2),
                                time=c(0, 1, 2),
                                time.out=1.5,
                                check=FALSE),
               1.5,
               info="Skipping the checks with an NA, but not bounding the interpolation gives the expected value.")

  ## ##############################
  ## Confirm errors that should happen

  ## Change documentation when this is not true
  expect_error(interpolate.conc(conc=0:1,
                                time=0:1,
                                time.out=0:1),
               regexp="Can only interpolate for one time point per function call",
               info="Confirm that more than one time.out requested is an error")

  expect_error(interpolate.conc(conc=0:1,
                                time=0:1,
                                time.out=0.5,
                                interp.method="this doesn't work"),
               regexp=tryCatch(expr={
                 match.arg("foo", choices=c("lin up/log down", "linear"))
               },
               error=function(e) e)$message,
               fixed=TRUE,
               info="Confirm that invalid interpolation methods are an error.")

  expect_error(interpolate.conc(conc=0:1,
                                time=0:1,
                                time.out=0.5,
                                conc.origin=1:2),
               regexp="conc.origin must be a scalar",
               info="conc.origin must be a scalar")
  expect_error(interpolate.conc(conc=0:1,
                                time=0:1,
                                time.out=0.5,
                                conc.origin="A"),
               regexp="conc.origin must be NA or a number \\(and not a factor\\)",
               info="conc.origin must be a number and not a factor (character)")
  expect_error(interpolate.conc(conc=0:1,
                                time=0:1,
                                time.out=0.5,
                                conc.origin=factor("A")),
               regexp="conc.origin must be NA or a number \\(and not a factor\\)",
               info="conc.origin must be a number and not a factor (factor)")
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

  # Check inputs
  expect_error(
    interp.extrap.conc(conc=c(0, 1, 0.5, 1, 0),
                       time=0:4,
                       time.out=c(),
                       interp.method="lin up/log down"),
    regexp="time.out must be a vector with at least one element")
})

test_that("interp.extrap.conc.dose handles all eventualities", {
  event_choices <- unlist(event_choices_interp.extrap.conc.dose, use.names=FALSE)
  eventualities <-
    expand.grid(event_before=setdiff(event_choices, "output_only"),
                event=setdiff(event_choices, "none"),
                event_after=setdiff(event_choices, "output_only"))
  eventualities$method <- NA_character_
  for (nm in names(interp.extrap.conc.dose.select)) {
    mask_selected <- do.call(interp.extrap.conc.dose.select[[nm]]$select, list(x=eventualities))
    expect_true(any(mask_selected),
                info=sprintf("interp.extrap.conc.dose.select[[%s]] matched at least one eventuality", nm))
    expect_true(!any(mask_selected & !is.na(eventualities$method)),
                info=sprintf("interp.extrap.conc.dose.select[[%s]] overlapped with another method.", nm))
    eventualities$method[mask_selected] <- nm
    # if (interactive()) {
    #   cat(sprintf("%d/%d rows filled at %s\n", sum(!is.na(eventualities$method)), nrow(eventualities), nm))
    # }
  }
  expect_false(any(is.na(eventualities$method)),
               info="interp.extrap.conc.dose.select matched all eventualities")
})

test_that("interp.extrap.conc.dose", {
  # Check inputs
  expect_error(interp.extrap.conc.dose(conc=c(0, 1, 2, 1, 0.5, 0.25),
                                       time=c(-1, 1:5),
                                       time.dose=0, route="foo", duration.dose=NA,
                                       time.out=c(-1, -0.1, 0, 0.1, 7), out.after=FALSE),
               regexp="route.dose must be either 'extravascular' or 'intravascular'",
               info="Route must be valid")
  expect_error(interp.extrap.conc.dose(conc=c(0, 1, 2, 1, 0.5, 0.25),
                                       time=c(-1, 1:5),
                                       time.dose=0, route=c("extravascular", "extravascular"), duration.dose=NA,
                                       time.out=c(-1, -0.1, 0, 0.1, 7), out.after=FALSE),
               regexp="route.dose must either be a scalar or the same length as time.dose",
               info="Route must have the correct length")
  expect_error(interp.extrap.conc.dose(conc=c(0, 1, 2, 1, 0.5, 0.25),
                                       time=c(-1, 1:5),
                                       time.dose=0, route="extravascular", duration.dose="A",
                                       time.out=c(-1, -0.1, 0, 0.1, 7), out.after=FALSE),
               regexp="duration.dose must be NA or a number.",
               info="duration.dose must be NA or a number (character).")
  expect_error(interp.extrap.conc.dose(conc=c(0, 1, 2, 1, 0.5, 0.25),
                                       time=c(-1, 1:5),
                                       time.dose=0, route="extravascular", duration.dose=factor("A"),
                                       time.out=c(-1, -0.1, 0, 0.1, 7), out.after=FALSE),
               regexp="duration.dose must be NA or a number.",
               info="duration.dose must be NA or a number (factor).")
  expect_error(interp.extrap.conc.dose(conc=c(0, 1, 2, 1, 0.5, 0.25),
                                       time=c(-1, 1:5),
                                       time.dose=0, route="extravascular", duration.dose=c(1, NA),
                                       time.out=c(-1, -0.1, 0, 0.1, 7), out.after=FALSE),
               regexp="duration.dose must either be a scalar or the same length as time.dose",
               info="duration.dose must match the length of time.dose or be a scalar.")

  expect_equal(interp.extrap.conc.dose(conc=c(0, 1, 2, 1, 0.5, 0.25),
                                       time=c(-1, 1:5),
                                       time.dose=0,
                                       time.out=-2,
                                       check=FALSE),
               interp.extrap.conc.dose(conc=c(0, 1, 2, 1, 0.5, 0.25),
                                       time=c(-1, 1:5),
                                       time.dose=0,
                                       time.out=-2),
               info="Check is respected")

  expect_equal(interp.extrap.conc.dose(conc=c(0, 1, 2, 1, 0.5, 0.25),
                                       time=c(-1, 1:5),
                                       time.dose=0,
                                       time.out=-2),
               structure(0, Method="Before all events"),
               info="Interpolation before all events yields conc.origin which defaults to zero.")

  expect_equal(interp.extrap.conc.dose(conc=c(0, 1, 2, 1, 0.5, 0.25),
                                       time=c(-1, 1:5),
                                       time.dose=0,
                                       conc.origin=NA,
                                       time.out=-2),
               structure(NA_real_, Method="Before all events"),
               info="Interpolation before all events yields conc.origin respecting its input.")

  expect_equal(interp.extrap.conc.dose(conc=c(0, 1, 2, 1, 0.5, 0.25),
                                       time=c(-1, 1:5),
                                       time.dose=0,
                                       time.out=-1),
               structure(0, Method="Observed concentration"),
               info="When there is a concentration measurement at a time point, it is returned.")
  
  expect_equal(interp.extrap.conc.dose(conc=c(0, 1, 2, 1, 0.5, 0.25),
                                       time=c(-1, 1:5),
                                       time.dose=0,
                                       time.out=-0.1),
               structure(0, Method="Extrapolation"),
               info="When the previous measurement is zero and there is no dose between, it is returned.")

  expect_equal(interp.extrap.conc.dose(conc=c(0, 1, 2, 1, 0.5, 0.25),
                                       time=c(-1, 1:5),
                                       time.dose=0,
                                       time.out=0),
               structure(0, Method="Extrapolation"),
               info="When the previous measurement is zero it is at the time of the dose, zero is returned.")

  expect_equal(interp.extrap.conc.dose(conc=c(0, 1, 2, 1, 0.5, 0.25),
                                       time=c(-1, 1:5),
                                       time.dose=0,
                                       time.out=0.1),
               structure(0.1, Method="Dose before, concentration after without a dose"),
               info="Extrapolation to a dose then interpolation between the dose and the next time works.")

  expect_warning(r_double_dose <-
                   interp.extrap.conc.dose(conc=c(0, 1, 2, 1, 0.5, 0.25),
                                           time=c(-1, 1:5),
                                           time.dose=c(0, 0.1),
                                           time.out=0.2),
                 regexp="Cannot interpolate between two doses or after a dose without a concentration after the first dose.",
                 fixed=TRUE,
                 info="Two doses in a row generates a warning")

  expect_equal(r_double_dose,
               structure(NA_real_, Method="Dose before, concentration after without a dose"),
               info="Extrapolation to a dose then interpolation between the dose and the next time gives NA when the dose is NA.")

  expect_equal(interp.extrap.conc.dose(conc=c(0, 1, 2, 1, 0.5, 0.25),
                                       time=c(-1, 1:5),
                                       time.dose=0,
                                       time.out=5),
               structure(0.25, Method="Observed concentration"),
               info="Copy from after the dose.")

  expect_equal(interp.extrap.conc.dose(conc=c(0, 1, 2, 1, 0.5, 0.25),
                                       time=c(-1, 1:5),
                                       time.dose=0,
                                       time.out=7),
               structure(NA_real_, Method="Extrapolation"),
               info="Extrapolation without lambda.z gives NA result")

  expect_equal(interp.extrap.conc.dose(conc=c(0, 1, 2, 1, 0.5, 0.25),
                                       time=c(-1, 1:5),
                                       time.dose=0,
                                       time.out=7,
                                       lambda.z=log(2)),
               structure(0.0625, Method="Extrapolation"),
               info="Extrapolation with lambda.z gives result")
  
  expect_equal(interp.extrap.conc.dose(conc=0:2,
                                       time=0:2,
                                       time.dose=0,
                                       time.out=0.5),
               structure(0.5, Method="Interpolation"),
               info="Interpolation works")

  expect_equal(interp.extrap.conc.dose(conc=c(0:2, 1),
                                       time=0:3,
                                       time.dose=0,
                                       time.out=2.5,
                                       method="linear"),
               structure(sqrt(2), Method="Interpolation"),
               info="Interpolation respects method")
  
  expect_equal(interp.extrap.conc.dose(conc=0:2,
                                       time=0:2,
                                       time.dose=0,
                                       route.dose="intravascular",
                                       time.out=0,
                                       duration.dose=0,
                                       out.after=FALSE),
               structure(0, Method="Observed concentration"),
               info="Observed before IV bolus")

  expect_equal(interp.extrap.conc.dose(conc=0:2,
                                       time=0:2,
                                       time.dose=0,
                                       route.dose="intravascular",
                                       time.out=0,
                                       duration.dose=0,
                                       out.after=TRUE),
               structure(1, Method="Immediately after an IV bolus with a concentration next"),
               info="Observed after IV bolus, one concentration")

  expect_equal(interp.extrap.conc.dose(conc=c(0, 2, 1),
                                       time=0:2,
                                       time.dose=0,
                                       route.dose="intravascular",
                                       time.out=0,
                                       duration.dose=0,
                                       out.after=TRUE),
               structure(4, Method="Immediately after an IV bolus with a concentration next"),
               info="Observed after IV bolus, two concentrations")

  expect_equal(interp.extrap.conc.dose(conc=c(2, 1),
                                       time=1:2,
                                       time.dose=0,
                                       route.dose="intravascular",
                                       time.out=0.5,
                                       duration.dose=0),
               structure(2*sqrt(2), Method="After an IV bolus with a concentration next"),
               info="After IV bolus, two concentrations")

  expect_equal(interp.extrap.conc.dose(conc=2,
                                       time=1,
                                       time.dose=0,
                                       route.dose="intravascular",
                                       time.out=0.5,
                                       duration.dose=0),
               structure(2, Method="After an IV bolus with a concentration next"),
               info="After IV bolus, one concentration")
  
  expect_equal(interp.extrap.conc.dose(conc=c(0, 1, 2, 1, 0.5, 0.25),
                                       time=c(-1, 1:5),
                                       time.dose=0,
                                       time.out=c(-2, 2)),
               structure(c(0, 2),
                         Method=c("Before all events", "Observed concentration")),
               info="Outputs are in the same order as inputs (initially sorted)")
  
  expect_equal(interp.extrap.conc.dose(conc=c(0, 1, 2, 1, 0.5, 0.25),
                                       time=c(-1, 1:5),
                                       time.dose=0,
                                       time.out=c(2, -2)),
               structure(c(2, 0),
                         Method=c("Observed concentration", "Before all events")),
               info="Outputs are in the same order as inputs (reverse sorted time.out)")
})

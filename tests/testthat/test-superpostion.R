context("Superposition")

## The requirements for superposition are defined below to enable
## test-driven development.

#' Dosing must be single dose
#' * A single dose within the time period
#'   * Time period may be a subset of the data from a multi-dose study
#' * Concentrations must be BLQ at the beginning of the interval
#'   * User can over-ride the requirement of BLQ at the beginning of the interval
#' * Half-life
#'   * If the number of dosing intervals times the duration of the dosing interval is less than the observed data, half-life is not required.  (This will likely be rare.)
#'   * Half-life can be specified on input.
#'   * Half-life can be sufficiently determined from the data given (see the half-life cleaning functions for the definition of "sufficient")
#'   * Mixing methods of half-life determination are supported.
#' * Extrapolation can occur via either Clast,obs or Clast,pred.
#'   * If using Clast,pred, it must be given if the half-life is given.
#' * Dosing interval can be specified as either a simple number for dosing interval (tau) or multiple doses within the interval
#'   * Example: Daily dosing would be: 24
#'   * Example: Breakfast/dinner dosing would be: c(10, 24)
#' * How many times to repeat the dosing interval can either be given as a number of taus to repeat or "steady-state"
#'   * If steady-state, then a relative tolerance will be provided to confirm that no concentration changes by a greater ratio from the previous to the current dose.
#' * Doses amount may be provided
#'   * If doses are not given, they are assumed to be the same as the input dose.
#'   * If doses are given, linear dose-proportionality will be assumed.
#'   * Dose amount can be given either as a scalar which will apply to all intervals or a vector which will apply to each interval.
#'     * If dose amount is given as a vector, it must be the same length as the number of intervals input.
#' * User can specify either integration as either linear or linear-up/log-down
#' * Time points in the output are automatically selected as all observed points modulo the dosing interval.
#'   * More time points can be added by user specification, but points cannot be removed.  (Note that this will ensure the highest possible Cmax in the output data.)
#' * The functions work with either individual or grouped data.

test_that("superposition inputs", {
  ## FIXME: Enforce single dose

  ## Enforce BLQ checking for the first data point
  expect_error(superposition(conc=c(1, 2), time=c(0, 1), tau=24),
               regexp="The first concentration must be 0")
  expect_error(superposition(conc=c(NA, 2), time=c(0, 1), tau=24),
               regexp="The first concentration must be 0")

  ## dose.input must be scalar
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), dose.input=c(0, 1)),
               regexp="dose.input must be a scalar")
  ## Ensure that dose.input must be a number
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), dose.input=NA),
               regexp="dose.input must be a number")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), dose.input=factor("1")),
               regexp="dose.input must be a number")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), dose.input="1"),
               regexp="dose.input must be a number")
  ## dose.input must be > 0
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), dose.input=-1),
               regexp="dose.input must be > 0")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), dose.input=0),
               regexp="dose.input must be > 0")
  
  ## tau must be scalar
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=c(0, 1)),
               regexp="tau must be a scalar")
  ## Ensure that tau must be a number
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=NA),
               regexp="tau must be a number")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=factor("1")),
               regexp="tau must be a number")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau="1"),
               regexp="tau must be a number")
  ## tau must be > 0
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=-1),
               regexp="tau must be > 0")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=0),
               regexp="tau must be > 0")

  ## >= 1 dose time
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             dose.times=c()),
               regexp="There must be at least one dose time")
  ## Dose times must be a number
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             dose.times="1"),
               regexp="dose.times must be a number")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             dose.times=factor("1")),
               regexp="dose.times must be a number")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             dose.times=NA),
               regexp="dose.times must be a number")
  ## Dose times nonnegative
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             dose.times=c(-1, 0)),
               regexp="dose.times must be non-negative")
  ## Dose times <= tau
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             dose.times=c(0, 25)),
               regexp="dose.times must be < tau")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             dose.times=25),
               regexp="dose.times must be < tau")

  ## dose.amount without dose.input
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             dose.times=0, dose.amount=1),
               regexp="must give dose.input to give dose.amount")
  ## dose.amount same length as dose.times
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), dose.input=1,
                             tau=24, dose.times=0, dose.amount=c(1, 2)),
               regexp="dose.amount must either be a scalar or match the length of dose.times")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), dose.input=1,
                             tau=24, dose.times=c(0, 1), dose.amount=c(1, 2, 3)),
               regexp="dose.amount must either be a scalar or match the length of dose.times")
  ## dose.amount numeric
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), dose.input=1,
                             tau=24, dose.times=0, dose.amount="1"),
               regexp="dose.amount must be a number")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), dose.input=1,
                             tau=24, dose.times=0, dose.amount=factor("1")),
               regexp="dose.amount must be a number")
  ## dose.amount finite/NA
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), dose.input=1,
                             tau=24, dose.times=0, dose.amount=as.numeric(NA)),
               regexp="dose.amount must be finite and not NA")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), dose.input=1,
                             tau=24, dose.times=0, dose.amount=Inf),
               regexp="dose.amount must be finite and not NA")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), dose.input=1,
                             tau=24, dose.times=0, dose.amount=-Inf),
               regexp="dose.amount must be finite and not NA")
  ## dose.amount non-negative
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), dose.input=1,
                             tau=24, dose.times=0, dose.amount=-1),
               regexp="dose.amount must be positive")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), dose.input=1,
                             tau=24, dose.times=0, dose.amount=0),
               regexp="dose.amount must be positive")

  ## n.tau scalar
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             n.tau=c(1, 2)),
               regexp="n.tau must be a scalar")
  ## n.tau number
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             n.tau="1"),
               regexp="n.tau must be a number")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             n.tau=factor("1")),
               regexp="n.tau must be a number")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             n.tau=NA),
               regexp="n.tau must be a number")
  ## n.tau >= 1
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             n.tau=0),
               regexp="n.tau must be >= 1")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             n.tau=-Inf),
               regexp="n.tau must be >= 1")
  ## n.tau integer
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             n.tau=1.1),
               regexp="n.tau must be an integer or Inf")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             n.tau=1.0000001),
               regexp="n.tau must be an integer or Inf")

  ## lambda.z given and not a scalar
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             n.tau=Inf, lambda.z=c(1, 2)),
               regexp="lambda.z must be a scalar")
  ## lambda.z given and not a number
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             lambda.z="1"),
               regexp="lambda.z must be a number")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             lambda.z=factor("1")),
               regexp="lambda.z must be a number")

  ## clast.pred must be a scalar
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             clast.pred=c(1, 2)),
               regexp="clast.pred must be a scalar")
  ## clast.pred must be a number or logical or NA
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             clast.pred="1"),
               regexp="clast.pred must either be a logical .* or numeric value")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             clast.pred=factor("1")),
               regexp="clast.pred must either be a logical .* or numeric value")

  ## tlast given and not a scalar
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             tlast=c(1, 2)),
               regexp="tlast must be a scalar")
  ## tlast given and not a number
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             tlast="1"),
               regexp="tlast must be a number")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             tlast=factor("1")),
               regexp="tlast must be a number")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             tlast=NA),
               regexp="tlast must be a number")
  ## TODO: test mixing clast.pred and lambda.z

  ## additional.times NA
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             additional.times=NA),
               regexp="No additional.times may be NA \\(to not include any additional.times, enter c\\(\\) as the function argument\\)")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             additional.times=c(2, NA)),
               regexp="No additional.times may be NA \\(to not include any additional.times, enter c\\(\\) as the function argument\\)")
  ## additional.times nonnumeric
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             additional.times="1"),
               regexp="additional.times must be a number")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             additional.times=factor("1")),
               regexp="additional.times must be a number")
  ## additional times < 0
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             additional.times=-1),
               regexp="All additional.times must be nonnegative")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             additional.times=c(-1, 0)),
               regexp="All additional.times must be nonnegative")
  ## Additional times > tau
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             additional.times=25),
               regexp="All additional.times must be <= tau")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             additional.times=c(0, 25)),
               regexp="All additional.times must be <= tau")

  ## steady.state.tol scalar
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             steady.state.tol=c(1, 2)),
               regexp="steady.state.tol must be a scalar")
  ## steady.state.tol numeric
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             steady.state.tol="1"),
               regexp="steady.state.tol must be a number")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             steady.state.tol="1"),
               regexp="steady.state.tol must be a number")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             steady.state.tol=factor("1")),
               regexp="steady.state.tol must be a number")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             steady.state.tol=NA),
               regexp="steady.state.tol must be a number")
  ## steady.state.tol range
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             steady.state.tol=0),
               regexp="steady.state.tol must be between 0 and 1, exclusive.")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             steady.state.tol=1),
               regexp="steady.state.tol must be between 0 and 1, exclusive.")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             steady.state.tol=-1),
               regexp="steady.state.tol must be between 0 and 1, exclusive.")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             steady.state.tol=2),
               regexp="steady.state.tol must be between 0 and 1, exclusive.")

  ## Combinations of lambda.z, clast.pred, tlast
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             clast.pred=1),
               regexp="Either give all or none of the values for these arguments: lambda.z, clast.pred, and tlast")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             lambda.z=1),
               regexp="Either give all or none of the values for these arguments: lambda.z, clast.pred, and tlast")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             tlast=1),
               regexp="Either give all or none of the values for these arguments: lambda.z, clast.pred, and tlast")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             clast.pred=1, lambda.z=1),
               regexp="Either give all or none of the values for these arguments: lambda.z, clast.pred, and tlast")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             clast.pred=1, tlast=1),
               regexp="Either give all or none of the values for these arguments: lambda.z, clast.pred, and tlast")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             lambda.z=1),
               regexp="Either give all or none of the values for these arguments: lambda.z, clast.pred, and tlast")
  expect_error(superposition(conc=c(0, 2), time=c(0, 1), tau=24,
                             lambda.z=1, tlast=1),
               regexp="Either give all or none of the values for these arguments: lambda.z, clast.pred, and tlast")
  
})

test_that("superposition math", {
  ## Superposition without half-life needed and no
  ## interpolation/extrapolation
  c1 <- c(0, 1, 2, 3, 2, 1, 0.5)
  t1 <- 0:6
  expect_equal(superposition(conc=c1, time=t1, tau=3, n.tau=2),
               data.frame(conc=c(3, 3, 3, 3.5),
                          time=0:3))
  ## With half-life extrapolation
  expect_warning(v1 <- superposition(conc=c1, time=t1, tau=3, n.tau=3),
                 regexp="essentially perfect fit")
  expect_equal(v1,
               data.frame(conc=c(3.5, 3.25, 3.125, 3.5625),
                          time=0:3))
  ## To steady-state
  expect_warning(v2 <- superposition(conc=c1, time=t1, tau=3, n.tau=Inf),
                 regexp="essentially perfect fit")
  expect_equal(v2,
               data.frame(conc=c(3.571, 3.286, 3.143, 3.571),
                          time=0:3),
               tol=0.001)

  ## all zeros input gives all zeros output
  expect_equal(superposition(conc=rep(0, 6), time=0:5, tau=24),
               data.frame(conc=0, time=c(0:5, 24)))

  ## Times for output are a collation of 0, tau, dose.times, and
  ## additional.times with time included and shifted for all the
  ## dose.times.
  expect_equal(superposition(conc=rep(0, 6), time=0:5, tau=24,
                             additional.times=c(2.5)),
               data.frame(conc=0, time=c(0, 1, 2, 2.5, 3, 4, 5, 24)))
  expect_equal(superposition(conc=rep(0, 6), time=0:5, tau=25,
                             additional.times=c(2.5)),
               data.frame(conc=0, time=c(0, 1, 2, 2.5, 3, 4, 5, 25)))
  expect_equal(superposition(conc=rep(0, 6), time=0:5, tau=24,
                             additional.times=c(2.5, 3.5)),
               data.frame(conc=0, time=c(0, 1, 2, 2.5, 3, 3.5, 4, 5, 24)))
  expect_equal(superposition(conc=rep(0, 6), time=0:5, tau=24, dose.times=c(0, 1),
                             additional.times=c(2.5, 3.5)),
               data.frame(conc=0, time=c(0, 1, 2, 2.5, 3, 3.5, 4, 5, 6, 24)))
  expect_equal(superposition(conc=rep(0, 6), time=0:5, tau=24,
                             dose.times=c(0, 0.5),
                             additional.times=c(2.5, 3.5)),
               data.frame(conc=0,
                          time=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5,
                            5, 5.5, 24)))

  ## Dose scaling
  expect_warning(v1 <- superposition(conc=c1, time=t1, dose.input=1, tau=24,
                                     dose.times=c(0, 0.5),
                                     dose.amount=1,
                                     additional.times=c(2.5, 3.5)),
                 regexp="essentially perfect fit")
  expect_equal(v1,
               data.frame(conc=c(
                            4.6047e-06,
                            0.5,
                            1.5,
                            2.5,
                            3.5,
                            4.5,
                            5.5,
                            5.4495,
                            4.4495,
                            3.4142,
                            2.4142,
                            1.7071,
                            1.2071,
                            0.85355,
                            4.6047e-06),
                          time=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5,
                            5, 5.5, 6, 6.5, 24)),
               tol=0.001,
               info="Dose scaling with matching input and output doses")

  expect_warning(v2 <- superposition(conc=c1, time=t1, dose.input=1, tau=24,
                                     dose.times=c(0, 0.5),
                                     dose.amount=c(0.5, 5),
                                     additional.times=c(2.5, 3.5)),
                 regexp="essentially perfect fit")
  expect_equal(v2,
               data.frame(conc=c(
                            1.444e-05,
                            0.25,
                            3,
                            5.75,
                            8.5,
                            11.25,
                            14,
                            16.22,
                            13.25,
                            10.71,
                            7.571,
                            5.354,
                            3.786,
                            2.677,
                            1.444e-05),
                          time=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5,
                            5, 5.5, 6, 6.5, 24)),
               tol=0.001,
               info="Dose scaling with different input and output doses")

  expect_warning(v3 <- superposition(conc=c(0, 2, 3, 5, 6, 3, 1, 0),
                                     time=c(0, 0.5, 1, 1.5, 2, 8, 12, 24),
                                     tau=24),
                 regexp="Too few points for half-life calculation \\(min.hl.points=3 with only 2 points\\)")
  expect_equal(v3,
               data.frame(conc=NA, time=c(0, 0.5, 1, 1.5, 2, 8, 12, 24)),
               info="Uncalculable lambda.z with extrapolation to steady-state gives NA conc")

  expect_equal(superposition(conc=c(0, 2, 3, 5, 6, 3, 1, 0),
                             time=c(0, 0.5, 1, 1.5, 2, 8, 12, 24),
                             tau=24,
                             lambda.z=1, clast.pred=1, tlast=12),
               data.frame(conc=c(6.144e-6, 2, 3, 5, 6, 3, 1, 6.144e-6),
                          time=c(0, 0.5, 1, 1.5, 2, 8, 12, 24)),
               tol=0.001,
               info="Uncalculable lambda.z with extrapolation to steady-state with lambda.z given, gives conc values")
})

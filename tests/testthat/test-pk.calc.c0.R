context("C0 Calculations")

test_that("pk.calc.c0", {
  ## Input checks
  expect_error(pk.calc.c0(5:1, 4:0),
               regexp="Time must be monotonically increasing",
               info="conc, time inputs are checked")
  expect_error(pk.calc.c0(5:1, 0:4, time.dose=1:2),
               regexp="time.dose must be a scalar")
  expect_error(pk.calc.c0(5:1, 0:4, time.dose="1"),
               regexp="time.dose must be a number")
  expect_error(pk.calc.c0(5:1, 0:4, method="blah"),
               regexp="should be one of",
               info="method must be valid")
  expect_warning(pk.calc.c0(5:1, 0:4, time.dose=30),
                 regexp="time.dose is after all available data")

  ## Simple calculations
  expect_equal(pk.calc.c0(1:5, 0:4), 1,
               info="When there is a nonzero concentration at the default time of dosing")
  expect_equal(pk.calc.c0(1:5, 0:4, time.dose=1), 2,
               info="When there is a nonzero concentration at the nondefault time of dosing")
  ## It will fall through to subsequent methods if applicable
  expect_equal(pk.calc.c0(c(0, 2, 1, 0.5), 0:3),
               pk.calc.c0.method.logslope(c(0, 2, 1, 0.5), 0:3),
               info="Falls through to second method")
  expect_equal(pk.calc.c0(c(0, 2, 0, 0), 0:3),
               pk.calc.c0.method.c1(c(0, 2, 0, 0), 0:3),
               info="Falls through to third method")
  expect_equal(pk.calc.c0(c(0, 2, 1, 0.5), 0:3, method="c1"),
               pk.calc.c0.method.c1(c(0, 2, 1, 0.5), 0:3),
               info="Respects method order")
  
})

test_that("pk.calc.c0.method.logslope", {
  expect_equal(pk.calc.c0.method.logslope(c(0, 2, 1, 0.5), 0:3), 4)
  expect_equal(pk.calc.c0.method.logslope(c(0, 2, 0, 0), 0:3), NA,
               info="returns NA when C2 is 0")
  expect_equal(pk.calc.c0.method.logslope(c(0, 2, 0, 0.5), 0:3), NA,
               info="returns NA when C2 is 0, even when subsequent are > 0")
  expect_equal(pk.calc.c0.method.logslope(c(0, 2, NA, 0.5), 0:3), 4,
               info="ignores missing concentrations")
  expect_equal(pk.calc.c0.method.logslope(c(0, 2, NA, 0.5), 2:5, time.dose=2), 4,
               info="handles nonzero time.dose")
  expect_equal(pk.calc.c0.method.logslope(c(2, NA, 0.5), 3:5, time.dose=2), 4,
               info="works when time.dose isn't an observed time")
})

test_that("pk.calc.c0.method.c0", {
  expect_equal(pk.calc.c0.method.c0(c(0, 2, 1, 0.5), 0:3), NA)
  expect_equal(pk.calc.c0.method.c0(c(2, NA, 0.5), 3:5, time.dose=3), 2,
               info="works when time.dose matches a nonzero conc")
})

test_that("pk.calc.c0.method.c1", {
  expect_equal(pk.calc.c0.method.c1(c(0, 2, 1, 0.5), 0:3), 2)
  expect_equal(pk.calc.c0.method.c1(c(2, NA, 0.5), 3:5, time.dose=3), 0.5,
               info="works when time.dose matches a nonzero conc and skips over NA")
  expect_equal(pk.calc.c0.method.c1(c(2, NA, 0.5), 3:5, time.dose=30), NA,
               info="returns NA when time.dose is after the last measurement")
  expect_warning(v1 <- pk.calc.c0.method.c1(rep(NA, 3), 3:5, time.dose=1),
                 regexp="All concentration data is missing")
  expect_equal(v1, NA,
               info="returns NA when all inputs are NA")
})

test_that("pk.calc.c0.method.set0", {
  expect_equal(pk.calc.c0.method.set0(), 0)
  expect_equal(pk.calc.c0.method.set0(1, 1), 0)
})

test_that("pk.calc.c0.method.cmin", {
  expect_equal(pk.calc.c0.method.cmin(c(0, 2, 1, 0.5), 0:3), 0)
  expect_equal(pk.calc.c0.method.cmin(c(2, NA, 0.5), 3:5, time.dose=3), 0.5)
})

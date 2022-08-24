test_that("time_above expected errors", {
  expect_error(
    pk.calc.time_above(conc=c(1, 1), time=c(1, 2), conc_above="X", method="linear"),
    regexp='conc_above must be numeric'
  )
  expect_error(
    pk.calc.time_above(conc=c(1, 1), time=c(1, 2), conc_above=1:2, method="linear"),
    regexp='conc_above must be a scalar'
  )
  expect_error(
    pk.calc.time_above(conc=c(1, 1), time=c(1, 2), conc_above=NA, method="linear"),
    regexp='conc_above must not be NA'
  )
  expect_error(
    pk.calc.time_above(time="X", conc_above=5, method="linear"),
    regexp='conc must be given'
  )
  expect_error(
    pk.calc.time_above(conc="X", conc_above=5, method="linear"),
    regexp='time must be given'
  )
  expect_error(
    pk.calc.time_above(conc="X", time="X", conc_above=5, method="linear"),
    regexp='Concentration data must be numeric and not a factor'
  )
  expect_error(
    pk.calc.time_above(conc=5, time="X", conc_above=5, method="linear"),
    regexp='Time data must be numeric and not a factor'
  )
  expect_error(
    pk.calc.time_above(conc=5, time=5, conc_above=5, method="foo"),
    regexp='should be one of'
  )
})

test_that("time_above simple scenarios", {
  expect_equal(
    pk.calc.time_above(conc=5, time=5, conc_above=5, method="linear"),
    structure(NA_real_, exclude="Too few measured concentrations to assess time_above")
  )
})

test_that("time_above linear", {
  expect_equal(
    pk.calc.time_above(conc=c(5, 5), time=1:2, conc_above=5, method="linear"),
    1,
    info="All above, single period"
  )
  expect_equal(
    pk.calc.time_above(conc=c(5, 6, 5), time=1:3, conc_above=5, method="linear"),
    2,
    info="All above, multiple periods"
  )
  expect_equal(
    pk.calc.time_above(conc=c(4, 4), time=1:2, conc_above=5, method="linear"),
    0,
    info="All below, single period"
  )
  expect_equal(
    pk.calc.time_above(conc=c(5, 4, 5), time=1:3, conc_above=5, method="linear"),
    0,
    info="Some equal to conc_above, but all below"
  )
  expect_equal(
    pk.calc.time_above(conc=c(6, 4, 5), time=1:3, conc_above=5, method="linear"),
    0.5
  )
  expect_equal(
    pk.calc.time_above(conc=c(6, NA, 5), time=1:3, conc_above=5, method="linear"),
    2,
    info="NA is ignored"
  )
})

test_that("time_above with 'lin up/log down'", {
  expect_equal(
    pk.calc.time_above(conc=c(6, 4, 5), time=1:3, conc_above=5, method='lin up/log down'),
    (log(6) - log(5))/(log(6) - log(4))
  )
  expect_equal(
    pk.calc.time_above(conc=c(6, 4, 5), time=(1:3)*2, conc_above=5, method='lin up/log down'),
    (log(6) - log(5))/(log(6) - log(4))*2
  )

  expect_equal(
    pk.calc.time_above(conc=c(6, 4, 6), time=1:3, conc_above=5, method='lin up/log down'),
    (log(6) - log(5))/(log(6) - log(4)) + 0.5
  )
  expect_equal(
    pk.calc.time_above(conc=c(6, 4, 6, 0), time=1:4, conc_above=5, method='lin up/log down'),
    (log(6) - log(5))/(log(6) - log(4)) + 0.5 + 1/6
  )
  expect_equal(
    pk.calc.time_above(conc=c(6, 0, 6, 0), time=1:4, conc_above=5, method='lin up/log down'),
    1/6 + 1/6 + 1/6
  )
})

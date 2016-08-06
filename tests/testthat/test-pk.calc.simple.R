context("Simple NCA functions")

test_that("adj.r.squared", {
  ## Ensure correct calculation
  expect_equal(adj.r.squared(1, 5), 1)
  expect_equal(adj.r.squared(0.5, 5), 1-0.5*4/3)

  ## Ensure that N must be an integer > 2
  expect_error(adj.r.squared(1, 2),
               regexp="n must be > 2")
})

test_that("pk.calc.cmax", {
  ## Confirm that all NAs give NA as an output
  expect_warning(v1 <- pk.calc.cmax(NA))
  expect_equal(v1, NA)
  expect_warning(v2 <- pk.calc.cmax(c(NA, NA)))
  expect_equal(v2, NA)

  ## Confirm that no NAs give the max value
  expect_equal(pk.calc.cmax(c(1, 2)), 2)
  expect_equal(pk.calc.cmax(c(1, 2, 3)), 3)

  ## Confirm that some NAs give the NA-removed maximum value
  expect_equal(pk.calc.cmax(c(1, NA, 3)), 3)
  expect_equal(pk.calc.cmax(c(NA, NA, 3)), 3)
  expect_equal(pk.calc.cmax(c(1, NA, NA)), 1)
  expect_equal(pk.calc.cmax(c(1, NA, 3, NA)), 3)

  ## Confirm that no data gives NA with a warning
  expect_warning(v3 <- pk.calc.cmax(c()))
  expect_equal(v3, NA)
})

test_that("pk.calc.cmin", {
  ## Confirm that all NAs give NA as an output
  expect_warning(v1 <- pk.calc.cmin(NA))
  expect_equal(v1, NA)
  expect_warning(v2 <- pk.calc.cmin(c(NA, NA)))
  expect_equal(v2, NA)

  ## Confirm that no NAs give the min value
  expect_equal(pk.calc.cmin(c(1, 2)), 1)
  expect_equal(pk.calc.cmin(c(1, 2, 3)), 1)

  ## Confirm that some NAs give the NA-removed minimum value
  expect_equal(pk.calc.cmin(c(1, NA, 3)), 1)
  expect_equal(pk.calc.cmin(c(NA, NA, 3)), 3)
  expect_equal(pk.calc.cmin(c(1, NA, NA)), 1)
  expect_equal(pk.calc.cmin(c(1, NA, 3, NA)), 1)

  ## Confirm that no data gives NA with a warning
  expect_warning(v3 <- pk.calc.cmin(c()))
  expect_equal(v3, NA)
})

test_that("pk.calc.tmax", {
  ## No data give a warning and NA
  expect_warning(v1 <- pk.calc.tmax(c(), c()))
  expect_equal(v1, NA)

  ## Either concentration or time is missing, give an error
  expect_error(pk.calc.tmax(conc=c()),
               regexp="time must be given")
  expect_error(pk.calc.tmax(time=c()),
               regexp="conc must be given")

  ## It calculates tmax correctly based on the use.first option
  expect_equal(pk.calc.tmax(c(1, 2), c(0, 1), first.tmax=TRUE),
               1)
  expect_equal(pk.calc.tmax(c(1, 2), c(0, 1), first.tmax=FALSE),
               1)
  expect_equal(pk.calc.tmax(c(1, 1), c(0, 1), first.tmax=TRUE),
               0)
  expect_equal(pk.calc.tmax(c(1, 1), c(0, 1), first.tmax=FALSE),
               1)
})

test_that("pk.calc.tlast", {
  ## Either concentration or time is missing, give an error
  expect_error(pk.calc.tlast(conc=c()),
               regexp="time must be given")
  expect_error(pk.calc.tlast(time=c()),
               regexp="conc must be given")

  ## It calculates tlast correctly
  expect_equal(pk.calc.tlast(c(1, 2), c(0, 1)),
               1)
  expect_equal(pk.calc.tlast(c(0, 0), c(0, 1)),
               NA)
  expect_equal(pk.calc.tlast(c(0, 1), c(0, 1)),
               1)
  expect_equal(pk.calc.tlast(c(1, 0), c(0, 1)),
               0)
  expect_equal(pk.calc.tlast(c(1, 1), c(0, 1)),
               1)
})

test_that("pk.calc.clast.obs", {
  ## Ensure that it handles BLQ (0) values correctly
  c1 <- c(0, 1, 2, 0)
  t1 <- c(0, 1, 2, 3)
  expect_equal(pk.calc.clast.obs(c1, t1), 2)

  ## Ensure that it handles all ALQ values correctly
  c1 <- c(0, 1, 2, 3)
  t1 <- c(0, 1, 2, 3)
  expect_equal(pk.calc.clast.obs(c1, t1), 3)

  ## Ensure that it handles NA values correctly
  c1 <- c(0, 1, 2, NA)
  t1 <- c(0, 1, 2, 3)
  expect_equal(pk.calc.clast.obs(c1, t1), 2)

  c1 <- c(0, 1, NA, 3)
  t1 <- c(0, 1, 2, 3)
  expect_equal(pk.calc.clast.obs(c1, t1), 3)

  c1 <- c(NA, NA, NA, NA)
  t1 <- c(0, 1, 2, 3)
  expect_warning(v1 <- pk.calc.clast.obs(c1, t1))
  expect_equal(v1, NA)

  c1 <- rep(0, 4)
  t1 <- c(0, 1, 2, 3)
  expect_equal(pk.calc.clast.obs(c1, t1), NA)
})

test_that("pk.calc.thalf.eff", {
  ## No input gives equivalent no output
  expect_equal(pk.calc.thalf.eff(c()),
               numeric())
  
  ## NA input gives equivalent NA output
  expect_equal(pk.calc.thalf.eff(NA),
               as.numeric(NA))

  ## Numbers mixed with NA give appropriate output
  d1 <- c(0, 1, NA, 3)
  r1 <- log(2)*d1
  expect_equal(pk.calc.thalf.eff(d1),
               r1)
})

test_that("pk.calc.kel", {
  ## No input gives equivalent no output
  expect_equal(pk.calc.kel(c()),
               numeric())
  
  ## NA input gives equivalent NA output
  expect_equal(pk.calc.kel(NA),
               as.numeric(NA))

  ## Numbers mixed with NA give appropriate output
  d1 <- c(0, 1, NA, 3)
  r1 <- 1/d1
  expect_equal(pk.calc.kel(d1),
               r1)
})

test_that("pk.calc.cl", {
  ## Ensure that dose and auc are required
  expect_error(pk.calc.cl(aucinf=NA),
               info="dose is required for clearance calculation")
  expect_error(pk.calc.cl(dose=NA),
               info="aucinf is required for clearance calculation")

  expect_equal(pk.calc.cl(dose=10, aucinf=100), 0.1,
               info="Normal clearance calcualtion works")
  expect_equal(pk.calc.cl(dose=c(10, 10, 10), aucinf=c(10, NA, 100)),
               c(1, NA, 0.1),
               info="Vectors for both dose and aucinf give vectors for clearance (with NAs)")
  expect_equal(pk.calc.cl(dose=c(50, 50), aucinf=100),
               1,
               info="Vector for dose and scalar for aucinf scalar output with the sum of doses")
})

test_that("pk.calc.f", {
  expect_equal(pk.calc.f(1, 1, 1, 2), 2)
})

test_that("pk.calc.aucpext", {
  expect_equal(pk.calc.aucpext(1, 2), 50)
  expect_equal(pk.calc.aucpext(1.8, 2), 10)
  expect_warning(v1 <- pk.calc.aucpext(2, 1))
  expect_equal(v1, -100)
  expect_warning(pk.calc.aucpext(2, 1),
                 regexp="auclast should be less than aucinf")
})

test_that("pk.calc.mrt", {
  expect_equal(pk.calc.mrt(1, 2), 2)
})

test_that("pk.calc.vz", {
  ## Ensure that cl and lambda.z are required
  expect_equal(pk.calc.vz(cl=NA, lambda.z=NA), NA_integer_)
  expect_error(pk.calc.vz(cl=NA),
               info="lambda.z required for Vz calculation")
  expect_error(pk.calc.vz(lambda.z=NA),
               info="CL required for Vz calculation")

  ## Ensure that length of cl and lambda.z are either 1 or the same length
  expect_error(pk.calc.vz(cl=1:2, lambda.z=1:3),
               regexp="'cl' and 'lambda.z' must be the same length",
               info="CL and lambda.z must be the same length (CL shorter)")
  expect_error(pk.calc.vz(cl=1:3, lambda.z=1:2),
               regexp="'cl' and 'lambda.z' must be the same length",
               info="CL and lambda.z must be the same length (lambda.z shorter)")
  
  ## Estimate a single Vz (with permutations to ensure the right math
  ## is happening)
  expect_equal(pk.calc.vz(cl=1, lambda.z=1), 1,
               info="vz math test 1")
  expect_equal(pk.calc.vz(cl=1, lambda.z=2), 0.5,
               info="vz math test 2")
  expect_equal(pk.calc.vz(cl=2, lambda.z=1), 2,
               info="vz math test 3")
  expect_equal(pk.calc.vz(cl=2, lambda.z=2), 1,
               info="vz math test 4")

  ## Ensure that NA can go into either position or both
  expect_equal(pk.calc.vz(cl=NA, lambda.z=1), NA_integer_,
               info="Vz with missing (NA) cl")
  expect_equal(pk.calc.vz(cl=1, lambda.z=NA), NA_integer_,
               info="Vz with missing (NA) lambda.z")
  expect_equal(pk.calc.vz(cl=NA, lambda.z=NA), NA_integer_,
               info="Vz with missing (NA) lambda.z and cl")

  ## vectorized vz calculation works
  expect_equal(pk.calc.vz(cl=
                            c(1, 1, 1, 2, 2, 2, NA, NA, NA),
                          lambda.z=
                            c(1, 2, NA, 1, 2, NA, 1, 2, NA)),
               c(1, 0.5, NA, 2, 1, NA, NA, NA, NA),
               info="Vz with vector inputs including missing values for both parameters (cl and lambda.z)")
})

test_that("pk.calc.vss", {
  expect_equal(pk.calc.vss(1, 1), 1)
  expect_equal(pk.calc.vss(2, 1), 2)
  expect_equal(pk.calc.vss(1, 2), 2)
})

test_that("pk.calc.vd", {
  expect_equal(pk.calc.vd(1, 2, 3), 1/6,
               info="Normal Vd calculation works")
  expect_equal(pk.calc.vd(NA, 2, 3), NA_integer_,
               info="Vd calculation returns NA when dose is NA")
  expect_equal(pk.calc.vd(1, NA, 3), NA_integer_,
               info="Vd calculation returns NA when aucinf is NA")
  expect_equal(pk.calc.vd(1, 2, NA), NA_integer_,
               info="Vd calculation returns NA when lambda.z is NA")
  
  expect_equal(pk.calc.vd(c(1, 2), c(2, 4), c(3, 6)), c(1/6, 1/12),
               info="Vd calculation works with three vector inputs returning a vector")
  expect_equal(pk.calc.vd(c(1, 2), 2, 3), 0.5,
               info="Vd calculation works with vector dose and scalar aucinf and lambda.z inputs returning a scalar with the sum of doses used.")
})

test_that("pk.calc.cav", {
  expect_equal(pk.calc.cav(2, 0, 1), 2)
  expect_equal(pk.calc.cav(NA, 0, 1), as.numeric(NA))
})

test_that("pk.calc.ctrough", {
  expect_equal(pk.calc.ctrough(1:5, 0:4, 0), 1,
               info="Found and it's the first time")
  expect_equal(pk.calc.ctrough(1:5, 0:4, 1), 2,
               info="Found and it's not the first time")
  expect_equal(pk.calc.ctrough(1:5, 0:4, 1.5), NA,
               info="Not found")
  expect_error(pk.calc.ctrough(1:5, c(0, 0:3), 0),
               regexp="Time must be monotonically increasing")
})

test_that("pk.calc.ptr", {
  expect_equal(pk.calc.ptr(cmax=2, cmin=1), 2,
               info="Confirm that the ratio goes the right way")
  expect_equal(pk.calc.ptr(2, 0), as.numeric(NA),
               info="Division by zero returns NA")
})

test_that("pk.calc.tlag", {
  expect_equal(pk.calc.tlag(1:5, 0:4), 0,
               info="find the first point")
  expect_equal(pk.calc.tlag(c(0, 0, 0, 0, 1), 0:4), 3,
               info="find the next to last point")
  expect_equal(pk.calc.tlag(c(0, 0, 0, 0, 0), 0:4), NA,
               info="No increase gives NA")
  expect_equal(pk.calc.tlag(5:1, 0:4), NA,
               info="No increase gives NA")
})

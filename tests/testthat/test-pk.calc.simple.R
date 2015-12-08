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
  expect_error(pk.calc.cl(auc=NA))
  expect_error(pk.calc.cl(dose=NA))

  ## Estimate a single CL
  expect_equal(pk.calc.cl(10, 100), 0.1)
  ## Estimate multiple CLs from a single dose level (including some
  ## NA)
  expect_equal(pk.calc.cl(10, c(10, NA, 100)),
               c(1, NA, 0.1))
  ## Do unit conversion
  expect_equal(pk.calc.cl(10, c(10, NA, 100), 1000),
               c(1, NA, 0.1)*1000)
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
  ## Ensure that dose, auc, and kel are required
  expect_error(pk.calc.vz(auc=NA, kel=NA))
  expect_error(pk.calc.vz(dose=NA, kel=NA))
  expect_error(pk.calc.vz(auc=NA, dose=NA))

  ## Ensure that dose is either 1 or the same length as AUC
  expect_error(pk.calc.vz(dose=c(1, 2), auc=1:4, kel=1:4),
               regexp="'dose' and 'auc' must be the same length")
  ## Ensure that auc and kel are the same length
  expect_error(pk.calc.vz(dose=1, auc=1:4, kel=1:3),
               regexp="'auc' and 'kel' must be the same length")
  
  ## Estimate a single Vz (with permutations to ensure the right math
  ## is happening)
  expect_equal(pk.calc.vz(dose=1, auc=1, kel=1), 1)
  expect_equal(pk.calc.vz(dose=1, auc=1, kel=2), 0.5)
  expect_equal(pk.calc.vz(dose=1, auc=2, kel=1), 0.5)
  expect_equal(pk.calc.vz(dose=1, auc=2, kel=2), 0.25)
  expect_equal(pk.calc.vz(dose=2, auc=1, kel=1), 2)
  expect_equal(pk.calc.vz(dose=2, auc=1, kel=2), 1)
  expect_equal(pk.calc.vz(dose=2, auc=2, kel=1), 1)
  expect_equal(pk.calc.vz(dose=2, auc=2, kel=2), 0.5)

  ## Ensure that NA can go into any position
  expect_equal(pk.calc.vz(dose=NA, auc=1, kel=1),
               as.numeric(NA))
  expect_equal(pk.calc.vz(dose=1, auc=NA, kel=1),
               as.numeric(NA))
  expect_equal(pk.calc.vz(dose=1, auc=1, kel=NA),
               as.numeric(NA))

  ## Estimate multiple Vzs from a single dose level (including some
  ## NA)
  expect_equal(pk.calc.vz(dose=10,
                          auc=c(10, NA, 100),
                          kel=c(2, 2, 2)),
               c(0.5, NA, 0.05))
  ## Do unit conversion
  expect_equal(pk.calc.vz(dose=10,
                          auc=c(10, NA, 100),
                          kel=c(2, 2, 2),
                          unitconv=1000),
               c(0.5, NA, 0.05)*1000)
})

test_that("pk.calc.vss", {
  expect_equal(pk.calc.vss(1, 1), 1)
  expect_equal(pk.calc.vss(2, 1), 2)
  expect_equal(pk.calc.vss(1, 2), 2)
})

test_that("pk.calc.vd", {
  expect_equal(pk.calc.vd(1, 2, 3), 1/6)
  expect_equal(pk.calc.vd(NA, 2, 3), as.numeric(NA))
  expect_equal(pk.calc.vd(1, NA, 3), as.numeric(NA))
  expect_equal(pk.calc.vd(1, 2, NA), as.numeric(NA))
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

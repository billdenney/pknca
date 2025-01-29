test_that("adj.r.squared", {
  # Ensure correct calculation
  expect_equal(adj.r.squared(1, 5), 1)
  expect_equal(adj.r.squared(0.5, 5), 1-0.5*4/3)

  # Ensure that N must be an integer > 2
  expect_warning(
    expect_equal(
      adj.r.squared(1, 2),
      structure(NA_real_, exclude="n must be > 2")
    ),
    regexp="n must be > 2"
  )
})

test_that("pk.calc.cmax", {
  # Confirm that all NAs give NA as an output
  expect_warning(v1 <- pk.calc.cmax(NA))
  expect_equal(v1, NA)
  expect_warning(v2 <- pk.calc.cmax(c(NA, NA)))
  expect_equal(v2, NA)

  # Confirm that no NAs give the max value
  expect_equal(pk.calc.cmax(c(1, 2)), 2)
  expect_equal(pk.calc.cmax(c(1, 2, 3)), 3)

  # Confirm that some NAs give the NA-removed maximum value
  expect_equal(pk.calc.cmax(c(1, NA, 3)), 3)
  expect_equal(pk.calc.cmax(c(NA, NA, 3)), 3)
  expect_equal(pk.calc.cmax(c(1, NA, NA)), 1)
  expect_equal(pk.calc.cmax(c(1, NA, 3, NA)), 3)

  # Confirm that no data gives NA with a warning
  expect_warning(v3 <- pk.calc.cmax(numeric()))
  expect_equal(v3, NA)
})

test_that("pk.calc.cmin", {
  # Confirm that all NAs give NA as an output
  expect_warning(v1 <- pk.calc.cmin(NA))
  expect_equal(v1, NA)
  expect_warning(v2 <- pk.calc.cmin(c(NA, NA)))
  expect_equal(v2, NA)

  # Confirm that no NAs give the min value
  expect_equal(pk.calc.cmin(c(1, 2)), 1)
  expect_equal(pk.calc.cmin(c(1, 2, 3)), 1)

  # Confirm that some NAs give the NA-removed minimum value
  expect_equal(pk.calc.cmin(c(1, NA, 3)), 1)
  expect_equal(pk.calc.cmin(c(NA, NA, 3)), 3)
  expect_equal(pk.calc.cmin(c(1, NA, NA)), 1)
  expect_equal(pk.calc.cmin(c(1, NA, 3, NA)), 1)

  # Confirm that no data gives NA with a warning
  expect_warning(v3 <- pk.calc.cmin(numeric()))
  expect_equal(v3, NA)
})

test_that("pk.calc.tmax", {
  # No data give a warning and NA
  expect_warning(expect_warning(
    v1 <- pk.calc.tmax(numeric(), numeric()),
    class = "pknca_conc_none"),
    class = "pknca_time_none"
  )
  expect_equal(v1, NA)

  # Either concentration or time is missing, give an error
  expect_error(
    suppressWarnings(pk.calc.tmax(conc = numeric())),
    regexp='argument "time" is missing, with no default'
  )
  expect_error(
    pk.calc.tmax(time=numeric()),
    regexp='argument "conc" is missing, with no default'
  )

  # It calculates tmax correctly based on the use.first option
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
  # Either concentration or time is missing, give an error
  expect_error(
    suppressWarnings(pk.calc.tlast(conc = numeric())),
    regexp='argument "time" is missing, with no default'
  )
  expect_error(
    pk.calc.tlast(time = numeric()),
    regexp = 'argument "conc" is missing, with no default'
  )

  # It calculates tlast correctly
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

test_that("pk.calc.tfirst", {
  expect_error(
    pk.calc.tfirst(),
    regexp='argument "conc" is missing, with no default'
  )
  expect_error(
    pk.calc.tfirst(conc=1),
    regexp='argument "time" is missing, with no default'
  )
  expect_equal(
    pk.calc.tfirst(conc=1, time=2, check=TRUE),
    2
  )
})

test_that("pk.calc.clast.obs", {
  # Ensure that it handles BLQ (0) values correctly
  c1 <- c(0, 1, 2, 0)
  t1 <- c(0, 1, 2, 3)
  expect_equal(pk.calc.clast.obs(c1, t1), 2)

  # Ensure that it handles all ALQ values correctly
  c1 <- c(0, 1, 2, 3)
  t1 <- c(0, 1, 2, 3)
  expect_equal(pk.calc.clast.obs(c1, t1), 3)

  # Ensure that it handles NA values correctly
  c1 <- c(0, 1, 2, NA)
  t1 <- c(0, 1, 2, 3)
  expect_equal(pk.calc.clast.obs(c1, t1), 2)

  c1 <- c(0, 1, NA, 3)
  t1 <- c(0, 1, 2, 3)
  expect_equal(pk.calc.clast.obs(c1, t1), 3)

  c1 <- c(NA, NA, NA, NA)
  t1 <- c(0, 1, 2, 3)
  expect_warning(
    v1 <- pk.calc.clast.obs(c1, t1),
    class = "pknca_conc_all_missing"
  )
  expect_equal(v1, NA_real_)

  c1 <- rep(0, 4)
  t1 <- c(0, 1, 2, 3)
  expect_equal(pk.calc.clast.obs(c1, t1), 0)
})

test_that("pk.calc.thalf.eff", {
  # No input gives equivalent no output
  expect_equal(
    pk.calc.thalf.eff(numeric()),
    numeric()
  )

  # NA input gives equivalent NA output
  expect_equal(pk.calc.thalf.eff(NA),
               as.numeric(NA))

  # Numbers mixed with NA give appropriate output
  d1 <- c(0, 1, NA, 3)
  r1 <- log(2)*d1
  expect_equal(pk.calc.thalf.eff(d1),
               r1)
})

test_that("pk.calc.kel", {
  # No input gives equivalent no output
  expect_equal(
    pk.calc.kel(numeric()),
    numeric()
  )

  # NA input gives equivalent NA output
  expect_equal(
    pk.calc.kel(NA),
    as.numeric(NA)
  )

  # Numbers mixed with NA give appropriate output
  d1 <- c(0, 1, NA, 3)
  r1 <- 1/d1
  expect_equal(pk.calc.kel(d1), r1)
})

test_that("pk.calc.cl", {
  # Ensure that dose and auc are required
  expect_error(pk.calc.cl(auc=NA),
               info="dose is required for clearance calculation")
  expect_error(pk.calc.cl(dose=NA),
               info="auc is required for clearance calculation")

  expect_equal(pk.calc.cl(dose=10, auc=100), 0.1,
               info="Normal clearance calculation works")
  expect_equal(pk.calc.cl(dose=c(10, 10, 10), auc=c(10, NA, 100)),
               c(1, NA, 0.1),
               info="Vectors for both dose and auc give vectors for clearance (with NAs)")
  expect_equal(pk.calc.cl(dose=c(50, 50), auc=100),
               1,
               info="Vector for dose and scalar for auc scalar output with the sum of doses")
  expect_equal(pk.calc.cl(dose=c(50, 50), auc=c(NA, 100)),
               c(NA_real_, 0.5),
               info="NA generation works with NA, including for a vector")
  expect_equal(pk.calc.cl(dose=c(50, 50), auc=c(0, 100)),
               c(NA_real_, 0.5),
               info="NA generation works with zero, including for a vector")
})

test_that("pk.calc.f", {
  expect_equal(pk.calc.f(1, 1, 1, 2), 2,
               info="Standard bioavailability calculation")
  expect_equal(pk.calc.f(NA, 1, 1, 1), NA_real_,
               info="NA handling per parameter in bioavailability (1)")
  expect_equal(pk.calc.f(1, NA, 1, 1), NA_real_,
               info="NA handling per parameter in bioavailability (2)")
  expect_equal(pk.calc.f(1, 1, NA, 1), NA_real_,
               info="NA handling per parameter in bioavailability (3)")
  expect_equal(pk.calc.f(1, 1, 1, NA), NA_real_,
               info="NA handling per parameter in bioavailability (4)")
  expect_equal(pk.calc.f(0, 1, 1, 1), NA_real_,
               info="Zero handling per parameter in bioavailability (1)")
  expect_equal(pk.calc.f(1, 0, 1, 1), NA_real_,
               info="Zero handling per parameter in bioavailability (2)")
  expect_equal(pk.calc.f(1, 1, 0, 1), NA_real_,
               info="Zero handling per parameter in bioavailability (3)")
  expect_equal(pk.calc.f(1, 1, 1, 0), 0,
               info="Zero handling per parameter in bioavailability (4)")
})

test_that("pk.calc.aucpext", {
  expect_equal(pk.calc.aucpext(1, 2), 50)
  expect_equal(pk.calc.aucpext(1.8, 2), 10)
  expect_warning(v1 <- pk.calc.aucpext(2, 1),
                 regexp="aucpext is typically only calculated when aucinf is greater than auclast.")
  expect_equal(v1, -100)
  expect_warning(expect_warning(
    v2 <- pk.calc.aucpext(auclast=0, aucinf=0),
    class = "pknca_aucpext_aucinf_le_auclast"),
    class = "pknca_aucpext_aucinf_auclast_positive"
  )
  expect_equal(v2, NA_real_,
               info="aucinf<=0 gives NA_real_ (not infinity)")
  expect_equal(pk.calc.aucpext(NA, NA),
               NA_real_,
               info="Percent extrapolated is NA when both inputs are NA.")
  expect_equal(pk.calc.aucpext(NA, 1),
               NA_real_,
               info="Percent extrapolated is NA when input auclast is NA.")
  expect_equal(pk.calc.aucpext(1, NA),
               NA_real_,
               info="Percent extrapolated is NA when input aucinf is NA.")
  expect_error(pk.calc.aucpext(1:2, 1:3),
               regexp="auclast and aucinf must either be a scalar or the same length.",
               info="aucpext input length checks require consistency.")
  expect_equal(pk.calc.aucpext(c(1, NA), 2),
               c(50, NA_real_),
               info="Percent extrapolated works with vector/scalar input.")
  expect_equal(pk.calc.aucpext(1, c(2, NA)),
               c(50, NA_real_),
               info="Percent extrapolated works with scalar/vector input.")
  expect_equal(pk.calc.aucpext(c(1, 1), c(2, NA)),
               c(50, NA_real_),
               info="Percent extrapolated works with vector/vector input.")
})

test_that("pk.calc.mrt", {
  expect_equal(pk.calc.mrt(auc=1, aumc=2),
               2,
               info="MRT is calculated correctly")
})

test_that("pk.calc.mrt.iv", {
  expect_equal(pk.calc.mrt.iv(auc=1, aumc=2, duration.dose=1),
               1.5,
               info="MRT.iv is calculated correctly")
  expect_equal(pk.calc.mrt.iv(auc=1, aumc=2, duration.dose=0),
               2,
               info="MRT.iv is calculated correctly when duration is 0")
  expect_equal(pk.calc.mrt.iv(auc=1, aumc=2, duration.dose=NA),
               NA_real_,
               info="MRT.iv is calculated correctly when duration is missing")
  expect_equal(pk.calc.mrt.iv(auc=0, aumc=2, duration.dose=NA),
               NA_real_,
               info="MRT.iv is calculated correctly when auc is zero")
})

test_that("pk.calc.mrt.md", {
  expect_equal(pk.calc.mrt.md(1, 2, 1.5, 24), 2 + 24*0.5)
  expect_equal(pk.calc.mrt.md(0, 2, 1.5, 24), NA_real_,
               info="auctau <= 0 becomes NA (not Inf)")
})

test_that("pk.calc.vz", {
  # Ensure that cl and lambda.z are required
  expect_equal(pk.calc.vz(cl=NA, lambda.z=NA), NA_integer_)
  expect_error(pk.calc.vz(cl=NA),
               info="lambda.z required for Vz calculation")
  expect_error(pk.calc.vz(lambda.z=NA),
               info="CL required for Vz calculation")

  # Ensure that length of cl and lambda.z are either 1 or the same length
  expect_error(
    pk.calc.vz(cl=1:2, lambda.z=1:3),
    regexp="'cl' and 'lambda.z' must be the same length",
    info="CL and lambda.z must be the same length (CL shorter)"
  )
  expect_error(
    pk.calc.vz(cl=1:3, lambda.z=1:2),
    regexp="'cl' and 'lambda.z' must be the same length",
    info="CL and lambda.z must be the same length (lambda.z shorter)"
  )

  # Estimate a single Vz (with permutations to ensure the right math
  # is happening)
  expect_equal(pk.calc.vz(cl=1, lambda.z=1), 1,
               info="vz math test 1")
  expect_equal(pk.calc.vz(cl=1, lambda.z=2), 0.5,
               info="vz math test 2")
  expect_equal(pk.calc.vz(cl=2, lambda.z=1), 2,
               info="vz math test 3")
  expect_equal(pk.calc.vz(cl=2, lambda.z=2), 1,
               info="vz math test 4")

  # Ensure that NA can go into either position or both
  expect_equal(pk.calc.vz(cl=NA, lambda.z=1), NA_integer_,
               info="Vz with missing (NA) cl")
  expect_equal(pk.calc.vz(cl=1, lambda.z=NA), NA_integer_,
               info="Vz with missing (NA) lambda.z")
  expect_equal(pk.calc.vz(cl=NA, lambda.z=NA), NA_integer_,
               info="Vz with missing (NA) lambda.z and cl")

  # vectorized vz calculation works
  expect_equal(pk.calc.vz(cl=
                            c(1, 1, 1, 2, 2, 2, NA, NA, NA),
                          lambda.z=
                            c(1, 2, NA, 1, 2, NA, 1, 2, NA)),
               c(1, 0.5, NA, 2, 1, NA, NA, NA, NA),
               info="Vz with vector inputs including missing values for both parameters (cl and lambda.z)")
})

test_that("pk.calc.vss and its wrappers", {
  expect_equal(pk.calc.vss(1, 1), 1)
  expect_equal(pk.calc.vss(2, 1), 2)
  expect_equal(pk.calc.vss(1, 2), 2)
})

test_that("pk.calc.cav", {
  expect_equal(pk.calc.cav(2, 0, 1), 2)
  expect_equal(pk.calc.cav(NA, 0, 1), NA_real_)
  expect_equal(pk.calc.cav(2, 1, 1), NA_real_,
               info="If end == start, return NA_real_")
})

test_that("pk.calc.ctrough", {
  expect_equal(pk.calc.ctrough(1:5, 0:4, 0), 1,
               info="Found and it's the first time")
  expect_equal(pk.calc.ctrough(1:5, 0:4, 1), 2,
               info="Found and it's not the first time")
  expect_equal(pk.calc.ctrough(1:5, 0:4, 1.5), NA_real_,
               info="Not found")
  expect_error(
    pk.calc.ctrough(1:5, c(0, 0:3), 0),
    regexp = "Assertion on 'time' failed: Contains duplicated values, position 2."
  )
})

test_that("pk.calc.ptr", {
  expect_equal(pk.calc.ptr(cmax=2, ctrough=1), 2,
               info="Confirm that the ratio goes the right way")
  expect_equal(pk.calc.ptr(2, 0), NA_real_,
               info="Division by zero returns NA")
})

test_that("pk.calc.tlag", {
  expect_equal(pk.calc.tlag(1:5, 0:4), 0,
               info="find the first point")
  expect_equal(pk.calc.tlag(c(0, 0, 0, 0, 1), 0:4), 3,
               info="find the next to last point")
  expect_equal(pk.calc.tlag(c(0, 0, 0, 0, 0), 0:4), NA_real_,
               info="No increase gives NA")
  expect_equal(pk.calc.tlag(5:1, 0:4), NA_real_,
               info="No increase gives NA")
})

test_that("pk.calc.deg.fluc", {
  expect_equal(pk.calc.deg.fluc(cmax=100, cmin=10, cav=45), 200,
               info="Degree of fluctuation math works")
  expect_equal(pk.calc.deg.fluc(cmax=100, cmin=10, cav=0), NA_real_,
               info="Degree of fluctuation returns NA when cav=0")
})

test_that("pk.calc.swing", {
  expect_equal(pk.calc.swing(100, 10), 900,
               info="Swing math works")
  expect_equal(pk.calc.swing(100, 0), Inf,
               info="Swing handle Ctrough=0")
  expect_equal(pk.calc.swing(100, -1), Inf,
               info="Swing handle Ctrough<0")
})

test_that("pk.calc.ceoi", {
  expect_equal(pk.calc.ceoi(conc=0:5, time=0:5, duration.dose=1),
               1,
               info="Ceoi returns the concentration at the end of the dosing duration")
  expect_equal(pk.calc.ceoi(conc=0:5, time=0:5, duration.dose=1.5),
               NA_real_,
               info="Ceoi returns NA if there is no measurement at the end of the dosing duration")
  expect_equal(pk.calc.ceoi(conc=0:5, time=0:5, duration.dose=NA),
               NA_real_,
               info="Ceoi returns NA if there is no dosing duration")
})

test_that("pk.calc.aucabove", {
  expect_equal(
    pk.calc.aucabove(conc = c(0:5, 1), time = 0:6, conc_above = 2),
    pk.calc.auc.all(conc = c(0, 0, 0, 1:3, 0), time = 0:6)
  )

  expect_equal(
  pk.calc.aucabove(conc = c(0:5, 1), time = 0:6, conc_above = NA_real_),
  structure(NA_real_, exclude = "Missing concentration to be above")
  )

  # Confirm that it works through NCA calculations
  d_conc <- data.frame(conc = c(2, 1:5, 3), time = 0:6)
  d_intervals <- data.frame(start = 0, end = 6, aucabove.trough.all = TRUE, aucabove.predose.all = TRUE)
  o_conc <- PKNCAconc(d_conc, conc~time)
  o_data <- PKNCAdata(o_conc, intervals = d_intervals)
  suppressMessages(
    o_nca <- pk.nca(o_data)
  )
  expect_equal(
    as.data.frame(o_nca),
    tibble::tibble(
      start = 0, end = 6,
      PPTESTCD = c("ctrough", "cstart", "aucabove.predose.all", "aucabove.trough.all"),
      PPORRES =
        c(
          3, 2,
          pk.calc.aucabove(conc = d_conc$conc, time = d_conc$time, conc_above = 2),
          pk.calc.aucabove(conc = d_conc$conc, time = d_conc$time, conc_above = 3)
        ),
      exclude = NA_character_
    )
  )
})

test_that("pk.calc.count_conc", {
  expect_equal(pk.calc.count_conc(1:5), 5)
  expect_equal(pk.calc.count_conc(c(1:2, NA)), 2)
  expect_equal(pk.calc.count_conc(c(1:2, NA, 0)), 3)
  expect_equal(suppressWarnings(pk.calc.count_conc(numeric())), 0)
  expect_equal(suppressWarnings(pk.calc.count_conc(NA)), 0)
})

test_that("pk.calc.count_conc_measured", {
  expect_equal(pk.calc.count_conc_measured(1:5), 5)
  expect_equal(pk.calc.count_conc_measured(c(1:2, NA)), 2)
  # including BLQ
  expect_equal(pk.calc.count_conc_measured(c(1:2, NA, 0)), 2)
  # Other
  expect_equal(suppressWarnings(pk.calc.count_conc_measured(numeric())), 0)
  expect_equal(suppressWarnings(pk.calc.count_conc_measured(NA)), 0)
})

test_that("pk.calc.totdose", {
  expect_equal(pk.calc.totdose(1), 1)
  expect_equal(pk.calc.totdose(c(1, 1)), 2)
})

test_that("pk.calc.cstart", {
  expect_equal(pk.calc.cstart(1:5, 0:4, 0), 1)
  expect_equal(pk.calc.cstart(1:5, 0:4, 1), 2)
  expect_equal(pk.calc.cstart(1:5, 0:4, 1.5), NA_real_)
  expect_error(
    pk.calc.cstart(1:5, c(0, 0:3), 0),
    regexp = "Assertion on 'time' failed: Contains duplicated values, position 2."
  )
})

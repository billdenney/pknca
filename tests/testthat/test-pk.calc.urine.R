test_that("pk.calc.ae", {
  expect_equal(
    pk.calc.ae(conc=1:5, volume=1:5),
    sum((1:5)^2)
  )
  expect_equal(
    pk.calc.ae(conc = NA, volume = NA),
    structure(NA_real_, exclude = "All concentrations and volumes are missing")
  )
  expect_equal(
    pk.calc.ae(conc = c(1, NA), volume = c(1, NA)),
    structure(NA_real_, exclude = "1 of 2 concentrations and volumes are missing")
  )
  expect_equal(
    pk.calc.ae(conc = c(NA, NA), volume = c(1, 1)),
    structure(NA_real_, exclude = "All concentrations are missing")
  )
  expect_equal(
    pk.calc.ae(conc = c(1, NA), volume = c(1, 1)),
    structure(NA_real_, exclude = "1 of 2 concentrations are missing")
  )
  expect_equal(
    pk.calc.ae(conc = c(1, 1), volume = c(NA, NA)),
    structure(NA_real_, exclude = "All volumes are missing")
  )
  expect_equal(
    pk.calc.ae(conc = c(1, 1), volume = c(NA, 1)),
    structure(NA_real_, exclude = "1 of 2 volumes are missing")
  )
  expect_equal(
    pk.calc.ae(conc = c(NA, NA, 1, 1), volume = c(NA, 1, NA, 1)),
    structure(NA_real_, exclude = "1 of 4 concentrations and volumes are missing; 1 of 4 concentrations are missing; 1 of 4 volumes are missing")
  )
})

test_that("pk.calc.clr", {
  expect_equal(pk.calc.clr(ae=1, auc=10),
               0.1,
               info="CLr is calculated correctly with both scalars")
  expect_equal(pk.calc.clr(ae=c(1, 2), auc=10),
               0.3,
               info="CLr is calculated correctly with a vector Ae and a scalar AUC")
  expect_equal(pk.calc.clr(ae=c(1, 2), auc=c(1, 10)),
               c(3, 0.3),
               info="CLr is calculated correctly with both vectors (but that is not the likely calculation method)")
})

test_that("pk.calc.fe", {
  expect_equal(pk.calc.fe(1, 10),
               0.1,
               info="fe is calculated correctly with both scalars")
  expect_equal(pk.calc.fe(c(1, 2), 10),
               0.3,
               info="fe is calculated correctly with both vector/scalar")
})

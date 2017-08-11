context("Urine")

test_that("pk.calc.ae", {
  expect_equal(pk.calc.ae(conc=1:5,
                          volume=1:5),
               sum((1:5)^2),
               info="Ae is the sum of concentration and volume")
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
test_that("pk.calc.auciv", {
  expect_equal(
    pk.calc.auciv(conc = 1:2, time = 1:2),
    structure(NA_real_, exclude="No time 0 in data")
  )
  expect_equal(
    # No check is done to confirm that the auc argument matches the data
    pk.calc.auciv(conc = 0:5, time = 0:5, c0 = 1, auc = 2.75),
    2.75 + 1 - 0.5
  )
})

test_that("pk.calc.auciv_pbext", {
  expect_equal(
    pk.calc.auciv_pbext(auc = 1, auciv = 2.1),
    100 * (1 - 1/2.1)
  )
})

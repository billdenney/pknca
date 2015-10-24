context("Half-life")

test_that("pk.calc.half.life", {
  ## Confirm that half-life is correctly calculated with a simple
  ## exponential decay
  expect_equal(pk.calc.half.life(conc=c(1, 0.5, 0.25),
                                 time=c(0, 1, 2),
                                 min.hl.points=3,
                                 allow.tmax.in.half.life=TRUE,
                                 adj.r.squared.factor=0.0001)$half.life,
               1)

  ## Ensure that when input data is not checked, the code works
  ## correctly.
  expect_equal(pk.calc.half.life(conc=c(1, 0.5, 0.25),
                                 time=c(0, 1, 2),
                                 min.hl.points=3,
                                 allow.tmax.in.half.life=TRUE,
                                 adj.r.squared.factor=0.0001,
                                 check=FALSE)$half.life,
               1)

  ## Ensure that min.hl.points is respected
  expect_warning(pk.calc.half.life(conc=c(1, 0.5, 0.25),
                                   time=c(0, 1, 2),
                                   min.hl.points=4),
                 regexp="Too few points for half-life calculation")

  ## Ensure that when there are more than one best models by adjusted
  ## r-squared, the one with the most points is used.
  expect_warning(
    expect_equal(pk.calc.half.life(conc=c(1, 0.5, 0.25, 0.1251),
                                   time=c(0, 1, 2, 3),
                                   min.hl.points=3,
                                   allow.tmax.in.half.life=TRUE,
                                   adj.r.squared.factor=0.1,
                                   check=FALSE)$half.life,
                 1.000346,
                 tolerance=0.0001))

  ## Ensure that the allow.tmax.in.half.life parameter is followed
  expect_warning(
    expect_equal(pk.calc.half.life(conc=c(1, 0.5, 0.25, 0.1251),
                                   time=c(0, 1, 2, 3),
                                   min.hl.points=3,
                                   adj.r.squared.factor=0.1,
                                   check=FALSE)$half.life,
                 1.000577,
                 tolerance=0.0001))

})

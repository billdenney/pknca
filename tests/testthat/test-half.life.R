test_that("pk.calc.half.life", {
  ## Confirm that half-life is correctly calculated with a simple
  ## exponential decay
  expect_equal(pk.calc.half.life(conc=c(1, 0.5, 0.25),
                                 time=c(0, 1, 2),
                                 min.hl.points=3,
                                 adj.r.squared.factor=0.0001)$half.life,
               1)

  ## Ensure that min.hl.points is respected
  expect_warning(pk.calc.half.life(conc=c(1, 0.5, 0.25),
                                   time=c(0, 1, 2),
                                   min.hl.points=4),
                 regexp="Too few points for half-life calculation \\(min.hl.points=4 with only 3 points\\)")
})

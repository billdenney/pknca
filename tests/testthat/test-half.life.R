context("Half-life")

test_that("pk.calc.half.life", {
  ## Confirm that half-life is correctly calculated with a simple
  ## exponential decay
  expect_warning(v1 <-
    pk.calc.half.life(conc=c(1, 0.5, 0.25),
                      time=c(0, 1, 2),
                      min.hl.points=3,
                      allow.tmax.in.half.life=TRUE,
                      adj.r.squared.factor=0.0001)$half.life,
                 regexp="essentially perfect fit: summary may be unreliable")
  expect_equal(v1, 1)

  ## Ensure that when input data is not checked, the code works
  ## correctly.
  expect_warning(v2 <-
    pk.calc.half.life(conc=c(1, 0.5, 0.25),
                      time=c(0, 1, 2),
                      min.hl.points=3,
                      allow.tmax.in.half.life=TRUE,
                      adj.r.squared.factor=0.0001,
                      check=FALSE)$half.life,
                 regexp="essentially perfect fit: summary may be unreliable")
                 
  expect_equal(v2, 1)

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
  expect_equal(pk.calc.half.life(conc=c(1, 0.5, 0.25, 0.1251),
                                 time=c(0, 1, 2, 3),
                                 min.hl.points=3,
                                 allow.tmax.in.half.life=FALSE,
                                 check=FALSE)$half.life,
               1.000577,
               tolerance=0.00001)
  expect_warning(
    expect_equal(pk.calc.half.life(conc=c(1, 0.5, 0.25, 0.1251),
                                   time=c(0, 1, 2, 3),
                                   min.hl.points=3,
                                   allow.tmax.in.half.life=TRUE,
                                   adj.r.squared.factor=0.1,
                                   check=FALSE)$half.life,
                 1.000346,
                 tolerance=0.00001))
  

  ## Make sure that when tmax and tlast are given that they are
  ## automatically included and not recalculated
  expect_equal(pk.calc.half.life(conc=c(1, 0.5, 0.25, 0.1251),
                                 time=c(0, 1, 2, 3),
                                 min.hl.points=3,
                                 tmax=0,
                                 tlast=3,
                                 allow.tmax.in.half.life=TRUE,
                                 adj.r.squared.factor=0.01,
                                 check=FALSE),
               data.frame(lambda.z=0.6929073,
                          r.squared=0.9999999,
                          adj.r.squared=0.9999999,
                          lambda.z.time.first=0,
                          lambda.z.n.points=4,
                          clast.pred=0.12507,
                          half.life=1.000346,
                          span.ratio=2.998962),
               tolerance=0.0001)
  ## It only gives tlast or tmax if you don't give them as inputs.
  expect_equal(pk.calc.half.life(conc=c(1, 0.5, 0.25, 0.1251),
                                 time=c(0, 1, 2, 3),
                                 min.hl.points=3,
                                 tlast=3,
                                 allow.tmax.in.half.life=TRUE,
                                 adj.r.squared.factor=0.01,
                                 check=FALSE),
               data.frame(lambda.z=0.6929073,
                          r.squared=0.9999999,
                          adj.r.squared=0.9999999,
                          lambda.z.time.first=0,
                          lambda.z.n.points=4,
                          clast.pred=0.12507,
                          half.life=1.000346,
                          span.ratio=2.998962,
                          tmax=0),
               tolerance=0.0001)
    expect_equal(pk.calc.half.life(conc=c(1, 0.5, 0.25, 0.1251),
                                 time=c(0, 1, 2, 3),
                                 min.hl.points=3,
                                 tmax=0,
                                 allow.tmax.in.half.life=TRUE,
                                 adj.r.squared.factor=0.01,
                                 check=FALSE),
               data.frame(lambda.z=0.6929073,
                          r.squared=0.9999999,
                          adj.r.squared=0.9999999,
                          lambda.z.time.first=0,
                          lambda.z.n.points=4,
                          clast.pred=0.12507,
                          half.life=1.000346,
                          span.ratio=2.998962,
                          tlast=3),
               tolerance=0.0001)

})

test_that("half-life manual point selection", {
  expect_equal(
      pk.calc.half.life(conc=c(3, 1, 0.5, 0.13, 0.12, 0.113),
                        time=c(0, 1, 2, 3, 4, 5),
                        manually.selected.points=TRUE,
                        min.hl.points=3,
                        allow.tmax.in.half.life=FALSE,
                        check=FALSE)$half.life,
      1.00653,
      tolerance=0.0001,
      info="manual selection uses the given points as is")
  expect_true(
    pk.calc.half.life(conc=c(3, 1, 0.5, 0.13, 0.12, 0.113),
                      time=c(0, 1, 2, 3, 4, 5),
                      manually.selected.points=TRUE,
                      min.hl.points=3,
                      allow.tmax.in.half.life=FALSE,
                      check=FALSE)$half.life !=
      pk.calc.half.life(conc=c(3, 1, 0.5, 0.13, 0.12, 0.113),
                        time=c(0, 1, 2, 3, 4, 5),
                        min.hl.points=3,
                        allow.tmax.in.half.life=FALSE,
                        check=FALSE)$half.life,
    info="manual selection is different than automatic selection")
  expect_true(
    pk.calc.half.life(conc=c(3, 1, 0.5, 0.13, 0.12, 0.113),
                      time=c(0, 1, 2, 3, 4, 5),
                      manually.selected.points=TRUE,
                      min.hl.points=3,
                      tlast=20,
                      allow.tmax.in.half.life=FALSE,
                      check=FALSE)$clast.pred <
      pk.calc.half.life(conc=c(3, 1, 0.5, 0.13, 0.12, 0.113),
                        time=c(0, 1, 2, 3, 4, 5),
                        manually.selected.points=TRUE,
                        min.hl.points=3,
                        allow.tmax.in.half.life=FALSE,
                        check=FALSE)$clast.pred,
    info="manually-selected half-life respects tlast and generates a different clast.pred")
  expect_warning(manual_blq <-
                   pk.calc.half.life(conc=rep(0, 6),
                                     time=c(0, 1, 2, 3, 4, 5),
                                     manually.selected.points=TRUE,
                                     min.hl.points=3,
                                     tlast=20,
                                     allow.tmax.in.half.life=FALSE,
                                     check=FALSE),
                 regexp="No data to manually fit for half-life (all concentrations may be 0)",
                 fixed=TRUE,
                 info="All BLQ with manual point selection gives a warning")
  expect_true(all(is.na(unlist(manual_blq))),
              info="All BLQ with manual point selection gives all NA results")
})
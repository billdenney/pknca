test_that("pk.calc.half.life", {
  # Confirm that half-life is correctly calculated with a simple
  # exponential decay
  v1 <-
    pk.calc.half.life(conc=c(1, 0.5, 0.25),
                      time=c(0, 1, 2),
                      min.hl.points=3,
                      allow.tmax.in.half.life=TRUE,
                      adj.r.squared.factor=0.0001)$half.life
  expect_equal(v1, 1)

  # Ensure that when input data is not checked, the code works
  # correctly.
  v2 <-
    pk.calc.half.life(conc=c(1, 0.5, 0.25),
                      time=c(0, 1, 2),
                      min.hl.points=3,
                      allow.tmax.in.half.life=TRUE,
                      adj.r.squared.factor=0.0001,
                      check=FALSE)$half.life

  expect_equal(v2, 1)

  # Ensure that min.hl.points is respected
  expect_warning(pk.calc.half.life(conc=c(1, 0.5, 0.25),
                                   time=c(0, 1, 2),
                                   min.hl.points=4),
                 regexp="Too few points for half-life calculation")

  # Ensure that when there are more than one best models by adjusted
  # r-squared, the one with the most points is used.
  expect_warning(
    expect_equal(pk.calc.half.life(conc=c(1, 0.5, 0.25, 0.1251),
                                   time=c(0, 1, 2, 3),
                                   min.hl.points=3,
                                   allow.tmax.in.half.life=TRUE,
                                   adj.r.squared.factor=0.1,
                                   check=FALSE)$half.life,
                 1.000346,
                 tolerance=0.0001))

  # Ensure that the allow.tmax.in.half.life parameter is followed
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


  # Make sure that when tmax and tlast are given that they are
  # automatically included and not recalculated
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
  # It only gives tlast or tmax if you don't give them as inputs.
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
  expect_warning(
    manual_blq <-
      pk.calc.half.life(
        conc=rep(0, 6),
        time=c(0, 1, 2, 3, 4, 5),
        manually.selected.points=TRUE,
        min.hl.points=3,
        tlast=20,
        allow.tmax.in.half.life=FALSE,
        check=FALSE
      ),
    regexp="No data to manually fit for half-life (all concentrations may be 0 or excluded)",
    fixed=TRUE,
    info="All BLQ with manual point selection gives a warning"
  )
  expect_true(all(is.na(unlist(manual_blq))),
              info="All BLQ with manual point selection gives all NA results")

  excluded_result <-
    data.frame(
      lambda.z = -log(2),
      r.squared = 1,
      adj.r.squared = 1,
      lambda.z.time.first = 1L,
      lambda.z.n.points = 3L,
      clast.pred = 8,
      half.life = -1,
      span.ratio = -2,
      tmax = 3L,
      tlast = 3L
    )
  attr(excluded_result, "exclude") <- "Negative half-life estimated with manually-selected points"
  expect_equal(
    pk.calc.half.life(conc = 2^(1:3), time = 1:3, manually.selected.points = TRUE),
    excluded_result
  )
})

test_that("two-point half-life succeeds (fix #114)", {
  expect_warning(expect_warning(
    expect_equal(
      pk.calc.half.life(
        conc=c(1, 0.5),
        time=c(0, 1),
        min.hl.points=2,
        allow.tmax.in.half.life=TRUE,
        check=FALSE
      ),
      data.frame(
        lambda.z=log(2),
        r.squared=1,
        adj.r.squared=NA_real_,
        lambda.z.time.first=0,
        lambda.z.n.points=2,
        clast.pred=0.5,
        half.life=1,
        span.ratio=1,
        tmax=0,
        tlast=1
      )
    ),
    class = "pknca_halflife_2points"),
    class = "pknca_adjr2_2points"
  )
})

test_that("get_halflife_points", {
  o_conc <- PKNCAconc(Theoph, conc~Time|Subject)
  o_data <- PKNCAdata(o_conc, intervals = data.frame(start = 0, end = Inf, half.life = TRUE))
  o_nca <- pk.nca(o_data)
  hl_points <- suppressMessages(get_halflife_points(o_nca))

  # Visualize the results
  # d_all <- Theoph
  # d_all$hl_points <- hl_points
  # ggplot(d_all, aes(x = Time, y = conc, colour = hl_points, groups = Subject)) +
  #   geom_point() + geom_line()
  expect_equal(
    hl_points,
    c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE,
      TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE,
      FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE,
      FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE,
      TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)
  )

  # Setup for the remaining tests
  d_conc <- as.data.frame(Theoph[Theoph$Subject %in% Theoph$Subject[1], ])
  d_conc$incl <- c(FALSE, TRUE, rep(NA, 8), TRUE)
  d_conc$excl_txt <- c(NA, "x", rep(NA, 8), "x")
  d_conc$excl <- c(NA, TRUE, rep(NA, 6), TRUE, FALSE, TRUE)

  # No modification
  o_conc <- PKNCAconc(d_conc, conc~Time|Subject)
  o_data <- PKNCAdata(o_conc, intervals = data.frame(start = 0, end = Inf, half.life = TRUE))
  o_nca <- suppressMessages(suppressWarnings(pk.nca(o_data)))
  expect_equal(
    suppressMessages(get_halflife_points(o_nca)),
    c(rep(FALSE, 8), rep(TRUE, 3))
  )

  # Test include_half.life
  o_conc <- PKNCAconc(d_conc, conc~Time|Subject, include_half.life = "incl")
  o_data <- PKNCAdata(o_conc, intervals = data.frame(start = 0, end = Inf, half.life = TRUE))
  o_nca <- suppressMessages(suppressWarnings(pk.nca(o_data)))
  expect_equal(
    suppressMessages(get_halflife_points(o_nca)),
    c(FALSE, TRUE, rep(FALSE, 8), TRUE)
  )

  # Test exclude
  o_conc <- PKNCAconc(d_conc, conc~Time|Subject, exclude = "excl_txt")
  o_data <- PKNCAdata(o_conc, intervals = data.frame(start = 0, end = Inf, half.life = TRUE))
  o_nca <- suppressMessages(pk.nca(o_data))
  expect_equal(
    suppressMessages(get_halflife_points(o_nca)),
    # NA values indicate that the point was excluded from all calculations
    c(FALSE, NA, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, NA)
  )

  # Test exclude_half.life
  o_conc <- PKNCAconc(d_conc, conc~Time|Subject, exclude_half.life = "excl")
  o_data <- PKNCAdata(o_conc, intervals = data.frame(start = 0, end = Inf, half.life = TRUE))
  o_nca <- suppressMessages(pk.nca(o_data))
  expect_equal(
    suppressMessages(get_halflife_points(o_nca)),
    # NA values indicate that the point was excluded from all calculations
    c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE)
  )

  # Multiple, overlapping half-life calculations gives an error
  o_conc <- PKNCAconc(d_conc, conc~Time|Subject)
  o_data <-
    PKNCAdata(
      o_conc,
      intervals =
        data.frame(
          start = 0, end = c(24, Inf),
          half.life = TRUE
        )
    )
  o_nca <- suppressMessages(suppressWarnings(pk.nca(o_data)))
  expect_error(
    suppressMessages(get_halflife_points(o_nca)),
    regexp = "More than one half-life calculation was attempted on the following rows: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11"
  )
})

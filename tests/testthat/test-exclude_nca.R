test_that("exclude_nca", {
  my_conc <- PKNCAconc(data.frame(conc=c(1.1^(3:0), 1.1), time=0:4, subject=1), conc~time|subject)
  my_data <- PKNCAdata(my_conc, intervals=data.frame(start=0, end=Inf, aucinf.obs=TRUE, aucpext.obs=TRUE))
  suppressMessages(
    my_result <- pk.nca(my_data)
  )

  my_result_excluded <- exclude(my_result, FUN=exclude_nca_max.aucinf.pext())
  expect_equal(as.data.frame(my_result_excluded)$exclude,
               c(rep(NA_character_, nrow(my_result_excluded$result)-2),
                 rep("AUC percent extrapolated > 20", 2)))

  my_result_excluded <- exclude(my_result, FUN=exclude_nca_max.aucinf.pext(50))
  expect_equal(as.data.frame(my_result_excluded)$exclude,
               c(rep(NA_character_, nrow(my_result_excluded$result)-2),
                 rep("AUC percent extrapolated > 50", 2)))

  my_result_excluded <- exclude(my_result, FUN=exclude_nca_span.ratio())
  expect_equal(as.data.frame(my_result_excluded)$exclude,
               c(rep(NA_character_, 4),
                 rep("Span ratio < 2", 11)))
  my_result_excluded <- exclude(my_result, FUN=exclude_nca_span.ratio(1))
  expect_equal(as.data.frame(my_result_excluded)$exclude,
               c(rep(NA_character_, 4),
                 rep("Span ratio < 1", 11)))

  my_result_excluded <- exclude(my_result, FUN=exclude_nca_min.hl.r.squared())
  expect_equal(as.data.frame(my_result_excluded)$exclude,
               c(rep(NA_character_, 4),
                 rep("Half-life r-squared < 0.9", 11)))
  my_result_excluded <- exclude(my_result, FUN=exclude_nca_min.hl.r.squared(0.95))
  expect_equal(as.data.frame(my_result_excluded)$exclude,
               c(rep(NA_character_, 4),
                 rep("Half-life r-squared < 0.95", 11)))

  my_data <- PKNCAdata(my_conc, intervals=data.frame(start=0, end=Inf, cmax=TRUE))
  suppressMessages(
    my_result <- pk.nca(my_data)
  )
  expect_equal(my_result,
               exclude(my_result, FUN=exclude_nca_max.aucinf.pext()),
               info="Result is ignored when not calculated")
  expect_equal(my_result,
               exclude(my_result, FUN=exclude_nca_span.ratio()),
               info="Result is ignored when not calculated")
  expect_equal(my_result,
               exclude(my_result, FUN=exclude_nca_min.hl.r.squared()),
               info="Result is ignored when not calculated")
})

test_that("exclude_nca_count_conc_measured", {
  my_conc <- PKNCAconc(data.frame(conc=c(1.1^(c(3:0, -Inf)), 1.1), time=0:5, subject = 1), conc~time|subject)
  my_data <- PKNCAdata(my_conc, intervals=data.frame(start=0, end=Inf, aucinf.obs=TRUE, aucpext.obs=TRUE, count_conc_measured = TRUE))
  suppressMessages(
    my_result <- pk.nca(my_data)
  )
  expect_equal(
    as.data.frame(my_result)$exclude,
    rep(NA_character_, 16)
  )
  my_result_exclude5 <- exclude(my_result, FUN = exclude_nca_count_conc_measured(min_count = 5))
  expect_equal(
    as.data.frame(my_result_exclude5)$exclude,
    rep(NA_character_, 16)
  )
  my_result_exclude10 <- exclude(my_result, FUN = exclude_nca_count_conc_measured(min_count = 10))
  expect_equal(
    as.data.frame(my_result_exclude10)$exclude,
    c("Number of measured concentrations is < 10", rep(NA_character_, 13), rep("Number of measured concentrations is < 10", 2))
  )
})

test_that("exclude_nca_tmax_early", {
  my_conc <- PKNCAconc(data.frame(conc=c(1.1^(c(3:0, -Inf)), 1.1), time=0:5, subject = 1), conc~time|subject)
  my_data <- PKNCAdata(my_conc, intervals=data.frame(start=0, end=Inf, auclast = TRUE, half.life = TRUE))
  suppressMessages(
    my_result <- pk.nca(my_data)
  )
  expect_equal(
    as.data.frame(my_result)$exclude,
    rep(NA_character_, nrow(as.data.frame(my_result)))
  )
  my_result_exclude_1 <- exclude(my_result, FUN = exclude_nca_tmax_early(tmax_early = 1))
  expect_equal(
    as.data.frame(my_result_exclude_1)$exclude,
    rep("Tmax is <=1 (likely missed dose, insufficient PK samples, or PK sample swap)", nrow(as.data.frame(my_result_exclude_1)))
  )
  my_result_exclude_0 <- exclude(my_result, FUN = exclude_nca_tmax_0())
  expect_equal(
    as.data.frame(my_result_exclude_0)$exclude,
    rep("Tmax is <=0 (likely missed dose, insufficient PK samples, or PK sample swap)", nrow(as.data.frame(my_result_exclude_0)))
  )
  # This should never happen in real code
  expect_error(
    exclude_nca_tmax_early()(data.frame(PPTESTCD = "tmax", PPORRES = 1:2)),
    regexp = "Should not see more than one tmax (please report this as a bug)",
    fixed = TRUE
  )
})

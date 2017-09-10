context("exclude_nca")

test_that("exclude_nca", {
  my_conc <- PKNCAconc(data.frame(conc=c(1.1^(3:0), 1.1), time=0:4, subject=1), conc~time|subject)
  my_data <- PKNCAdata(my_conc, intervals=data.frame(start=0, end=Inf, aucinf.obs=TRUE, aucpext.obs=TRUE))
  my_result <- pk.nca(my_data)

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
                 rep("Span ratio < 2", 10)))
  my_result_excluded <- exclude(my_result, FUN=exclude_nca_span.ratio(1))
  expect_equal(as.data.frame(my_result_excluded)$exclude,
               c(rep(NA_character_, 4),
                 rep("Span ratio < 1", 10)))
  
  my_result_excluded <- exclude(my_result, FUN=exclude_nca_min.hl.r.squared())
  expect_equal(as.data.frame(my_result_excluded)$exclude,
               c(rep(NA_character_, 4),
                 rep("Half-life r-squared < 0.9", 10)))
  my_result_excluded <- exclude(my_result, FUN=exclude_nca_min.hl.r.squared(0.95))
  expect_equal(as.data.frame(my_result_excluded)$exclude,
               c(rep(NA_character_, 4),
                 rep("Half-life r-squared < 0.95", 10)))
  
  my_data <- PKNCAdata(my_conc, intervals=data.frame(start=0, end=Inf, cmax=TRUE))
  my_result <- pk.nca(my_data)
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
test_that("assert_intervals works with valid intervals", {
  o_conc <- PKNCAconc(as.data.frame(datasets::Theoph), conc~Time|Subject)
  o_data <- PKNCAdata(o_conc, intervals = data.frame(start = 0, end = 1, cmax = TRUE))
  
  result <- assert_intervals(intervals = data.frame(start = 0, end = 1, cmax = TRUE), o_data)
  expect_equal(result, expected = data.frame(start = 0, end = 1, cmax = TRUE))
})

test_that("assert_intervals works with valid intervals (ungrouped)", {
  o_conc <- PKNCAconc(as.data.frame(datasets::Theoph)[datasets::Theoph$Subject == 1,], conc~Time)
  o_data <- PKNCAdata(o_conc, intervals = data.frame(start = 0, end = 1, cmax = TRUE))
  
  result <- assert_intervals(intervals = data.frame(start = 0, end = 1, cmax = TRUE), o_data)
  expect_equal(result, expected = data.frame(start = 0, end = 1, cmax = TRUE))
})


test_that("assert_intervals errors with non-data frame intervals", {
  o_conc <- PKNCAconc(as.data.frame(datasets::Theoph), conc~Time|Subject)
  o_data <- PKNCAdata(o_conc, intervals = data.frame(start = 0, end = 1, cmax = TRUE))
  non_df_intervals <- list(a = 1, b = 2)
  
  expect_error(assert_intervals(non_df_intervals, o_data), 
               "The 'intervals' argument must be a data frame or a data frame-like object.")
})

test_that("assert_intervals errors with non-data frame intervals (ungrouped)", {
  o_conc <- PKNCAconc(as.data.frame(datasets::Theoph)[datasets::Theoph$Subject == 1,], conc~Time)
  o_data <- PKNCAdata(o_conc, intervals = data.frame(start = 0, end = 1, cmax = TRUE))
  non_df_intervals <- list(a = 1, b = 2)
  
  expect_error(assert_intervals(non_df_intervals, o_data), 
               "The 'intervals' argument must be a data frame or a data frame-like object.")
})

test_that("assert_intervals errors with non-PKNCAdata data object", {
  expect_error(assert_intervals(intervals = data.frame(start = 0, end = 1, cmax = TRUE), 
               data = data.frame(a = 1, b = 2)),
               "The 'data' argument must be a PKNCAdata object.")
})

test_that("assert_intervals errors with invalid columns", {
  o_conc <- PKNCAconc(as.data.frame(datasets::Theoph), conc~Time|Subject)
  o_data <- PKNCAdata(o_conc, intervals = data.frame(start = 0, end = 1, cmax = TRUE))
  
  invalid_intervals <- data.frame(
    mean = TRUE,  # Not allowed NCA params
    median = TRUE
  )
  
  expect_error(assert_intervals(invalid_intervals, o_data), 
               "The following columns in 'intervals' are not allowed:")
})

test_that("assert_intervals errors with invalid columns", {
  o_conc <- PKNCAconc(as.data.frame(datasets::Theoph)[datasets::Theoph$Subject == 1,], conc~Time)
  o_data <- PKNCAdata(o_conc, intervals = data.frame(start = 0, end = 1, cmax = TRUE))
  
  invalid_intervals <- data.frame(
    mean = TRUE,  # Not allowed NCA params
    median = TRUE
  )
  
  expect_error(assert_intervals(invalid_intervals, o_data), 
               "The following columns in 'intervals' are not allowed:")
})

test_that("set_intervals works with valid intervals", {
  o_conc <- PKNCAconc(as.data.frame(datasets::Theoph), conc~Time|Subject)
  o_data <- PKNCAdata(o_conc, intervals = data.frame(start = 0, end = 1, cmax = TRUE))
  
  result <- set_intervals(o_data, intervals = data.frame(start = 0, end = 1, cmin = TRUE))
  
  expect_equal(result$intervals, data.frame(start = 0, end = 1, cmin = TRUE))
})

test_that("set_intervals works with valid intervals (ungrouped)", {
  o_conc <- PKNCAconc(as.data.frame(datasets::Theoph)[datasets::Theoph$Subject == 1,], conc~Time)
  o_data <- PKNCAdata(o_conc, intervals = data.frame(start = 0, end = 1, cmax = TRUE))
  
  result <- set_intervals(o_data, intervals = data.frame(start = 0, end = 1, cmin = TRUE))
  
  expect_equal(result$intervals, data.frame(start = 0, end = 1, cmin = TRUE))
})

test_that("set_intervals fails with invalid intervals", {
  o_conc <- PKNCAconc(as.data.frame(datasets::Theoph), conc~Time|Subject)
  o_data <- PKNCAdata(o_conc, intervals = data.frame(start = 0, end = 1, cmax = TRUE))
  
  expect_error(set_intervals(o_data, intervals = data.frame(start = 0, end = 1, cmedian = TRUE)), 
               "The following columns in 'intervals' are not allowed:")
})

test_that("set_intervals fails when not using PKNCAdata", {
  o_conc <- PKNCAconc(as.data.frame(datasets::Theoph), conc~Time|Subject)
  o_data <- PKNCAdata(o_conc, intervals = data.frame(start = 0, end = 1, cmax = TRUE))
  
  expect_error(set_intervals(o_conc, intervals = data.frame(start = 0, end = 1, cmin = TRUE)), 
               "The 'data' argument must be a PKNCAdata object.")
})

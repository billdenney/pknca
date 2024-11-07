test_that("assert_intervals works with valid intervals", {
  o_conc <- PKNCAconc(as.data.frame(datasets::Theoph), conc~Time|Subject)
  o_data <- PKNCAdata(o_conc, intervals = data.frame(start = 0, end = 1, cmax = TRUE))
  
  valid_intervals <- data.frame(
  )
  
  result <- assert_intervals(o_data, valid_intervals)
  expect_equal(result, valid_intervals)
})

test_that("assert_intervals errors with non-data frame intervals", {
  o_conc <- PKNCAconc(as.data.frame(datasets::Theoph), conc~Time|Subject)
  o_data <- PKNCAdata(o_conc, intervals = data.frame(start = 0, end = 1, cmax = TRUE))
  non_df_intervals <- list(a = 1, b = 2)
  
  expect_error(assert_intervals(o_data, non_df_intervals), 
               "The 'intervals' argument must be a data frame or a data frame-like object.")
})

test_that("assert_intervals errors with invalid columns", {
  o_conc <- PKNCAconc(as.data.frame(datasets::Theoph), conc~Time|Subject)
  o_data <- PKNCAdata(o_conc, intervals = data.frame(start = 0, end = 1, cmax = TRUE))
  
  invalid_intervals <- data.frame(
    mean = TRUE,  # Not allowed NCA params
    median = TRUE
  )
  
  expect_error(assert_intervals(o_data, invalid_intervals), 
               "The following columns in 'intervals' are not allowed:")
})

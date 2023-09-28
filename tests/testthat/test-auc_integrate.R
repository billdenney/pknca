test_that("choose_interval_method", {
  # Increasing, no zeros, AUCinf ####
  expect_equal(
    choose_interval_method(conc = 1:2, time = 1:2, method = "linear", auc.type = "AUCinf"),
    c("linear", "extrap_log")
  )
  expect_equal(
    choose_interval_method(conc = 1:3, time = 1:3, method = "linear", auc.type = "AUCinf"),
    c("linear", "linear", "extrap_log")
  )
  expect_equal(
    choose_interval_method(conc = 1:2, time = 1:2, method = "lin up/log down", auc.type = "AUCinf"),
    c("linear", "extrap_log")
  )
  expect_equal(
    choose_interval_method(conc = 1:3, time = 1:3, method = "lin up/log down", auc.type = "AUCinf"),
    c("linear", "linear", "extrap_log")
  )
  # Increasing, no zeros, AUCall ####
  expect_equal(
    choose_interval_method(conc = 1:2, time = 1:2, method = "linear", auc.type = "AUCall"),
    c("linear", "zero")
  )
  expect_equal(
    choose_interval_method(conc = 1:3, time = 1:3, method = "linear", auc.type = "AUCall"),
    c("linear", "linear", "zero")
  )
  # Decreasing, no zeros ####
  expect_equal(
    choose_interval_method(conc = 2:1, time = 1:2, method = "linear", auc.type = "AUCinf"),
    c("linear", "extrap_log")
  )
  expect_equal(
    choose_interval_method(conc = 3:1, time = 1:3, method = "linear", auc.type = "AUCinf"),
    c("linear", "linear", "extrap_log")
  )
  expect_equal(
    choose_interval_method(conc = 2:1, time = 1:2, method = "lin up/log down", auc.type = "AUCinf"),
    c("log", "extrap_log")
  )
  expect_equal(
    choose_interval_method(conc = 3:1, time = 1:3, method = "lin up/log down", auc.type = "AUCinf"),
    c("log", "log", "extrap_log")
  )
  # Increasing, one initial zero ####
  expect_equal(
    choose_interval_method(conc = c(0, 1:2), time = 1:3, method = "linear", auc.type = "AUCinf"),
    c("linear", "linear", "extrap_log")
  )
  expect_equal(
    choose_interval_method(conc = c(0, 1:2), time = 1:3, method = "lin up/log down", auc.type = "AUCinf"),
    c("linear", "linear", "extrap_log")
  )
  # Increasing, two initial zeros ####
  expect_equal(
    choose_interval_method(conc = c(0, 0, 1:2), time = 1:4, method = "linear", auc.type = "AUCinf"),
    c("zero", "linear", "linear", "extrap_log")
  )
  expect_equal(
    choose_interval_method(conc = c(0, 0, 1:2), time = 1:4, method = "lin up/log down", auc.type = "AUCinf"),
    c("zero", "linear", "linear", "extrap_log")
  )
  # Increasing and decreasing, two initial zeros ####
  expect_equal(
    choose_interval_method(conc = c(0, 0, 1:2, 1), time = 1:5, method = "linear", auc.type = "AUCinf"),
    c("zero", "linear", "linear", "linear", "extrap_log")
  )
  expect_equal(
    choose_interval_method(conc = c(0, 0, 1:2, 1), time = 1:5, method = "lin up/log down", auc.type = "AUCinf"),
    c("zero", "linear", "linear", "log", "extrap_log")
  )
  # Increasing and decreasing, two initial zeros, one final zero, AUCinf ####
  expect_equal(
    choose_interval_method(conc = c(0, 0, 1:2, 1, 0), time = 1:6, method = "linear", auc.type = "AUCinf"),
    c("zero", "linear", "linear", "linear", "zero", "extrap_log")
  )
  expect_equal(
    choose_interval_method(conc = c(0, 0, 1:2, 1, 0), time = 1:6, method = "lin up/log down", auc.type = "AUCinf"),
    c("zero", "linear", "linear", "log", "zero", "extrap_log")
  )
  # Increasing and decreasing, two initial zeros, one final zero, AUClast ####
  expect_equal(
    choose_interval_method(conc = c(0, 0, 1:2, 1, 0), time = 1:6, method = "linear", auc.type = "AUClast"),
    c("zero", "linear", "linear", "linear", "zero", "zero")
  )
  expect_equal(
    choose_interval_method(conc = c(0, 0, 1:2, 1, 0), time = 1:6, method = "lin up/log down", auc.type = "AUClast"),
    c("zero", "linear", "linear", "log", "zero", "zero")
  )
  # Increasing and decreasing, two initial zeros, one final zero, AUCall ####
  expect_equal(
    choose_interval_method(conc = c(0, 0, 1:2, 1, 0), time = 1:6, method = "linear", auc.type = "AUCall"),
    c("zero", "linear", "linear", "linear", "linear", "zero")
  )
  expect_equal(
    choose_interval_method(conc = c(0, 0, 1:2, 1, 0), time = 1:6, method = "lin up/log down", auc.type = "AUCall"),
    c("zero", "linear", "linear", "log", "linear", "zero")
  )
  # Increasing and decreasing, two initial zeros, two final zeros, AUCinf ####
  expect_equal(
    choose_interval_method(conc = c(0, 0, 1:2, 1, 0, 0), time = 1:7, method = "linear", auc.type = "AUCinf"),
    c("zero", "linear", "linear", "linear", "zero", "zero", "extrap_log")
  )
  expect_equal(
    choose_interval_method(conc = c(0, 0, 1:2, 1, 0, 0), time = 1:7, method = "lin up/log down", auc.type = "AUCinf"),
    c("zero", "linear", "linear", "log", "zero", "zero", "extrap_log")
  )

  # Increasing and decreasing, two initial zeros, one middle zero, two final zeros, AUCinf ####
  expect_equal(
    choose_interval_method(conc = c(0, 0, 1:2, 0, 2:1, 0, 0), time = 1:9, method = "linear", auc.type = "AUCinf"),
    c("zero", "linear", "linear", "linear", "linear", "linear", "zero", "zero", "extrap_log")
  )
  expect_equal(
    choose_interval_method(conc = c(0, 0, 1:2, 0, 2:1, 0, 0), time = 1:9, method = "lin up/log down", auc.type = "AUCinf"),
    c("zero", "linear", "linear", "linear", "linear", "log", "zero", "zero", "extrap_log")
  )
  # Increasing and decreasing, two initial zeros, one middle zero, two final zeros, AUCall ####
  expect_equal(
    choose_interval_method(conc = c(0, 0, 1:2, 0, 2:1, 0, 0), time = 1:9, method = "linear", auc.type = "AUCall"),
    c("zero", "linear", "linear", "linear", "linear", "linear", "linear", "zero", "zero")
  )
  expect_equal(
    choose_interval_method(conc = c(0, 0, 1:2, 0, 2:1, 0, 0), time = 1:9, method = "lin up/log down", auc.type = "AUCall"),
    c("zero", "linear", "linear", "linear", "linear", "log", "linear", "zero", "zero")
  )
  # Increasing and decreasing, two initial zeros, one middle zero, two final zeros, AUClast ####
  expect_equal(
    choose_interval_method(conc = c(0, 0, 1:2, 0, 2:1, 0, 0), time = 1:9, method = "linear", auc.type = "AUClast"),
    c("zero", "linear", "linear", "linear", "linear", "linear", "zero", "zero", "zero")
  )
  expect_equal(
    choose_interval_method(conc = c(0, 0, 1:2, 0, 2:1, 0, 0), time = 1:9, method = "lin up/log down", auc.type = "AUClast"),
    c("zero", "linear", "linear", "linear", "linear", "log", "zero", "zero", "zero")
  )
})

test_that("choose_interval_method expected errors", {
  expect_error(choose_interval_method())
  expect_error(choose_interval_method(conc = "A"))
  expect_error(choose_interval_method(conc = 1, time = "A"))
  expect_error(choose_interval_method(conc = NA_real_, time = 1))
  expect_error(choose_interval_method(conc = 1, time = NA_real_))
  expect_error(choose_interval_method(conc = 1, time = 1:2))
  expect_error(choose_interval_method(conc = 1:2, time = 1))
  expect_error(choose_interval_method(conc = 1:2, time = 1:2))
  expect_error(choose_interval_method(conc = 1:2, time = 1:2, method = c("foo", "bar")))
  expect_error(choose_interval_method(conc = 1:2, time = 1:2, method = "foo"))
  expect_error(choose_interval_method(conc = 1:2, time = 1:2, method = "linear"))
  expect_error(choose_interval_method(conc = 1:2, time = 1:2, method = "linear", auc.type = c("foo", "bar")))
  expect_error(choose_interval_method(conc = 1:2, time = 1:2, method = "linear", auc.type = "foo"))
  expect_error(choose_interval_method(conc = 1:2, time = 1:2, method = "linear", auc.type = "AUCinf", tlast = "A"))
  expect_error(choose_interval_method(conc = 1:2, time = 1:2, method = "linear", auc.type = "AUCinf", tlast = 1:2))
})

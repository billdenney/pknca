# Create sample data for testing
d_conc <- data.frame(
  conc = c(1, 0.6, 0.2, 0.1, 0.9, 0.4, 1.2, 0.8, 0.3, 0.2, 1.1, 0.5),
  time = rep(0:5, 2),
  analyte = rep(c("Analyte1", "Analyte2"), each = 6),
  include_hl = c(FALSE, NA, TRUE, TRUE, TRUE, TRUE, FALSE, NA, TRUE, TRUE, TRUE, TRUE),
  ID = rep(1:2, each = 6)
)

d_dose <- data.frame(
  dose = c(100, 200),
  time = c(0, 0),
  ID = c(1, 2)
)

intervals <- data.frame(
  start = c(0, 0, 0),
  end = c(24, 48, Inf),
  half.life = c(TRUE, TRUE, TRUE),
  cmax = c(TRUE, TRUE, TRUE),
  impute = c("start_conc0,start_predose", "start_predose", "start_conc0"),
  analyte = c("Analyte1", "Analyte2", "Analyte1"),
  ID = c(1, 2, 1)
)

o_conc <- PKNCAconc(d_conc, conc ~ time | ID / analyte, include_half.life = "include_hl")
o_dose <- PKNCAdose(d_dose, dose ~ time | ID)
o_data <- PKNCAdata(o_conc, o_dose, intervals = intervals)


### Test interval_add_impute

test_that("interval_add_impute throws an error if either data or target_impute is missing", {
  expect_error(interval_add_impute(), "Both 'data' and 'target_impute' must be provided.")
  expect_error(interval_add_impute(o_data), "Both 'data' and 'target_impute' must be provided.")
  expect_error(interval_add_impute(target_impute = "start_conc0"), "Both 'data' and 'target_impute' must be provided.")
})

test_that("interval_add_impute throws an error for non-character target_impute", {
  expect_error(interval_add_impute(o_data, target_impute = 123), "'target_impute' must be a character string.")
})

test_that("interval_add_impute throws an error when input data is a non PKNCAdata object or has no intervals", {
  expect_error(interval_add_impute(o_data$conc, target_impute = "start_conc0"), "The 'data' object must be a PKNCAdata object or a data frame")
  expect_no_error(interval_add_impute(data = o_data, target_impute = "start_conc0"))
})

test_that("interval_add_impute throws an error for unknown target_params", {
  expect_error(interval_add_impute(o_data, target_impute = "start_conc0", target_params = "unknown_param"), 
               "The following target_params are not interval columns and/or known PKNCA parameters: unknown_param")
})

test_that("interval_add_impute handles impute column with different names", {
  o_data_changed_impute_name <- o_data
  o_data_changed_impute_name$impute <- "impute_col"
  o_data_changed_impute_name$intervals <- o_data_changed_impute_name$intervals %>% rename(impute_col = impute)
  result <- interval_add_impute(o_data_changed_impute_name, target_impute = "new_impute")
  expect_equal(result$intervals %>% select(analyte, half.life, cmax, impute_col),
               data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1"),
                          half.life = c(TRUE, TRUE, TRUE),
                          cmax = c(TRUE, TRUE, TRUE),
                          impute_col = c("start_conc0,start_predose,new_impute", "start_predose,new_impute", "start_conc0,new_impute")))
})

test_that("interval_add_impute handles impute column with NA values correctly", {
  o_data_with_na_impute <- o_data
  o_data_with_na_impute$intervals <- o_data_with_na_impute$intervals %>% mutate(impute = NA_character_)
  result <- interval_add_impute(o_data_with_na_impute, target_impute = "new_impute")
  expect_equal(result$intervals %>% select(analyte, half.life, cmax, impute),
               data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1"),
                          half.life = c(TRUE, TRUE, TRUE),
                          cmax = c(TRUE, TRUE, TRUE),
                          impute = c("new_impute", "new_impute", "new_impute")))
})

test_that("interval_add_impute with no optional parameters uses all relevant cases, with new intervals below", {
  result <- interval_add_impute(o_data, target_impute = "new_impute")
  expect_equal(result$intervals %>% select(analyte, half.life, cmax, impute),
               data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1"),
                          half.life = c(TRUE, TRUE, TRUE),
                          cmax = c(TRUE, TRUE, TRUE),
                          impute = c("start_conc0,start_predose,new_impute", "start_predose,new_impute", "start_conc0,new_impute")))
})

test_that("interval_add_impute handles specified target_params correctly", {
  result <- interval_add_impute(o_data, target_impute = "new_impute", target_params = "half.life")
  expect_equal(result$intervals %>% filter(half.life) %>% select(analyte, half.life, impute),
               data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1"),
                          half.life = c(TRUE, TRUE, TRUE),
                          impute = c("start_conc0,start_predose,new_impute", "start_predose,new_impute", "start_conc0,new_impute")))
  expect_equal(result$intervals %>% filter(cmax) %>% select(analyte, cmax, impute),
               o_data$intervals %>% filter(cmax) %>% select(analyte, cmax, impute))
})

test_that("interval_add_impute handles target_groups correctly", {
  result <- interval_add_impute(o_data, target_impute = "new_impute", target_groups = data.frame(analyte = "Analyte1"))
  expect_equal(result$intervals %>% filter(analyte == "Analyte1") %>% select(analyte, half.life, cmax, impute),
               data.frame(analyte = c("Analyte1", "Analyte1"),
                          half.life = c(TRUE, TRUE),
                          cmax = c(TRUE, TRUE),
                          impute = c("start_conc0,start_predose,new_impute", "start_conc0,new_impute")))
  expect_equal(result$intervals %>% filter(analyte == "Analyte2") %>% select(analyte, half.life, cmax, impute),
               o_data$intervals %>% filter(analyte == "Analyte2") %>% select(analyte, half.life, cmax, impute))
})

test_that("interval_add_impute handles multiple target_params correctly", {
  result <- interval_add_impute(o_data, target_impute = "new_impute", target_params = c("half.life", "cmax"))
  expect_equal(result$intervals %>% select(analyte, half.life, cmax, impute),
               data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1"),
                          half.life = c(TRUE, TRUE, TRUE),
                          cmax = c(TRUE, TRUE, TRUE),
                          impute = c("start_conc0,start_predose,new_impute", "start_predose,new_impute", "start_conc0,new_impute")))
})

test_that("interval_add_impute handles allow_duplication correctly", {
  
  # When allow_duplication is FALSE, intervals with already the same impute method do not add it
  result1 <- interval_add_impute(o_data, target_impute = "start_conc0", allow_duplication = FALSE)
  expect_equal(result1$intervals %>% select(analyte, half.life, cmax, impute),
               data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1"),
                          half.life = c(TRUE, TRUE, TRUE),
                          cmax = c(TRUE, TRUE, TRUE),
                          impute = c("start_conc0,start_predose", "start_predose,start_conc0", "start_conc0")))
  
  # When allow_duplication is TRUE, intervals with already the same impute method add it
  result2 <- interval_add_impute(o_data, target_impute = "start_conc0", allow_duplication = TRUE)
  expect_equal(result2$intervals %>% select(analyte, half.life, cmax, impute),
               data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1"),
                          half.life = c(TRUE, TRUE, TRUE),
                          cmax = c(TRUE, TRUE, TRUE),
                          impute = c("start_conc0,start_predose,start_conc0", "start_predose,start_conc0", "start_conc0,start_conc0")))
  
})

test_that("interval_add_impute handles correctly argument new_rows_after_original", {
  
  # When true the new rows are added after the original rows
  result1 <- interval_add_impute(o_data, target_impute = "new_impute", target_param = "cmax", new_rows_after_original = TRUE)
  expect_equal(result1$intervals %>% select(analyte, half.life, cmax, impute),
               data.frame(analyte = c("Analyte1", "Analyte1", "Analyte2", "Analyte2", "Analyte1", "Analyte1"),
                          half.life = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
                          cmax = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE),
                          impute = c("start_conc0,start_predose", 
                                     "start_conc0,start_predose,new_impute", 
                                     "start_predose", 
                                     "start_predose,new_impute", 
                                     "start_conc0", 
                                     "start_conc0,new_impute"))
  )
  
  
  # When false the new rows are added at the end of the data frame
  result2 <- interval_add_impute(o_data, target_impute = "new_impute", target_param = "cmax", new_rows_after_original = FALSE)
  expect_equal(result2$intervals %>% select(analyte, half.life, cmax, impute),
               data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1", "Analyte1", "Analyte2", "Analyte1"),
                          half.life = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
                          cmax = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
                          impute = c("start_conc0,start_predose", 
                                     "start_predose", 
                                     "start_conc0", 
                                     "start_conc0,start_predose,new_impute", 
                                     "start_predose,new_impute", 
                                     "start_conc0,new_impute"))
  )
})


### Test interval_remove_impute


test_that("interval_remove_impute throws an error if either data or target_impute is missing", {
  expect_error(interval_remove_impute(), "Both 'data' and 'target_impute' must be provided.")
  expect_error(interval_remove_impute(o_data), "Both 'data' and 'target_impute' must be provided.")
  expect_error(interval_remove_impute(target_impute = "start_conc0"), "Both 'data' and 'target_impute' must be provided.")
})

test_that("interval_remove_impute throws an error for non-character target_impute", {
  expect_error(interval_remove_impute(o_data, target_impute = 123), "'target_impute' must be a character string.")
})

test_that("interval_remove_impute throws an error when input data is a non PKNCAdata object or has no intervals", {
  expect_error(interval_remove_impute(o_data$conc, target_impute = "start_conc0"), "The 'data' object must be a PKNCAdata object or a data frame")
  expect_no_error(interval_remove_impute(data = o_data, target_impute = "start_conc0"))
})

test_that("interval_remove_impute throws an error for unknown target_params", {
  expect_error(interval_remove_impute(o_data, target_impute = "start_conc0", target_params = "unknown_param"), 
               "The following target_params are not interval columns and/or known PKNCA parameters: unknown_param")
})

test_that("interval_remove_impute handles impute column with different names", {
  o_data_changed_impute_name <- o_data
  o_data_changed_impute_name$impute <- "impute_col"
  o_data_changed_impute_name$intervals <- o_data_changed_impute_name$intervals %>% rename(impute_col = impute)
  result <- interval_remove_impute(o_data_changed_impute_name, target_impute = "start_conc0")
  expect_equal(result$intervals %>% select(analyte, half.life, cmax, impute_col),
               data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1"),
                          half.life = c(TRUE, TRUE, TRUE),
                          cmax = c(TRUE, TRUE, TRUE),
                          impute_col = c("start_predose", "start_predose", NA)))
})

test_that("interval_remove_impute handles impute column with NA values correctly", {
  o_data_with_na_impute <- o_data
  o_data_with_na_impute$intervals <- o_data_with_na_impute$intervals %>% mutate(impute = NA_character_)
  result <- interval_remove_impute(o_data_with_na_impute, target_impute = "start_conc0")
  expect_equal(result$intervals %>% select(analyte, half.life, cmax, impute),
               data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1"),
                          half.life = c(TRUE, TRUE, TRUE),
                          cmax = c(TRUE, TRUE, TRUE),
                          impute = c(NA_character_, NA_character_, NA_character_)))
})

test_that("interval_remove_impute handles missing impute column by not modifying the dataset and warns the user", {
  o_data_no_impute <- o_data
  o_data_no_impute$intervals <- o_data_no_impute$intervals %>% select(-impute)
  result <- interval_remove_impute(o_data_no_impute, target_impute = "start_conc0")
  expect_equal(result, o_data_no_impute)
  expect_warning(interval_remove_impute(o_data_no_impute, target_impute = "start_conc0"),
                 "No default impute column identified. No impute methods to remove. If there is an impute column, please specify it in argument 'impute_column'")
})

# Test intervals for expected outputs with different inputs

test_that("interval_remove_impute with no optional parameters uses all relevant cases, with new intervals below", {
  result <- interval_remove_impute(o_data, target_impute = "start_conc0")
  expect_equal(result$intervals %>% select(analyte, half.life, cmax, impute),
               data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1"),
                          half.life = c(TRUE, TRUE, TRUE),
                          cmax = c(TRUE, TRUE, TRUE),
                          impute = c("start_predose", "start_predose", NA)))
})

test_that("interval_remove_impute handles specified target_params correctly", {
  result <- interval_remove_impute(o_data, target_impute = "start_conc0", target_params = "half.life")
  # half.life has no start_conc0 imputations
  expect_equal(result$intervals %>% filter(half.life) %>% select(analyte, half.life, impute),
               data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1"),
                          half.life = c(TRUE, TRUE, TRUE),
                          impute = c("start_predose", "start_predose", NA)))
  # cmax has the same exact imputations as before
  expect_equal(result$intervals %>% filter(cmax) %>% select(analyte, cmax, impute),
               o_data$intervals %>% filter(cmax) %>% select(analyte, cmax, impute))
})

test_that("interval_remove_impute handles target_groups correctly", {
  result <- interval_remove_impute(o_data, target_impute = "start_conc0", target_groups = data.frame(analyte = "Analyte1"))
  # Analyte1 has no start_conc0 imputations
  expect_equal(result$intervals %>% filter(analyte == "Analyte1") %>% select(analyte, half.life, cmax, impute),
               data.frame(analyte = c("Analyte1", "Analyte1"),
                          half.life = c(TRUE, TRUE),
                          cmax = c(TRUE, TRUE),
                          impute = c("start_predose", NA_character_)))  
  
  # Analyte2 has the same exact imputations as before
  expect_equal(result$intervals %>% filter(analyte == "Analyte2") %>% select(analyte, half.life, cmax, impute),
               o_data$intervals %>% filter(analyte == "Analyte2") %>% select(analyte, half.life, cmax, impute))
})

test_that("interval_remove_impute handles multiple target_params correctly", {
  result <- interval_remove_impute(o_data, target_impute = "start_conc0", target_params = c("half.life", "cmax"))
  expect_equal(result$intervals %>% select(analyte, half.life, cmax, impute),
               data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1"),
                          half.life = c(TRUE, TRUE, TRUE),
                          cmax = c(TRUE, TRUE, TRUE),
                          impute = c("start_predose", "start_predose", NA)))
})

test_that("interval_remove_impute handles with specificity impute character method with multiple imputes", {
  o_data_multiple_imputes <- o_data
  o_data_multiple_imputes$intervals <- o_data_multiple_imputes$intervals %>% mutate(impute = "start_conc0,start_predose")
  result <- interval_remove_impute(o_data_multiple_imputes, target_impute = "start_conc0")
  expect_equal(result$intervals %>% select(analyte, half.life, cmax, impute),
               data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1"),
                          half.life = c(TRUE, TRUE, TRUE),
                          cmax = c(TRUE, TRUE, TRUE),
                          impute = c("start_predose", "start_predose", "start_predose")))
})


test_that("interval_remove_impute handles correctly argument new_rows_after_original", {
  
  # When true the new rows are added after the original rows
  result1 <- interval_remove_impute(o_data, target_impute = "start_conc0", target_params = "cmax", new_rows_after_original = TRUE)
  expect_equal(result1$intervals %>% select(analyte, half.life, cmax, impute),
               data.frame(analyte = c("Analyte1", "Analyte1", "Analyte2", "Analyte1", "Analyte1"),
                          half.life = c(FALSE, TRUE, TRUE, FALSE, TRUE),
                          cmax = c(TRUE, FALSE, TRUE, TRUE, FALSE),
                          impute = c("start_predose", "start_conc0,start_predose", "start_predose", NA, "start_conc0")))
  
  # When false the new rows are added at the end of the data frame
  result2 <- interval_remove_impute(o_data, target_impute = "start_conc0", target_params = "cmax", new_rows_after_original = FALSE)
  expect_equal(result2$intervals %>% select(analyte, half.life, cmax, impute),
               data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1", "Analyte1", "Analyte1"),
                          half.life = c(FALSE, TRUE, FALSE, TRUE, TRUE),
                          cmax = c(TRUE, TRUE, TRUE, FALSE, FALSE),
                          impute = c("start_predose", "start_predose", NA, "start_conc0,start_predose", "start_conc0")))
})


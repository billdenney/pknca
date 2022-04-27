test_that("choose_interp_extrap_method chooses correctly", {
  conc_values <-
    list(
      all_zero=rep(0, 5),
      start_zero=c(0, 1, 2, 1, 0),
      middle_zero_before_tmax=c(1, 0, 2, 1, 0),
      middle_zero_after_tmax=c(1, 2, 0, 1, 0),
      tlast_not_end=c(1, 2, 1, 0, 0),
      tlast_end=c(1, 2, 3, 2, 1),
      single_zero=0,
      single_nonzero=1
    )
  conc_status <-
    list(
      all_zero=rep("zero", 4),
      start_zero=c("increase", "increase", "decrease", "zero"),
      middle_zero_before_tmax=c("zero", "increase", "decrease", "zero"),
      middle_zero_after_tmax=c("increase", "zero", "increase", "zero"),
      tlast_not_end=c("increase", "decrease", "zero", "zero"),
      tlast_end=c("increase", "increase", "decrease", "decrease"),
      single_zero=character(),
      single_nonzero=character()
    )
  conc_tlast <-
    list(
      all_zero=NA,
      start_zero=4,
      middle_zero_before_tmax=4,
      middle_zero_after_tmax=4,
      tlast_not_end=3,
      tlast_end=NA,
      single_zero=NA,
      single_nonzero=NA
    )
  
  interp_method_choices <- c("linear", "lin up/log down", "log")
  extrap_method_choices <- c("aucinf.obs", "aucinf.pred", "auclast", "aucall")
  for (current_conc in names(conc_values)) {
    for (current_interp in interp_method_choices) {
      for (current_extrap in extrap_method_choices) {
        result <- conc_status[[current_conc]]
        result[result %in% "increase"] <-
          if (current_interp %in% c("linear", "lin up/log down")) {
            "linear"
          } else if (current_interp %in% "log") {
            "log"
          }
        result[result %in% "decrease"] <-
          if (current_interp %in% c("linear")) {
            "linear"
          } else if (current_interp %in% c("lin up/log down", "log")) {
            "log"
          }
        result[result %in% "zero"] <- "linear"
        current_tlast <- conc_tlast[[current_conc]]
        if (!is.na(current_tlast)) {
          result[current_tlast:length(result)] <- "zero"
          if (current_extrap == "aucall") {
            result[current_tlast] <- "linear"
          }
        }
        if (current_conc %in% "all_zero") {
          result[] <- "zero"
        }
        result[length(result) + 1] <-
          if (current_extrap %in% "aucinf.obs" & !all(conc_values[[current_conc]] == 0)) {
            "clastobs"
          } else if (current_extrap %in% "aucinf.pred" & !all(conc_values[[current_conc]] == 0)) {
            "clastpred"
          } else {
            "zero"
          }
        expect_equal(
          choose_interp_extrap_method(
            conc=conc_values[[current_conc]],
            time=seq_along(conc_values[[current_conc]]),
            interp_method=current_interp,
            extrap_method=current_extrap
          ),
          result,
          info=paste(current_conc, current_interp, current_extrap)
        )
      }
    }
  }
})

test_that("choose_interp_extrap_method expected errors", {
  expect_error(choose_interp_extrap_method(conc="A"))
  expect_error(choose_interp_extrap_method(conc=1))
  expect_error(choose_interp_extrap_method(conc=1:2, time="A"))
  expect_error(choose_interp_extrap_method(conc=1:2, time=2:1))
  expect_error(choose_interp_extrap_method(conc=1:2, time=c(1, 1)))
  expect_error(choose_interp_extrap_method(conc=1:2, time=c(1, NA)))
  expect_error(choose_interp_extrap_method(conc=1:2, time=1:2, interp_method="A"))
  expect_error(choose_interp_extrap_method(conc=1:2, time=1:2, interp_method="linear", extrap_method="A"))
})

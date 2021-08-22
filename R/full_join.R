#' Combine PKNCAconc and PKNCAdose objects
#' 
#' The function is inspired by \code{dplyr::full_join}, but it has different
#' semantics.
#' 
#' @param conc a PKNCAconc object
#' @param dose a PKNCAdose object or \code{NA}
#' @return A tibble with columns for the groups, "data_conc" (the concentration
#'   data), and "data_dose" (the dosing data).  If \code{is.na(dose)},
#'   "data_dose" will be \code{NA}.
#' @family Combine PKNCA objects
#' @keywords Internal
#' @noRd
#' @importFrom dplyr full_join
#' @importFrom tidyr crossing
full_join_PKNCAconc_PKNCAdose <- function(conc, dose) {
  stopifnot(inherits(x=conc, what="PKNCAconc"))
  if (identical(dose, NA)) {
    n_dose <- tibble(data_dose=list(NA))
  } else {
    stopifnot(inherits(x=dose, what="PKNCAdose"))
    n_dose <- nest(dose)
  }
  n_conc <- nest(conc)
  shared_groups <- intersect(names(n_conc), names(n_dose))
  if (length(shared_groups) > 0) {
    dplyr::full_join(n_conc, n_dose, by=shared_groups)
  } else {
    tidyr::crossing(n_conc, n_dose)
  }
}

#' Convert a PKNCAdata object into a data.frame for analysis
#' 
#' The function is inspired by \code{dplyr::full_join}, but it has different
#' semantics.
#' 
#' @param x The PKNCAdata object
#' @return A tibble with columns the grouping variables, "data_conc" for
#'   concentration data, "data_dose" for dosing data, and "data_intervals" for
#'   intervals data.
#' @family Combine PKNCA objects
#' @keywords Internal
#' @noRd
#' @importFrom dplyr full_join
#' @importFrom tidyr crossing
full_join_PKNCAdata <- function(x) {
  conc_dose <- full_join_PKNCAconc_PKNCAdose(x$conc, x$dose)
  n_i <-
    nest_PKNCAintervals(
      .dat=x$intervals,
      vars=setdiff(names(conc_dose), c("data_conc", "data_dose", "data_intervals"))
    )
  shared_groups <- intersect(names(conc_dose), names(n_i))
  if (length(shared_groups) > 0) {
    ret <- dplyr::full_join(conc_dose, n_i, by=shared_groups)
  } else {
    ret <- tidyr::crossing(conc_dose, n_i)
  }
  ret
}

#' Nest a PKNCA object
#' 
#' @param .dat The PKNCA object to nest
#' @param ...,.names_sep,.key Ignored
#' @return A nested data.frame
#' @family Combine PKNCA objects
#' @keywords Internal
#' @exportS3Method tidyr::nest
#' @importFrom dplyr grouped_df
#' @importFrom tidyr nest
nest.PKNCAconc <- function(.dat, ..., .names_sep = NULL, .key = deprecated()) {
  check_reserved_column_names(.dat)
  as_nest <- tidyr::nest(.dat$data, data_conc=!group_vars(.dat))
  as_nest
}

#' @describeIn nest.PKNCAconc Nest a PKNCAdose object
#' @exportS3Method tidyr::nest
#' @family Combine PKNCA objects
#' @keywords Internal
#' @importFrom dplyr grouped_df
#' @importFrom tidyr nest
nest.PKNCAdose <- function(.dat, ..., .names_sep = NULL, .key = deprecated()) {
  check_reserved_column_names(.dat)
  as_nest <- tidyr::nest(.dat$data, data_dose=!group_vars(.dat))
  as_nest
}

#' @describeIn nest.PKNCAconc Nest a PKNCAdose object
#' @export
#' @family Combine PKNCA objects
#' @keywords Internal
#' @importFrom dplyr grouped_df
#' @importFrom tidyr nest
# Note that this is not a generic for nest because intervals are a data.frame
# and we do not want to use the generic data.frame method here.
nest_PKNCAintervals <- function(.dat, vars=character(0)) {
  check_reserved_column_names(.dat)
  .dat <- as_tibble(.dat)
  vars <- intersect(vars, names(.dat))
  if (length(vars) == 0) {
    as_nest <- tibble(data_intervals=list(.dat))
  } else {
    as_nest <- tidyr::nest(.dat, data_intervals=!vars)
  }
  as_nest
}

#' Confirm that PKNCA reserved column names are not in a data.frame
#'
#' @param x A data.frame or similar object
#' @return NULL (generate an error if a reserved name is present)
#' @keywords Internal
#' @noRd
check_reserved_column_names <- function(x) {
  reserved_names <- c("data_conc", "data_dose", "data_intervals")
  overlap <- intersect(reserved_names, names(x))
  if (length(overlap) > 0) {
    msg <-
      paste(
        ngettext(length(overlap), msg1="The column", msg2="The columns"),
        paste0("'", overlap, "'", collapse=", "),
        ngettext(length(overlap), msg1="is", msg2="are"),
        "reserved for internal use in PKNCA.  Change the",
        ngettext(length(overlap), msg1="name", msg2="names"),
        "and retry."
      )
    stop(msg)
  }
}

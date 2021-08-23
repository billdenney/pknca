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
#' @importFrom tibble tibble
#' @importFrom tidyr crossing
full_join_PKNCAconc_PKNCAdose <- function(conc, dose) {
  stopifnot(inherits(x=conc, what="PKNCAconc"))
  if (identical(dose, NA)) {
    message("No dose information provided, calculations requiring dose will return NA.")
    n_dose <- tibble::tibble(data_dose=list(NA))
  } else {
    stopifnot(inherits(x=dose, what="PKNCAdose"))
    n_dose <- prepare_PKNCAdose(dose)
  }
  n_conc <- prepare_PKNCAconc(conc)
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
    prepare_PKNCAintervals(
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

#' Prepare a PKNCA object and drop unnecessary columns
#' 
#' @param .dat The PKNCA object to prepare as a nested tibble
#' @param ...,.names_sep,.key Ignored
#' @return A nested tibble with a column named "data_conc" containing the concentration data and a column 
#' @family Combine PKNCA objects
#' @keywords Internal
#' @noRd
#' @importFrom dplyr grouped_df
#' @importFrom tidyr nest
prepare_PKNCA_general <- function(.dat, cols, exclude, group_cols, data_name, insert_if_missing=list()) {
  check_reserved_column_names(.dat)
  intermediate_group_cols <-
    if (length(group_cols) > 0) {
      paste0("group", seq_along(group_cols))
    } else {
      character(0)
    }
  data_no_exclude <-
    .dat[
      is.na(normalize_exclude(.dat[[exclude]])),,
      drop=FALSE
    ]
  data_standard <-
    standardize_column_names(
      x=data_no_exclude,
      cols=cols,
      group_cols=group_cols,
      insert_if_missing=insert_if_missing
    )
  # data_conc is used since it is reserved, and it will be replaced on the next
  # line.
  as_nest <- tidyr::nest(data_standard, data_conc=!intermediate_group_cols)
  names(as_nest)[names(as_nest) %in% "data_conc"] <- data_name
  ret <- restore_group_col_names(as_nest, group_cols=group_cols)
  ret
}

prepare_PKNCAconc <- function(.dat) {
  # Remove rows to be excluded from all calculations
  # Drop unnecessary column names
  pformula_conc <- parseFormula(.dat)
  needed_cols <-
    list(
      conc=all.vars(pformula_conc$lhs),
      time=all.vars(pformula_conc$rhs),
      volume=.dat$columns$volume,
      duration=.dat$columns$duration,
      include_half.life=.dat$columns$include_half.life,
      exclude_half.life=.dat$columns$exclude_half.life
    )
  ret <-
    prepare_PKNCA_general(
      .dat=.dat$data,
      exclude=.dat$exclude,
      cols=needed_cols,
      data_name="data_conc",
      group_cols=all.vars(pformula_conc$groups)
    )
  ret
}

#' @describeIn prepare_PKNCAconc Nest a PKNCAdose object
#' @noRd
#' @family Combine PKNCA objects
#' @keywords Internal
#' @importFrom dplyr grouped_df
#' @importFrom tidyr nest
prepare_PKNCAdose <- function(.dat) {
  pformula_dose <- parseFormula(.dat)
  dose_col <- all.vars(pformula_dose$lhs)
  time_col <- all.vars(pformula_dose$rhs)
  if (length(dose_col) == 0 || dose_col %in% ".") dose_col <- NULL
  if (time_col %in% ".") time_col <- NULL
  needed_cols <-
    list(
      dose=dose_col,
      time=time_col,
      duration=.dat$columns$duration,
      route=.dat$columns$route
    )
  ret <-
    prepare_PKNCA_general(
      .dat=.dat$data,
      exclude=.dat$exclude,
      cols=needed_cols,
      data_name="data_dose",
      group_cols=all.vars(pformula_dose$groups),
      insert_if_missing=list(dose=NA, time=NA)
    )
  ret
}

#' @describeIn prepare_PKNCAconc Nest a PKNCAdose object
#' @noRd
#' @family Combine PKNCA objects
#' @keywords Internal
#' @importFrom dplyr grouped_df
#' @importFrom tidyr nest
#' @importFrom tibble as_tibble tibble
prepare_PKNCAintervals <- function(.dat, vars=character(0)) {
  check_reserved_column_names(.dat)
  .dat <- tibble::as_tibble(.dat)
  vars <- intersect(vars, names(.dat))
  if (length(vars) == 0) {
    as_nest <- tibble::tibble(data_intervals=list(.dat))
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
  reserved_names <- c("data_conc", "data_dose", "data_intervals", "data_results")
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

#' Standardize column names and drop unnecessary columns from a data.frame or tibble
#' 
#' @param x The data.frame or tibble
#' @param cols A named list where the names are the standardized column names and the values are the original column names
#' @return A data.frame or tibble with columns cleaned of unlisted columns and with names set to the expected names.
#' @noRd
#' @keywords Internal
standardize_column_names <- function(x, cols, group_cols=NULL, insert_if_missing=list()) {
  stopifnot("cols must be a list"=is.list(cols))
  stopifnot("cols must be named"=!is.null(names(cols)))
  stopifnot("all cols must be named"=!any(names(cols) %in% ""))
  stopifnot("all original cols names must be names of x"=all(unlist(cols) %in% names(x)))
  stopifnot("group_cols must be NULL or a character vector"=is.null(group_cols) || is.character(group_cols))
  if (!is.null(group_cols)) {
    stopifnot("group_cols must not overlap with other column names"=!any(group_cols %in% unlist(cols)))
    stopifnot("group_cols must not overlap with standardized column names"=!any(group_cols %in% names(cols)))
    new_group_cols <- paste0("group", seq_along(group_cols))
  } else {
    new_group_cols <- NULL
  }
  cols_clean <- cols[!sapply(X=cols, FUN=is.null)]
  ret <-
    setNames(
      # Keep only columns of interest
      x[, c(group_cols, unlist(cols_clean)), drop=FALSE],
      nm=c(new_group_cols, names(cols_clean))
    )
  for (current_nm in names(insert_if_missing)) {
    if (!(current_nm %in% names(ret))) {
      ret[[current_nm]] <- insert_if_missing[[current_nm]]
    }
  }
  ret
}

restore_group_col_names <- function(x, group_cols=NULL) {
  if (is.null(group_cols)) {
    return(x)
  }
  new_group_cols <- paste0("group", seq_along(group_cols))
  stopifnot("missing intermediate group_cols names"=all(new_group_cols %in% names(x)))
  stopifnot(
    "Intermediate group_cols are out of order"=
      all(names(x)[names(x) %in% new_group_cols] == new_group_cols)
  )
  names(x)[names(x) %in% new_group_cols] <- group_cols
  x
}

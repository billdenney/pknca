#' Combine PKNCAconc and PKNCAdose objects
#'
#' The function is inspired by `dplyr::full_join`, but it has different
#' semantics.
#'
#' @param o_conc a PKNCAconc object
#' @param o_dose a PKNCAdose object or `NA`
#' @param extra_cols_conc Additional columns to include in the concentration data
#' @returns A tibble with columns for the groups, "data_conc" (the concentration
#'   data), and "data_dose" (the dosing data).  If `is.na(o_dose)`, "data_dose"
#'   will be `NA`.
#' @family Combine PKNCA objects
#' @keywords Internal
#' @noRd
full_join_PKNCAconc_PKNCAdose <- function(o_conc, o_dose, extra_cols_conc = character()) {
  stopifnot(inherits(x=o_conc, what="PKNCAconc"))
  if (identical(o_dose, NA)) {
    message("No dose information provided, calculations requiring dose will return NA.")
    n_dose <- tibble::tibble(data_dose=list(NA))
  } else {
    stopifnot(inherits(x=o_dose, what="PKNCAdose"))
    n_dose <- prepare_PKNCAdose(o_dose, sparse=is_sparse_pk(o_conc), subject_col=o_conc$columns$subject)
  }
  n_conc <- prepare_PKNCAconc(o_conc, extra_cols = extra_cols_conc)
  shared_groups <- intersect(names(n_conc), names(n_dose))
  if (length(shared_groups) > 0) {
    dplyr::full_join(n_conc, n_dose, by=shared_groups)
  } else {
    tidyr::crossing(n_conc, n_dose)
  }
}

#' Convert a PKNCAdata object into a data.frame for analysis
#'
#' The function is inspired by `dplyr::full_join`, but it has different
#' semantics.
#'
#' @param x The PKNCAdata object
#' @inheritParams full_join_PKNCAconc_PKNCAdose
#' @returns A tibble with columns the grouping variables, "data_conc" for
#'   concentration data, "data_dose" for dosing data, and "data_intervals" for
#'   intervals data.
#' @family Combine PKNCA objects
#' @keywords Internal
#' @noRd
full_join_PKNCAdata <- function(x, extra_conc_cols = character()) {
  conc_dose <- full_join_PKNCAconc_PKNCAdose(o_conc = x$conc, o_dose = x$dose, extra_cols_conc = extra_conc_cols)
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
#' @return A nested tibble with a column named "data_conc" containing the
#'   concentration data and a column
#' @family Combine PKNCA objects
#' @keywords Internal
#' @noRd
prepare_PKNCA_general <- function(.dat, cols, exclude, group_cols, data_name, insert_if_missing=list()) {
  check_reserved_column_names(.dat)
  intermediate_group_cols <-
    if (length(group_cols) > 0) {
      paste0("group", seq_along(group_cols))
    } else {
      character(0)
    }
  data_no_exclude <-
    as.data.frame(.dat[
      is.na(normalize_exclude(.dat[[exclude]])),,
      drop=FALSE
    ])
  data_standard <-
    standardize_column_names(
      x=data_no_exclude,
      cols=cols,
      group_cols=group_cols,
      insert_if_missing=insert_if_missing
    )
  # data_conc is used since it is reserved, and it will be replaced on the next
  # line.
  as_nest <- tidyr::nest(data_standard, data_conc=!dplyr::all_of(intermediate_group_cols))
  names(as_nest)[names(as_nest) %in% "data_conc"] <- data_name
  ret <- restore_group_col_names(as_nest, group_cols=group_cols)
  ret
}

prepare_PKNCAconc_dense <- function(.dat, needed_cols, group_cols_selected) {
  ret <-
    prepare_PKNCA_general(
      .dat=.dat$data,
      exclude=.dat$columns$exclude,
      cols=needed_cols,
      data_name="data_conc",
      group_cols=group_cols_selected
    )
  ret
}

prepare_PKNCAconc_sparse <- function(.dat, needed_cols, group_cols_selected) {
  needed_cols$subject <- .dat$columns$subject
  ret <-
    prepare_PKNCA_general(
      .dat=.dat$data_sparse,
      exclude=.dat$columns$exclude,
      cols=needed_cols,
      data_name="data_sparse_conc",
      group_cols=setdiff(group_cols_selected, .dat$columns$subject)
    )
  # Generate the mean profile for non-sparse parameters
  ret$data_conc <-
    ret$data_sparse_conc %>%
    lapply(FUN=as_sparse_pk) %>%
    lapply(FUN=sparse_mean) %>%
    lapply(FUN=sparse_to_dense_pk)
  ret
}

prepare_PKNCAconc <- function(.dat, extra_cols = character()) {
  # Remove rows to be excluded from all calculations
  # Drop unnecessary column names
  needed_cols <-
    list(
      conc=.dat$columns$concentration,
      time=.dat$columns$time,
      volume=.dat$columns$volume,
      duration=.dat$columns$duration,
      include_half.life=.dat$columns$include_half.life,
      exclude_half.life=.dat$columns$exclude_half.life
    )
  needed_cols <- append(needed_cols, stats::setNames(nm = extra_cols))
  data_name <- getDataName(.dat)
  group_cols_selected <- unlist(.dat$columns$groups)
  if (is_sparse_pk(.dat)) {
    ret <-
      prepare_PKNCAconc_sparse(
        .dat=.dat,
        needed_cols=needed_cols,
        group_cols_selected=group_cols_selected
      )
  } else if (data_name == "data") {
    ret <-
      prepare_PKNCAconc_dense(
        .dat=.dat,
        needed_cols=needed_cols,
        group_cols_selected=group_cols_selected
      )
  } else {
    stop("Please report this as a bug: Invalid data_name") # nocov
  }
  ret
}

#' @describeIn prepare_PKNCAconc Nest a PKNCAdose object
#'
#' @param sparse Is the data for sparse PK?
#' @param subject_col The column name indicating the subject identifier (to be
#'   dropped from groups with sparse PK)
#' @family Combine PKNCA objects
#' @keywords Internal
#' @noRd
prepare_PKNCAdose <- function(.dat, sparse, subject_col) {
  ret <- prepare_PKNCAdose_general(.dat)
  if (sparse && (length(subject_col) == 1) && (subject_col %in% names(ret))) {
    # Verify that all subjects in a group had the same data_dose and then drop
    # the subjects
    # ret_grp will have one column named "sparse_group_check" with one row per ID and all of the dosing information within a group and
    ret_grp <-
      tidyr::nest(
        ret,
        sparse_group_check=!setdiff(names(ret), c("data_dose", subject_col))
      )
    ret_grp$data_dose <- rep(list(NULL), nrow(ret_grp))
    for (current_row in seq_len(nrow(ret_grp))) {
      # Dosing information for all subjects are identical to the first subject
      all_match <-
        all(vapply(
          X=ret_grp$sparse_group_check[[current_row]]$data_dose,
          FUN=identical,
          y=ret_grp$sparse_group_check[[current_row]]$data_dose[[1]],
          FUN.VALUE = TRUE
        ))
      if (all_match) {
        # Drop the subject identifier from the dosing
        ret_grp$data_dose[[current_row]] <- ret_grp$sparse_group_check[[current_row]]$data_dose[[1]]
      } else {
        names_to_print <- setdiff(names(ret_grp), c("sparse_group_check", "data_dose"))
        msg_error_row <- paste(names_to_print, unlist(ret_grp[current_row, names_to_print]), sep="=", collapse="; ")
        msg_error <-
          if (length(names_to_print) > 0) {
            paste(
              "Not all subjects have the same dosing information for this group: ",
              msg_error_row
            )
          } else {
            "Not all subjects have the same dosing information."
          }
        stop(
          "With sparse PK, all subjects in a group must have the same dosing information.\n",
          msg_error
        )
      }
    }
    ret <- ret_grp[, setdiff(names(ret_grp), "sparse_group_check")]
  }
  ret
}

#' @describeIn prepare_PKNCAconc Nest a PKNCAdose object
#' @family Combine PKNCA objects
#' @keywords Internal
#' @noRd
prepare_PKNCAdose_general <- function(.dat) {
  dose_col <- .dat$columns$dose
  time_col <- .dat$columns$time
  if (length(dose_col) == 0) dose_col <- NULL
  if (length(time_col) == 0) time_col <- NULL
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
      exclude=.dat$columns$exclude,
      cols=needed_cols,
      data_name="data_dose",
      group_cols=unlist(.dat$columns$groups),
      insert_if_missing=list(dose=NA, time=NA)
    )
  ret
}

#' @describeIn prepare_PKNCAconc Nest a PKNCAdose object
#' @noRd
#' @family Combine PKNCA objects
#' @keywords Internal
prepare_PKNCAintervals <- function(.dat, vars=character(0)) {
  check_reserved_column_names(.dat)
  .dat <- tibble::as_tibble(.dat)
  vars <- intersect(vars, names(.dat))
  if (length(vars) == 0) {
    as_nest <- tibble::tibble(data_intervals=list(.dat))
  } else {
    as_nest <- tidyr::nest(.dat, data_intervals=!dplyr::all_of(vars))
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

#' Standardize column names and drop unnecessary columns from a data.frame or
#' tibble
#'
#' @param x The data.frame or tibble
#' @param cols A named list where the names are the standardized column names
#'   and the values are the original column names
#' @returns A data.frame or tibble with columns cleaned of unlisted columns and
#'   with names set to the expected names.
#' @noRd
#' @keywords Internal
standardize_column_names <- function(x, cols, group_cols=NULL, insert_if_missing=list()) {
  stopifnot("cols must be a list"=is.list(cols))
  stopifnot("cols must be named"=!is.null(names(cols)))
  stopifnot("all cols must be named"=!any(names(cols) %in% ""))
  stopifnot("all original cols names must be names of x"=all(unlist(cols) %in% names(x)))
  stopifnot("group_cols must be NULL or a character vector"=is.null(group_cols) || is.character(group_cols))
  if (!is.null(group_cols) && (length(group_cols) > 0)) {
    # Give a clear error message if group columns overlap
    mask_overlap_colvalues <- group_cols %in% unlist(cols)
    mask_overlap_colnames <- group_cols %in% names(cols)
    if (any(mask_overlap_colvalues)) {
      stop(
        "group_cols must not overlap with other column names.  Change the name of the following groups: ",
        paste(group_cols[mask_overlap_colvalues], collapse=", ")
      )
    }
    if (any(mask_overlap_colnames)) {
      stop(
        "group_cols must not overlap with standardized column names.  Change the name of the following groups: ",
        paste(group_cols[mask_overlap_colnames], collapse=", ")
      )
    }
    new_group_cols <- paste0("group", seq_along(group_cols))
  } else {
    new_group_cols <- NULL
  }
  cols_clean <- cols[!vapply(X = cols, FUN = is.null, FUN.VALUE = TRUE)]
  ret <-
    stats::setNames(
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
  if (is.null(group_cols) || (length(group_cols) == 0)) {
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

#' Create a PKNCAdata object.
#'
#' `PKNCAdata()` combines `PKNCAconc` and `PKNCAdose` objects and adds in the
#' intervals for PK calculations.
#'
#' @inheritParams PKNCA.choose.option
#' @param data.conc Concentration data as a `PKNCAconc` object or a data frame
#' @param data.dose Dosing data as a `PKNCAdose` object (see details)
#' @param impute Methods for imputation.  `NA` for to search for the column
#'   named "impute" in the intervals or no imputation if that column does not
#'   exist, a comma-or space-separated list of names, or the name of a column in
#'   the `intervals` data.frame.  See
#'   `vignette("v08-data-imputation", package="PKNCA")` for more details.
#' @param formula.conc Formula for making a `PKNCAconc` object with `data.conc`.
#'   This must be given if `data.conc` is a data.frame, and it must not be given
#'   if `data.conc` is a `PKNCAconc` object.
#' @param formula.dose Formula for making a `PKNCAdose` object with `data.dose`.
#'   This must be given if `data.dose` is a data.frame, and it must not be given
#'   if `data.dose` is a `PKNCAdose` object.
#' @param intervals A data frame with the AUC interval specifications as defined
#'   in [check.interval.specification()].  If missing, this will be
#'   automatically chosen by [choose.auc.intervals()]. (see details)
#' @param units A data.frame of unit assignments and conversions as created by
#'   [pknca_units_table()]
#' @param ... arguments passed to `PKNCAdata.default`
#' @returns A PKNCAdata object with concentration, dose, interval, and
#'   calculation options stored (note that PKNCAdata objects can also have
#'   results after a NCA calculations are done to the data).
#' @details If `data.dose` is not given or is `NA`, then the `intervals` must be
#'   given.  At least one of `data.dose` and `intervals` must be given.
#' @family PKNCA objects
#' @seealso [choose.auc.intervals()], [pk.nca()], [pknca_units_table()]
#' @export
PKNCAdata <- function(data.conc, data.dose, ...) {
  UseMethod("PKNCAdata", data.conc)
}

# Ensure that arguments are reversible
#' @rdname PKNCAdata
#' @export
PKNCAdata.PKNCAconc <- function(data.conc, data.dose, ...) {
  PKNCAdata.default(data.conc=data.conc, data.dose=data.dose, ...)
}

#' @rdname PKNCAdata
#' @export
PKNCAdata.PKNCAdose <- function(data.conc, data.dose, ...) {
  # Swap the arguments
  PKNCAdata.default(data.dose=data.conc, data.conc=data.dose, ...)
}

#' @rdname PKNCAdata
#' @export
PKNCAdata.default <- function(data.conc, data.dose, ...,
                              formula.conc, formula.dose,
                              impute = NA_character_,
                              intervals, units = NULL, options=list()) {
  if (length(list(...))) {
    stop("Unknown argument provided to PKNCAdata.  All arguments other than `data.conc` and `data.dose` must be named.")
  }
  ret <- list()
  # Generate the conc element
  if (inherits(data.conc, "PKNCAconc")) {
    if (!missing(formula.conc)) {
      rlang::warn(
        message = "data.conc was given as a PKNCAconc object.  Ignoring formula.conc",
        class = "pknca_dataconc_formulaconc"
      )
    }
    ret$conc <- data.conc
  } else {
    ret$conc <- PKNCAconc(data.conc, formula=formula.conc)
  }
  # Generate the dose element
  if (missing(data.dose)) {
    ret$dose <- NA
  } else if (identical(data.dose, NA)) {
    ret$dose <- NA
  } else if (inherits(data.dose, "PKNCAdose")) {
    if (!missing(formula.dose))
      rlang::warn(
        message = "data.dose was given as a PKNCAdose object.  Ignoring formula.dose",
        class = "pknca_dataconc_formuladose"
      )
    ret$dose <- data.dose
  } else {
    ret$dose <- PKNCAdose(data.dose, formula.dose)
  }
  # Check the options
  if (!is.list(options)) {
    stop("options must be a list.")
  }
  if (length(options) > 0) {
    if (is.null(names(options)))
      stop("options must have names.")
    for (n in names(options)) {
      tmp.opt <- list(options[[n]], TRUE)
      names(tmp.opt) <- c(n, "check")
      do.call(PKNCA.options, tmp.opt)
    }
  }
  ret$options <- options

  # Assign the class and give it all back to the user.
  class(ret) <- c("PKNCAdata", class(ret))

  # Check the intervals
  if (missing(intervals) & identical(ret$dose, NA)) {
    stop("If data.dose is not given, intervals must be given")
  } else if (missing(intervals)) {
    # Generate the intervals for each grouping of concentration and
    # dosing.
    if (length(ret$dose$columns$time) == 0) {
      stop("Dose times were not given, so intervals must be manually specified.")
    }
    n_conc_dose <-
      full_join_PKNCAconc_PKNCAdose(
        o_conc = ret$conc,
        o_dose = ret$dose
      )
    n_conc_dose$data_intervals <- rep(list(NULL), nrow(n_conc_dose))
    for (idx in seq_len(nrow(n_conc_dose))) {
      current_conc <- n_conc_dose$data_conc[[idx]]
      current_dose <- n_conc_dose$data_dose[[idx]]
      current_group <-
        n_conc_dose[
          idx,
          setdiff(names(n_conc_dose), c("data_conc", "data_dose")),
          drop=FALSE
        ]
      warning_prefix <-
        if (ncol(current_group) > 0) {
          paste0(
            paste(names(current_group), unlist(lapply(current_group, as.character)), sep="=", collapse="; "),
            ": "
          )
        } else {
          ""
        }
      if (!is.null(current_conc)) {
        generated_intervals <-
          choose.auc.intervals(
            current_conc$time,
            current_dose$time,
            options=options
          )
        if (nrow(generated_intervals) > 0) {
          n_conc_dose$data_intervals[[idx]] <- generated_intervals
        } else {
          warning(warning_prefix, "No intervals generated likely due to limited concentration data")
        }
      } else {
        rlang::warn(
          message = paste(warning_prefix, "No intervals generated due to no concentration data"),
          class = "pknca_no_intervals_generated"
        )
      }
    }
    intervals <-
      tidyr::unnest(
        n_conc_dose[, setdiff(names(n_conc_dose), c("data_conc", "data_dose")), drop=FALSE],
        cols="data_intervals"
      )
  }
  ret <- set_intervals(data = ret, intervals = intervals)
  ret$intervals <- check.interval.specification(intervals)
  # Verify that either everything or nothing is using units
  units_interval_start <- inherits(ret$intervals$start, "units")
  units_interval_end <- inherits(ret$intervals$end, "units")

  # Insert the unit conversion table
  if (is.null(units)) {
    # What unit types are recognized?
    possible_units <-
      setdiff(
        grep(x = names(formals(pknca_units_table)), pattern = "_", invert = TRUE, value = TRUE),
        "conversions"
      )
    possible_units_pref <- paste0(possible_units, "_pref")
    # Accumulate available units
    conc_units_values <- ret$conc$units
    conc_units_cols <- ret$conc$columns[names(ret$conc$columns) %in% possible_units]

    unit_args <- conc_units_values
    for (nm in names(conc_units_cols)) {
      unit_args[[nm]] <- unique(stats::na.omit(ret$conc$data[[conc_units_cols[[nm]]]]))
    }

    if (!identical(ret$dose, NA)) {
      unit_args <- append(unit_args, ret$dose$units)
      dose_units_cols <- ret$dose$columns[names(ret$dose$columns) %in% possible_units]
      for (nm in names(dose_units_cols)) {
        unit_args[[nm]] <- unique(stats::na.omit(ret$dose$data[[dose_units_cols[[nm]]]]))
      }
    }
    # If there are any units to set, set them here
    if (length(unit_args) > 0) {
      unit_args <- lapply(X = unit_args, FUN = drop_attributes)
      ret$units <- do.call(pknca_units_table, args = unit_args)
    }
  } else {
    stopifnot("`units` must be a data.frame"=is.data.frame(units))
    stopifnot(
      "`units` data.frame must have at least names 'PPTESTCD' and 'PPORRESU'"=
        all(c("PPTESTCD", "PPORRESU") %in% names(units))
    )
    stopifnot("`units` must have at least one row"=nrow(units) > 0)
    ret$units <- units
  }

  # Insert the imputation methods, if applicable
  if (!identical(NA, impute)) {
    checkmate::assert_character(impute, len = 1)
    ret$impute <- impute
  }

  ret
}

drop_attributes <- function(x) {
  attributes(x) <- NULL
  x
}

#' @rdname is_sparse_pk
#' @export
is_sparse_pk.PKNCAdata <- function(object) {
  is_sparse_pk(object$conc)
}

#' Print a PKNCAdata object
#' @param x The object to print
#' @param ... Arguments passed on to [print.PKNCAconc()] and [print.PKNCAdose()]
#' @export
print.PKNCAdata <- function(x, ...) {
  print.PKNCAconc(x$conc, ...)
  if (identical(NA, x$dose)) {
    cat("No dosing information.\n")
  } else {
    print.PKNCAdose(x$dose, ...)
  }
  cat(sprintf("\nWith %d rows of interval specifications.\n",
              nrow(x$intervals)))
  if (!is.null(x$units)) {
    cat("With units\n")
  }
  if (!is.null(x$impute)) {
    cat(sprintf("With imputation: %s\n", x$impute))
  }
  if (length(x$options) == 0) {
    cat("No options are set differently than default.\n")
  } else {
    cat("Options changed from default are:\n")
    print(x$options)
  }
}

#' Summarize a PKNCAdata object showing important details about the
#' concentration, dosing, and interval information.
#' @param object The PKNCAdata object to summarize.
#' @param ... arguments passed on to [print.PKNCAdata()]
#' @export
summary.PKNCAdata <- function(object, ...) {
  print.PKNCAdata(object, summarize=TRUE, ...)
}

#' Get the groups (right hand side after the `|` from a PKNCA
#' object).
#'
#' @rdname getGroups.PKNCAconc
#' @param object The object to extract the data from
#' @param ... Arguments passed to other getGroups functions
#' @returns A data frame with the (selected) group columns.
#' @export
getGroups.PKNCAdata <- function(object, ...) {
  getGroups(as_PKNCAconc(object), ...)
}

#' @describeIn group_vars.PKNCAconc Get group_vars for a PKNCAdata object
#'   from the PKNCAconc object within
#' @exportS3Method dplyr::group_vars
group_vars.PKNCAdata <- function(x) {
  group_vars.PKNCAconc(as_PKNCAconc(x))
}

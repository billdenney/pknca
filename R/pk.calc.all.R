#' Compute NCA parameters for each interval for each subject.
#'
#' The `pk.nca` function computes the NCA parameters from a `PKNCAdata` object.
#' All options for the calculation and input data are set in prior functions
#' (`PKNCAconc`, `PKNCAdose`, and `PKNCAdata`).  Options for calculations are
#' set either in `PKNCAdata` or with the current default options in
#' `PKNCA.options`.
#'
#' When performing calculations, all time results are relative to the start of
#' the interval.  For example, if an interval starts at 168 hours, ends at 192
#' hours, and and the maximum concentration is at 169 hours, `tmax=169-168=1`.
#'
#' @param data A PKNCAdata object
#' @param verbose Indicate, by `message()`, the current state of calculation.
#' @returns A `PKNCAresults` object.
#' @seealso [PKNCAdata()], [PKNCA.options()], [summary.PKNCAresults()],
#'   [as.data.frame.PKNCAresults()], [exclude()]
#' @export
pk.nca <- function(data, verbose=FALSE) {
  assert_PKNCAdata(data)
  results <- data.frame()
  if (nrow(data$intervals) > 0) {
    if (verbose) message("Setting up options")
    # Merge the options into the default options.
    tmp_options <- PKNCA.options()
    tmp_options[names(data$options)] <- data$options
    data$options <- tmp_options
    splitdata <- full_join_PKNCAdata(data)
    group_info <-
      splitdata[
        ,
        setdiff(names(splitdata), c("data_conc", "data_sparse_conc", "data_dose", "data_intervals")),
        drop=FALSE
      ]
    # Calculate the results
    if (verbose) message("Starting dense PK NCA calculations.")
    results_dense <-
      purrr::pmap(
        .l = list(
          data_conc = splitdata$data_conc,
          data_dose = splitdata$data_dose,
          data_intervals = splitdata$data_intervals
        ),
        .f = pk.nca.intervals,
        options = data$options,
        impute = data$impute,
        verbose = verbose,
        sparse = FALSE,
        .progress = data$options$progress
      )
    if (verbose) message("Combining completed dense PK calculation results.")
    results <- pk_nca_result_to_df(group_info, results_dense)
    if (is_sparse_pk(data)) {
      if (verbose) message("Starting sparse PK NCA calculations.")
      results_sparse <-
        purrr::pmap(
          .l=list(
            data_conc=splitdata$data_sparse_conc,
            data_dose=splitdata$data_dose,
            data_intervals=splitdata$data_intervals
          ),
          .f=pk.nca.intervals,
          options=data$options,
          impute=data$impute,
          verbose=verbose,
          sparse=TRUE
        )
      if (verbose) message("Combining completed sparse PK calculation results.")
      results <-
        dplyr::bind_rows(
          results,
          pk_nca_result_to_df(group_info, results_sparse)
        )
    }
  }
  PKNCAresults(
    result=results,
    data=data,
    exclude="exclude"
  )
}

#' Convert the grouping info and list of results for each group into a results
#' data.frame
#'
#' @param group_info A data.frame of grouping columns
#' @param result A list of data.frames with the results from NCA parameter
#'   calculations
#' @return A data.frame with group_info and result combined, warnings filtered
#'   out, and results unnested.
#' @keywords Internal
pk_nca_result_to_df <- function(group_info, result) {
  ret <- group_info
  ret$data_result <- result
  # Gather, report, and remove warnings
  mask_warning <- vapply(X=ret$data_result, inherits, what="warning", TRUE)
  ret_warnings <- ret[mask_warning, ]
  if (nrow(ret_warnings) > 0) {
    group_names <- setdiff(names(ret_warnings), "data_result")
    # Tell the user where the warning comes from
    warning_preamble <-
      do.call(
        what=paste,
        args=
          append(
            lapply(
              X=group_names,
              FUN=function(x) paste(x, ret_warnings[[x]], sep="=")
            ),
            list(sep="; ")
          )
      )
    invisible(lapply(
      X=seq_along(warning_preamble),
      FUN=function(idx) {
        warning_prep <- ret_warnings$data_result[[idx]]
        warning_prep$message <- paste(warning_preamble[idx], warning_prep$message, sep=": ")
        warning(warning_prep)
      }
    ))
  }
  ret_nowarning <- ret[!mask_warning, ]
  # Generate the outputs
  if (nrow(ret_nowarning) == 0) {
    rlang::warn(
      message = "All results generated warnings or errors; no results generated",
      class = "pknca_all_warnings_no_results"
    )
    results <- data.frame()
  } else {
    results <- tidyr::unnest(ret_nowarning, cols="data_result")
    rownames(results) <- NULL
  }
  results
}

filter_interval <- function(data, start, end, include_na=FALSE, include_end=TRUE) {
  mask_na <- include_na & is.na(data$time)
  mask_keep_start <- start <= data$time
  mask_keep_end <-
    if (include_end) {
      data$time <= end
    } else {
      data$time < end
    }
  mask_time <- mask_keep_start & mask_keep_end
  data[mask_na | mask_time, ]
}

#' Determine if there are any sparse or dense calculations requested within an interval
#'
#' @param interval An interval specification
#' @inheritParams PKNCAconc
#' @return A logical value indicating if the interval requests any sparse (if
#'   `sparse=TRUE`) or dense (if `sparse=FALSE`) calculations.
#' @keywords Internal
any_sparse_dense_in_interval <- function(interval, sparse) {
  all_intervals <- get.interval.cols()
  interval_subset <- interval[, names(interval) %in% names(all_intervals)]
  requested <- vapply(X = interval_subset, FUN = isTRUE, FUN.VALUE = TRUE)
  # Extract if the parameters to be calculated (`names(requested[requested])`)
  # are sparse, and compare that to if the request is for sparse or dense
  any(
    vapply(
      X=all_intervals[names(requested[requested])],
      FUN="[[",
      "sparse",
      FUN.VALUE = TRUE
    ) %in% sparse
  )
}

# Subset data down to just the times of interest and then pass it
# further to the calculation routines.
#
# This is simply a helper for pk.nca
#' Compute NCA for multiple intervals
#'
#' @param data_conc A data.frame or tibble with standardized column names as
#'   output from `prepare_PKNCAconc()`
#' @param data_dose A data.frame or tibble with standardized column names as
#'   output from `prepare_PKNCAdose()`
#' @param data_intervals A data.frame or tibble with standardized column names
#'   as output from `prepare_PKNCAintervals()`
#' @param impute The column name in `data_intervals` to use for imputation
#' @inheritParams PKNCAdata
#' @inheritParams pk.nca
#' @inheritParams pk.nca.interval
#' @return A data.frame with all NCA results
pk.nca.intervals <- function(data_conc, data_dose, data_intervals, sparse,
                             options, impute, verbose=FALSE) {
  if (is.null(data_conc) || (nrow(data_conc) == 0)) {
    # No concentration data; potentially placebo data
    return(rlang::warning_cnd(class="pknca_no_conc_data", message="No concentration data"))
  } else if (is.null(data_intervals) || (nrow(data_intervals) == 0)) {
    # No intervals; potentially placebo data
    return(rlang::warning_cnd(class="pknca_no_intervals", message="No intervals for data"))
  }
  ret <- data.frame()
  for (i in seq_len(nrow(data_intervals))) {
    current_interval <- data_intervals[i, , drop=FALSE]
    has_calc_sparse_dense <- any_sparse_dense_in_interval(current_interval, sparse=sparse)
    # Choose only times between the start and end.
    conc_data_interval <- filter_interval(data_conc, start=data_intervals$start[i], end=data_intervals$end[i])
    # Sort the data in time order
    conc_data_interval <- conc_data_interval[order(conc_data_interval$time),]
    NA_data_dose_ <- data.frame(dose=NA_real_, time=NA_real_, duration=NA_real_, route=NA_real_)
    if (is.null(data_dose) || identical(data_dose, NA)) {
      data_dose <- dose_data_interval <- NA_data_dose_
    } else {
      # include_end=FALSE so that a dose at the end of an interval is not included
      dose_data_interval <-
        filter_interval(
          data_dose,
          start=data_intervals$start[i],
          end=data_intervals$end[i],
          include_na=TRUE, include_end=FALSE
        )
    }
    if (nrow(dose_data_interval) > 0) {
      dose_data_interval <- dose_data_interval[order(dose_data_interval$time),]
    } else {
      # When all data are filtered out
      dose_data_interval <- NA_data_dose_
    }
    # Setup for detailed error reporting in case it's needed
    error_preamble <-
      paste(
        "Error with interval",
        paste(
          c("start", "end"),
          unlist(current_interval[, c("start", "end")]),
          sep="=", collapse=", ")
      )
    if (nrow(conc_data_interval) == 0) {
      warning(paste(error_preamble, "No data for interval", sep=": "))
    } else if (!has_calc_sparse_dense) {
      if (verbose) message("No ", ifelse(sparse, "sparse", "dense"), " calculations requested for an interval")
    } else {
      impute_method <- get_impute_method(intervals = current_interval, impute = impute)
      args <- list(
        # Interval-level data
        conc=conc_data_interval$conc,
        time=conc_data_interval$time,
        volume=conc_data_interval$volume,
        duration.conc=conc_data_interval$duration,
        dose=dose_data_interval$dose,
        time.dose=dose_data_interval$time,
        duration.dose=dose_data_interval$duration,
        route=dose_data_interval$route,
        impute_method=impute_method,
        # Group-level data
        conc.group=data_conc$conc,
        time.group=data_conc$time,
        volume.group=data_conc$volume,
        duration.conc.group=data_conc$duration,
        dose.group=data_dose$dose,
        time.dose.group=data_dose$time,
        duration.dose.group=data_dose$duration,
        route.group=data_dose$route,
        # Generic data
        sparse=sparse,
        interval=current_interval,
        options=options)
      if ("subject" %in% names(conc_data_interval)) {
        args$subject <- conc_data_interval$subject
      }
      if ("include_half.life" %in% names(conc_data_interval)) {
        args$include_half.life <- conc_data_interval$include_half.life
      }
      if ("exclude_half.life" %in% names(conc_data_interval)) {
        args$exclude_half.life <- conc_data_interval$exclude_half.life
      }
      # Try the calculation
      calculated_interval <-
        tryCatch(
          do.call(pk.nca.interval, args),
          error=function(e) {
            e$message <- paste("Please report a bug.\n", error_preamble, e$message, sep=": ") # nocov
            stop(e) # nocov
          }
        )
      # Add all the new data into the output
      ret <-
        rbind(
          ret,
          cbind(
            current_interval[, c("start", "end", options$keep_interval_cols)],
            calculated_interval,
            row.names=NULL
          )
        )
    }
  }
  ret
}

#' Compute all PK parameters for a single concentration-time data set
#'
#' For one subject/time range, compute all available PK parameters. All the
#' internal options should be set by [PKNCA.options()] prior to running.  The
#' only part that changes with a call to this function is the `conc`entration
#' and `time`.
#'
#' @inheritParams assert_conc_time
#' @inheritParams PKNCA.choose.option
#' @param conc.group All concentrations measured for the group
#' @param time.group Time of all concentrations measured for the group
#' @param volume,volume.group The volume (or mass) of the concentration
#'   measurement for the current interval or all data for the group (typically
#'   for urine and fecal measurements)
#' @param duration.conc,duration.conc.group The duration of the concentration
#'   measurement for the current interval or all data for the group (typically
#'   for urine and fecal measurements)
#' @param dose,dose.group Dose amount (may be a scalar or vector) for the
#'   current interval or all data for the group
#' @param time.dose,time.dose.group Time of the dose for the current interval or
#'   all data for the group (must be the same length as `dose` or `dose.group`)
#' @param duration.dose,duration.dose.group The duration of the dose
#'   administration for the current interval or all data for the group
#'   (typically zero for extravascular and intravascular bolus and nonzero for
#'   intravascular infusion)
#' @param route,route.group The route of dosing for the current interval or all
#'   data for the group
#' @param impute_method The method to use for imputation as a character string
#' @param interval One row of an interval definition (see
#'   [check.interval.specification()] for how to define the interval.
#' @param include_half.life An optional boolean vector of the concentration
#'   measurements to include in the half-life calculation. If given, no
#'   half-life point selection will occur.
#' @param exclude_half.life An optional boolean vector of the concentration
#'   measurements to exclude from the half-life calculation.
#' @param subject Subject identifiers (used for sparse calculations)
#' @param sparse Should only sparse calculations be performed (TRUE) or only
#'   dense calculations (FALSE)?
#' @returns A data frame with the start and end time along with all PK
#'   parameters for the `interval`
#'
#' @seealso [check.interval.specification()]
#' @export
pk.nca.interval <- function(conc, time, volume, duration.conc,
                            dose, time.dose, duration.dose, route,
                            conc.group=NULL, time.group=NULL, volume.group=NULL, duration.conc.group=NULL,
                            dose.group=NULL, time.dose.group=NULL, duration.dose.group=NULL, route.group=NULL,
                            impute_method=NA_character_,
                            include_half.life=NULL, exclude_half.life=NULL,
                            subject, sparse, interval, options=list()) {
  if (!is.data.frame(interval)) {
    stop("Please report a bug.  Interval must be a data.frame")
  }
  if (nrow(interval) != 1) {
    stop("Please report a bug.  Interval must be a one-row data.frame")
  }
  if (!all(is.na(impute_method))) {
    impute_funs <- PKNCA_impute_fun_list(impute_method)
    stopifnot(length(impute_funs) == 1)
    impute_data <- data.frame(conc=conc, time=time)
    for (current_fun_nm in impute_funs[[1]]) {
      impute_args <- as.list(impute_data)
      impute_args$start <- interval$start[1]
      impute_args$end <- interval$end[1]
      impute_args$conc.group <- conc.group
      impute_args$time.group <- time.group
      impute_args$options <- options
      impute_data <- do.call(current_fun_nm, args=impute_args)
    }
    conc <- impute_data$conc
    time <- impute_data$time
  }
  # Prepare the return value using SDTM names
  ret <- data.frame(PPTESTCD=NA, PPORRES=NA)[-1,]
  # Determine exactly what needs to be calculated in what order. Start with the
  # interval specification and find any dependencies that are not listed for
  # calculation.  Then loop over the calculations in order confirming what needs
  # to be passed from a previous calculation to a later calculation.
  all_intervals <- get.interval.cols()
  # Set the dose to NA if its length is zero
  if (length(dose) == 0) {
    dose <- NA
    time.dose <- NA
    duration.dose <- NA
  }
  # Make sure that we calculate all of the dependencies.  Do this in
  # reverse order for dependencies of dependencies.
  for (n in rev(names(all_intervals))) {
    if (interval[[n]]) {
      interval[all_intervals[[n]]$depends] <- TRUE
    }
  }
  # Do the calculations
  for (n in names(all_intervals)) {
    request_to_calculate <- as.logical(interval[[n]])
    has_calculation_function <- !is.na(all_intervals[[n]]$FUN)
    is_correct_sparse_dense <- all_intervals[[n]]$sparse == sparse
    if (request_to_calculate & has_calculation_function & is_correct_sparse_dense) {
      call_args <- list()
      exclude_from_argument <- character(0)
      # Prepare to call the function by setting up its arguments.
      # Define the required arguments (arglist), and ignore the "..." argument
      # if it exists.
      arglist <- setdiff(names(formals(get(all_intervals[[n]]$FUN))),
                         "...")
      arglist <- stats::setNames(object=as.list(arglist), arglist)
      arglist[names(all_intervals[[n]]$formalsmap)] <- all_intervals[[n]]$formalsmap
      # Drop arguments that were set to NULL by the formalsmap
      arglist <- arglist[!vapply(X = arglist, FUN = is.null, FUN.VALUE = TRUE)]
      for (arg_formal in names(arglist)) {
        arg_mapped <- arglist[[arg_formal]]
        if (arg_mapped == "conc") {
          call_args[[arg_formal]] <- conc
        } else if (arg_mapped == "time") {
          # Realign the time to be relative to the start of the
          # interval
          call_args[[arg_formal]] <- time - interval$start[1]
        } else if (arg_mapped == "volume") {
          call_args[[arg_formal]] <- volume
        } else if (arg_mapped == "duration.conc") {
          call_args[[arg_formal]] <- duration.conc
        } else if (arg_mapped == "dose") {
          call_args[[arg_formal]] <- dose
        } else if (arg_mapped == "time.dose") {
          # Realign the time to be relative to the start of the
          # interval
          call_args[[arg_formal]] <- time.dose - interval$start[1]
        } else if (arg_mapped == "duration.dose") {
          call_args[[arg_formal]] <- duration.dose
        } else if (arg_mapped == "route") {
          call_args[[arg_formal]] <- route
        } else if (arg_mapped == "conc.group") {
          call_args[[arg_formal]] <- conc.group
        } else if (arg_mapped == "time.group") {
          # Don't realign the time to be relative to the start of the
          # interval
          call_args[[arg_formal]] <- time.group
        } else if (arg_mapped == "volume.group") {
          call_args[[arg_formal]] <- volume.group
        } else if (arg_mapped == "duration.conc.group") {
          call_args[[arg_formal]] <- duration.conc.group
        } else if (arg_mapped == "dose.group") {
          call_args[[arg_formal]] <- dose.group
        } else if (arg_mapped == "time.dose.group") {
          # Realign the time to be relative to the start of the
          # interval
          call_args[[arg_formal]] <- time.dose.group
        } else if (arg_mapped == "duration.dose.group") {
          call_args[[arg_formal]] <- duration.dose.group
        } else if (arg_mapped == "route.group") {
          call_args[[arg_formal]] <- route.group
        } else if (arg_mapped == "subject") {
          call_args[[arg_formal]] <- subject
        } else if (arg_mapped %in% c("start", "end")) {
          # Provide the start and end of the interval if they are requested
          call_args[[arg_formal]] <- interval[[arg_mapped]]
        } else if (arg_mapped == "options") {
          call_args[[arg_formal]] <- options
        } else if (any(mask_arg <- ret$PPTESTCD %in% arg_mapped)) {
          call_args[[arg_formal]] <- ret$PPORRES[mask_arg]
          exclude_from_argument <-
            c(exclude_from_argument, ret$exclude[mask_arg])
        } else if (!is.null(interval[[arg_mapped]])) {
          call_args[[arg_formal]] <- interval[[arg_mapped]]
        } else {
          # Give an error if there is not a default argument.
          if (inherits(formals(get(all_intervals[[n]]$FUN))[[arg_formal]], "name")) {
            arg_text <-
              if (arg_formal == arg_mapped) {
                sprintf("'%s'", arg_formal)
              } else {
                sprintf("'%s' mapped to '%s'", arg_formal, arg_mapped)
              }
            stop(sprintf( # nocov
              "Cannot find argument %s for NCA function '%s'", # nocov
              arg_text, all_intervals[[n]]$FUN) # nocov
            ) # nocov
          }
        }
      }
      # Apply manual inclusion and exclusion
      if (n %in% "half.life") {
        if (!is.null(include_half.life)) {
          call_args$conc <- call_args$conc[include_half.life]
          call_args$time <- call_args$time[include_half.life]
          call_args$manually.selected.points <- TRUE
        } else if (!is.null(exclude_half.life)) {
          call_args$conc <- call_args$conc[!exclude_half.life]
          call_args$time <- call_args$time[!exclude_half.life]
        }
      }
      # Do the calculation
      tmp_result <- do.call(all_intervals[[n]]$FUN, call_args)
      # The handling of the exclude column is documented in the
      # "Writing-Parameter-Functions.Rmd" vignette.  Document any changes to
      # this section of code there.
      exclude_reason <-
        stats::na.omit(c(
          exclude_from_argument, attr(tmp_result, "exclude")
        ))
      exclude_reason <-
        if (identical(attr(tmp_result, "exclude"), "DO NOT EXCLUDE")) {
          NA_character_
        } else if (length(exclude_reason) > 0) {
          paste(exclude_reason, collapse="; ")
        } else {
          NA_character_
        }
      # If the function returns a data frame, save all the returned values,
      # otherwise, save the value returned.
      if (is.data.frame(tmp_result)) {
        # if (uses_units) {
        #   # Convert to mixed_units so that rbind will work
        #   for (nm in names(tmp_result)) {
        #     if (inherits(tmp_result[[nm]], "units")) {
        #       tmp_result[[nm]] <- units::mixed_units(tmp_result[[nm]])
        #     } else {
        #       # unitless
        #       tmp_result[[nm]] <- units::mixed_units(tmp_result[[nm]], "")
        #     }
        #   }
        # }
        tmp_testcd <- names(tmp_result)
        tmp_result <- unlist(tmp_result, use.names=FALSE, recursive=FALSE)
      } else {
        # if (uses_units) {
        #   if (inherits(tmp_result, "units")) {
        #     # I() due to https://github.com/r-quantities/units/issues/309
        #     tmp_result <- I(units::mixed_units(tmp_result))
        #   } else {
        #     # unitless
        #     # I() due to https://github.com/r-quantities/units/issues/309
        #     tmp_result <- I(units::mixed_units(tmp_result, ""))
        #   }
        # }
        tmp_testcd <- n
      }
      single_result <-
        data.frame(
          PPTESTCD=tmp_testcd,
          PPORRES=tmp_result,
          exclude=exclude_reason,
          stringsAsFactors=FALSE
        )
      ret <- rbind(ret, single_result)
    }
  }
  ret
}

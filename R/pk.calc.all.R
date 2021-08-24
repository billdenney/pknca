#' Compute NCA parameters for each interval for each subject.
#'
#' The \code{pk.nca} function computes the NCA parameters from a
#' \code{PKNCAdata} object.  All options for the calculation and input data are
#' set in prior functions (\code{PKNCAconc}, \code{PKNCAdose}, and
#' \code{PKNCAdata}).  Options for calculations are set either in
#' \code{PKNCAdata} or with the current default options in \code{PKNCA.options}.
#'
#' When performing calculations, all time results are relative to the start of
#' the interval.  For example, if an interval starts at 168 hours, ends at 192
#' hours, and and the maximum concentration is at 169 hours,
#' \code{tmax=169-168=1}.
#' 
#' @param data A PKNCAdata object
#' @param verbose Indicate, by \code{message()}, the current state of
#'   calculation.
#' @return A \code{PKNCAresults} object.
#' @seealso \code{\link{PKNCAdata}}, \code{\link{PKNCA.options}},
#'  \code{\link{summary.PKNCAresults}}, \code{\link{as.data.frame.PKNCAresults}},
#'  \code{\link{exclude}}
#' @export
#' @importFrom dplyr bind_rows
#' @importFrom parallel mclapply
#' @importFrom stats as.formula update.formula
#' @importFrom utils capture.output
#' @importFrom purrr pmap
pk.nca <- function(data, verbose=FALSE) {
  if (nrow(data$intervals) == 0) {
    warning("No intervals given; no calculations done.")
    results <- data.frame()
  } else {
    if (verbose) message("Setting up options")
    ## Merge the options into the default options.
    tmp.opt <- PKNCA.options()
    tmp.opt[names(data$options)] <- data$options
    data$options <- tmp.opt
    splitdata <- full_join_PKNCAdata(data)
    ## Calculate the results
    if (verbose) message("Starting NCA calculations.")
    tmp_results <-
      purrr::pmap(
        .l=list(
          data_conc=splitdata$data_conc,
          data_dose=splitdata$data_dose,
          data_intervals=splitdata$data_intervals
        ),
        .f=pk.nca.intervals,
        options=data$options,
        verbose=verbose
      )
    if (verbose) message("Combining completed results.")
    ret_prep <-
      splitdata[
        ,
        setdiff(names(splitdata), c("data_conc", "data_dose", "data_intervals")),
        drop=FALSE
      ]
    ret_prep$data_result <- tmp_results
    # Gather, report, and remove warnings
    mask_warning <- sapply(X=ret_prep$data_result, inherits, what="warning")
    ret_warnings <- ret_prep[mask_warning, ]
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
    ret_nowarning <- ret_prep[!mask_warning, ]
    # Generate the outputs
    if (nrow(ret_nowarning) == 0) {
      warning("All results generated warnings or errors; no results generated")
      results <- data.frame()
    } else {
      results <- tidyr::unnest(ret_nowarning, cols="data_result")
      rownames(results) <- NULL
    }
  }
  PKNCAresults(
    result=results,
    data=data,
    exclude="exclude"
  )
}

## Subset data down to just the times of interest and then pass it
## further to the calculation routines.
##
## This is simply a helper for pk.nca
#' Compute NCA for multiple intervals
#' 
#' @param data_conc A data.frame or tibble with standardized column names as
#'   output from \code{prepare_PKNCAconc()}
#' @param data_dose A data.frame or tibble with standardized column names as
#'   output from \code{prepare_PKNCAdose()}
#' @param data_intervals A data.frame or tibble with standardized column names
#'   as output from \code{prepare_PKNCAintervals()}
#' @inheritParams PKNCAdata
#' @inheritParams pk.nca
#' @return A data.frame with all NCA results
#' @importFrom rlang warning_cnd
pk.nca.intervals <- function(data_conc, data_dose, data_intervals,
                             options, verbose=FALSE) {
  if (is.null(data_conc) || (nrow(data_conc) == 0)) {
    ## No concentration data; potentially placebo data
    return(rlang::warning_cnd(class="PKNCA_no_conc_data", message="No concentration data"))
  } else if (is.null(data_intervals) || (nrow(data_intervals) == 0)) {
    ## No intervals; potentially placebo data
    return(rlang::warning_cnd(class="PKNCA_no_intervals", message="No intervals for data"))
  }
  ret <- data.frame()
  for (i in seq_len(nrow(data_intervals))) {
    ## Choose only times between the start and end.
    mask.keep.conc <-
      (
        data_intervals$start[i] <= data_conc$time &
          data_conc$time <= data_intervals$end[i]
      )
    conc_data_interval <- data_conc[mask.keep.conc,]
    NA_data_dose_ <- data.frame(dose=NA_real_, time=NA_real_, duration=NA_real_, route=NA_real_)
    if (is.null(data_dose) || identical(data_dose, NA)) {
      data_dose <- dose_data_interval <- NA_data_dose_
    } else {
      mask.keep.dose <-
        (
          is.na(data_dose$time) |
            (data_intervals$start[i] <= data_dose$time &
               data_dose$time < data_intervals$end[i])
        )
      dose_data_interval <- data_dose[mask.keep.dose,]
    }
    ## Sort the data in time order
    conc_data_interval <- conc_data_interval[order(conc_data_interval$time),]
    if (nrow(dose_data_interval) > 0) {
      dose_data_interval <- dose_data_interval[order(dose_data_interval$time),]
    } else {
      # When all data are filtered out
      dose_data_interval <- NA_data_dose_
    }
    ## Setup for detailed error reporting in case it's needed
    error.preamble <-
      paste(
        "Error with interval",
        paste(
          c("start", "end"),
          unlist(data_intervals[i, c("start", "end")]),
          sep="=", collapse=", ")
      )
    if (nrow(conc_data_interval) == 0) {
      warning(paste(error.preamble, "No data for interval", sep=": "))
    } else {
      tryCatch(
        {
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
            interval=data_intervals[i, , drop=FALSE],
            options=options)
          if ("include_half.life" %in% names(conc_data_interval)) {
            args$include_half.life <- conc_data_interval$include_half.life
          }
          if ("exclude_half.life" %in% names(conc_data_interval)) {
            args$exclude_half.life <- conc_data_interval$exclude_half.life
          }
          ## Try the calculation
          calculated.interval <- do.call(pk.nca.interval, args)
        },
        error=function(e) {
          e$message <- paste(error.preamble, e$message, sep=": ")
          stop(e)
        })
      ## Add all the new data into the output
      ret <-
        rbind(
          ret,
          cbind(
            data_intervals[i, c("start", "end")],
            calculated.interval,
            row.names=NULL
          )
        )
    }
  }
  ret
}

#' Compute all PK parameters for a single concentration-time data set
#'
#' For one subject/time range, compute all available PK parameters. All
#' the internal options should be set by \code{\link{PKNCA.options}}
#' prior to running.  The only part that changes with a call to this
#' function is the \code{conc}entration and \code{time}.
#'
#' @param conc,conc.group Concentration measured for the current
#'   interval or all data for the group
#' @param time,time.group Time of concentration measurement for the
#'   current interval or all data for the group
#' @param volume,volume.group The volume (or mass) of the concentration
#'   measurement for the current interval or all data for the group
#'   (typically for urine and fecal measurements)
#' @param duration.conc,duration.conc.group The duration of the
#'   concentration measurement for the current interval or all data for
#'   the group (typically for urine and fecal measurements)
#' @param dose,dose.group Dose amount (may be a scalar or vector) for
#'   the current interval or all data for the group
#' @param time.dose,time.dose.group Time of the dose for the current
#'   interval or all data for the group (must be the same length as
#'   \code{dose} or \code{dose.group})
#' @param duration.dose,duration.dose.group The duration of the dose
#'   administration for the current interval or all data for the group
#'   (typically zero for extravascular and intravascular bolus and
#'   nonzero for intravascular infusion)
#' @param route,route.group The route of dosing for the current interval
#'   or all data for the group
#' @param interval One row of an interval definition (see
#'   \code{\link{check.interval.specification}} for how to define the
#'   interval.
#' @param include_half.life An optional boolean vector of the
#'   concentration measurements to include in the half-life calculation.
#'   If given, no half-life point selection will occur.
#' @param exclude_half.life An optional boolean vector of the
#'   concentration measurements to exclude from the half-life
#'   calculation.
#' @param options List of changes to the default
#'   \code{\link{PKNCA.options}} for calculations.
#' @return A data frame with the start and end time along with all PK
#'   parameters for the \code{interval}
#'
#' @seealso \code{\link{check.interval.specification}}
#' @importFrom stats na.omit setNames
#' @export
pk.nca.interval <- function(conc, time, volume, duration.conc,
                            dose, time.dose, duration.dose, route,
                            conc.group=NULL, time.group=NULL, volume.group=NULL, duration.conc.group=NULL,
                            dose.group=NULL, time.dose.group=NULL, duration.dose.group=NULL, route.group=NULL,
                            include_half.life=NULL, exclude_half.life=NULL,
                            interval, options=list()) {
  if (!is.data.frame(interval)) {
    stop("interval must be a data.frame")
  }
  if (nrow(interval) != 1)
    stop("interval must be a one-row data.frame")
  ## Prepare the return value using SDTM names
  ret <- data.frame(PPTESTCD=NA, PPORRES=NA)[-1,]
  ## Determine exactly what needs to be calculated in what order.
  ## Start with the interval specification and find any dependencies
  ## that are not listed for calculation.  Then loop over the
  ## calculations in order confirming what needs to be passed from a
  ## previous calculation to a later calculation.
  all.intervals <- get.interval.cols()
  ## Set the dose to NA if its length is zero
  if (length(dose) == 0) {
    dose <- NA
    time.dose <- NA
    duration.dose <- NA
  }
  ## Make sure that we calculate all of the dependencies.  Do this in
  ## reverse order for dependencies of dependencies.
  for (n in rev(names(all.intervals))) {
    if (interval[[1,n]]) {
      for (deps in all.intervals[[n]]$depends) {
        interval[1,deps] <- TRUE
      }
    }
  }
  ## Do the calculations
  for (n in names(all.intervals))
    if (interval[[1,n]] & !is.na(all.intervals[[n]]$FUN)) {
      call.args <- list()
      exclude_from_argument <- character(0)
      ## Prepare to call the function by setting up its arguments.
      ## Ignore the "..." argument if it exists.
      arglist <- setdiff(names(formals(get(all.intervals[[n]]$FUN))),
                         "...")
      arglist <- stats::setNames(object=as.list(arglist), arglist)
      arglist[names(all.intervals[[n]]$formalsmap)] <- all.intervals[[n]]$formalsmap
      # Drop arguments that were set to NULL by the formalsmap
      arglist <- arglist[!sapply(arglist, is.null)]
      for (arg_formal in names(arglist)) {
        arg_mapped <- arglist[[arg_formal]]
        if (arg_mapped == "conc") {
          call.args[[arg_formal]] <- conc
        } else if (arg_mapped == "time") {
          ## Realign the time to be relative to the start of the
          ## interval
          call.args[[arg_formal]] <- time - interval$start[1]
        } else if (arg_mapped == "volume") {
          call.args[[arg_formal]] <- volume
        } else if (arg_mapped == "duration.conc") {
          call.args[[arg_formal]] <- duration.conc
        } else if (arg_mapped == "dose") {
          call.args[[arg_formal]] <- dose
        } else if (arg_mapped == "time.dose") {
          ## Realign the time to be relative to the start of the
          ## interval
          call.args[[arg_formal]] <- time.dose - interval$start[1]
        } else if (arg_mapped == "duration.dose") {
          call.args[[arg_formal]] <- duration.dose
        } else if (arg_mapped == "route") {
          call.args[[arg_formal]] <- route
        } else if (arg_mapped == "conc.group") {
          call.args[[arg_formal]] <- conc.group
        } else if (arg_mapped == "time.group") {
          ## Realign the time to be relative to the start of the
          ## interval
          call.args[[arg_formal]] <- time.group
        } else if (arg_mapped == "volume.group") {
          call.args[[arg_formal]] <- volume.group
        } else if (arg_mapped == "duration.conc.group") {
          call.args[[arg_formal]] <- duration.conc.group
        } else if (arg_mapped == "dose.group") {
          call.args[[arg_formal]] <- dose.group
        } else if (arg_mapped == "time.dose.group") {
          ## Realign the time to be relative to the start of the
          ## interval
          call.args[[arg_formal]] <- time.dose.group
        } else if (arg_mapped == "duration.dose.group") {
          call.args[[arg_formal]] <- duration.dose.group
        } else if (arg_mapped == "route.group") {
          call.args[[arg_formal]] <- route.group
        } else if (arg_mapped %in% c("start", "end")) {
          ## Provide the start and end of the interval if they are requested
          call.args[[arg_formal]] <- interval[[arg_mapped]]
        } else if (arg_mapped == "options") {
          call.args[[arg_formal]] <- options
        } else if (any(mask.arg <- ret$PPTESTCD %in% arg_mapped)) {
          call.args[[arg_formal]] <- ret$PPORRES[mask.arg]
          exclude_from_argument <-
            c(exclude_from_argument, ret$exclude[mask.arg])
        } else {
          ## Give an error if there is not a default argument.
          ## FIXME: checking if the class is a name isn't perfect.  
          if (class(formals(get(all.intervals[[n]]$FUN))[[arg_formal]]) == "name") {
            arg_text <-
              if (arg_formal == arg_mapped) {
                sprintf("'%s'", arg_formal)
              } else {
                sprintf("'%s' mapped to '%s'", arg_formal, arg_mapped)
              }
            stop(sprintf(
              "Cannot find argument %s for NCA function '%s'",
              arg_text, all.intervals[[n]]$FUN))
          }
        }
      }
      # Apply manual inclusion and exclusion
      if (n %in% "half.life") {
        if (!is.null(include_half.life)) {
          call.args$conc <- call.args$conc[include_half.life]
          call.args$time <- call.args$time[include_half.life]
          call.args$manually.selected.points <- TRUE
        } else if (!is.null(exclude_half.life)) {
          call.args$conc <- call.args$conc[!exclude_half.life]
          call.args$time <- call.args$time[!exclude_half.life]
        }
      }
      # Do the calculation
      tmp.result <- do.call(all.intervals[[n]]$FUN, call.args)
      # The handling of the exclude column is documented in the
      # "Writing-Parameter-Functions.Rmd" vignette.  Document any changes to
      # this section of code there.
      exclude_reason <-
        stats::na.omit(c(
          exclude_from_argument, attr(tmp.result, "exclude")
        ))
      exclude_reason <-
        if (identical(attr(tmp.result, "exclude"), "DO NOT EXCLUDE")) {
          NA_character_
        } else if (length(exclude_reason) > 0) {
          paste(exclude_reason, collapse="; ")
        } else {
          NA_character_
        }
      ## If the function returns a data frame, save all the returned
      ## values, otherwise, save the value returned.
      if (is.data.frame(tmp.result)) {
        ret <- rbind(ret,
                     data.frame(PPTESTCD=names(tmp.result),
                                PPORRES=unlist(tmp.result, use.names=FALSE),
                                exclude=exclude_reason,
                                stringsAsFactors=FALSE))
      } else {
        ret <- rbind(ret,
                     data.frame(PPTESTCD=n,
                                PPORRES=tmp.result,
                                exclude=exclude_reason,
                                stringsAsFactors=FALSE))
      }
    }
  ret
}

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
pk.nca <- function(data, verbose=FALSE) {
  if (nrow(data$intervals) == 0) {
    warning("No intervals given; no calculations done.")
    results <- data.frame()
  } else {
    if (verbose) message("Setting up dosing information")
    if (identical(NA, data$dose)) {
      # If no dose information is given, add NULL dose information.
      message("No dose information provided, calculations requiring dose will return NA.")
      tmp.dose.data <- unique(getGroups(data$conc))
      data$dose <-
        PKNCAdose(
          data=tmp.dose.data,
          formula=stats::as.formula(
            paste0(
              ".~.|",
              paste(names(tmp.dose.data), collapse="+")
            )
          )
        )
    }
    if (identical(all.vars(parseFormula(data$dose)$lhs), character())) {
      ## If dose amount is not given, give a false dose column with all 
      ## NAs to simplify dose handling in subsequent steps.
      col.dose <- paste0(max(names(data$dose$data)), "X")
      data$dose$data[,col.dose] <- NA
      data$dose$formula <-
        stats::update.formula(data$dose$formula, paste0(col.dose, "~."))
    }
    if (verbose) message("Setting up options")
    ## Merge the options into the default options.
    tmp.opt <- PKNCA.options()
    tmp.opt[names(data$options)] <- data$options
    data$options <- tmp.opt
    splitdata <- split.PKNCAdata(data)
    # Calculations will only be performed when an interval is requested
    if (verbose) message("Checking that intervals have concentration and dose data.")
    mask_has_interval <-
      sapply(splitdata,
             FUN=function(x) {
               (!is.null(x$intervals)) &&
                 (nrow(x$intervals) > 0)
             })
    mask_has_dose <-
      sapply(splitdata,
             FUN=function(x) {
               !is.null(x$dose)
             })
    mask_has_conc <-
      sapply(splitdata,
             FUN=function(x) {
               !is.null(x$conc)
             })
    if (any(!mask_has_interval)) {
      message(sum(!mask_has_interval), " groups have no interval calculations requested.")
    }
    if (any(mask_missing_dose <- !mask_has_dose & mask_has_conc & mask_has_interval)) {
      missing_groups <- list()
      for (current_idx in which(mask_missing_dose)) {
        tmp_dose_data <- unique(getGroups(splitdata[[current_idx]]$conc))
        splitdata[[current_idx]]$dose <-
          PKNCAdose(
            data=tmp_dose_data,
            formula=stats::as.formula(
              paste0(
                ".~.|",
                paste(names(tmp_dose_data), collapse="+"))
            )
          )
        missing_groups <- append(missing_groups, tmp_dose_data)
      }
      warning("The following intervals are missing dosing data:\n",
              paste(
                capture.output(
                  print(as.data.frame(dplyr::bind_rows(missing_groups)))),
                collapse="\n"))
    }
    ## Calculate the results
    if (verbose) message("Starting NCA calculations.")
    tmp.results <- list()
    tmp.results[mask_has_interval] <-
      parallel::mclapply(
        X=splitdata[mask_has_interval],
        FUN=pk.nca.intervals,
        options=data$options,
        verbose=verbose
      )
    if (verbose) message("Combining completed results.")
    ## Put the group parameters with the results
    for (i in seq_len(length(tmp.results))) {
      ## If no calculations were performed, the results are NULL.
      if (!is.null(tmp.results[[i]])) {
        ## If calculations were performed, the results are non-NULL, add
        ## the grouping information to the results, if applicable.
        keep.group.names <- setdiff(names(attr(splitdata, "groupid")),
                                    names(tmp.results[[i]]))
        if (length(keep.group.names) > 0) {
          tmp.results[[i]] <-
            cbind(attr(splitdata, "groupid")[i,keep.group.names,drop=FALSE],
                  tmp.results[[i]])
        }
      }
    }
    ## Generate the outputs
    results <- dplyr::bind_rows(tmp.results)
    rownames(results) <- NULL
  }
  PKNCAresults(result=results,
               data=data,
               exclude="exclude")
}

## Subset data down to just the times of interest and then pass it
## further to the calculation routines.
##
## This is simply a helper for pk.nca
pk.nca.intervals <- function(conc.dose, intervals, options, verbose=FALSE) {
  if (is.null(conc.dose$conc)) {
    ## No data; potentially placebo data (the warning would have
    ## already been generated from making the PKNCAdata object.
    return(NULL)
  }
  ret <- data.frame()
  all.intervals <- conc.dose$intervals
  pformula.conc <- parseFormula(conc.dose$conc)
  pformula.dose <- parseFormula(conc.dose$dose)
  shared.names <- all.vars(pformula.conc$groups)
  ## Column names to use
  col.conc <- all.vars(pformula.conc$lhs)
  col.time <- all.vars(pformula.conc$rhs)
  col.volume <- conc.dose$conc$columns$volume
  col.duration.conc <- conc.dose$conc$columns$duration
  col.include_half.life <- conc.dose$conc$columns$include_half.life
  col.exclude_half.life <- conc.dose$conc$columns$exclude_half.life
  col.dose <- all.vars(pformula.dose$lhs)
  col.time.dose <- all.vars(pformula.dose$rhs)
  col.duration.dose <- conc.dose$dose$columns$duration
  col.route <- conc.dose$dose$columns$route
  # Insert NA doses and dose times if they are not given
  if (!(col.dose %in% names(conc.dose$dose$data))) {
    col.dose <- paste0(max(names(conc.dose$dose$data)), "X")
    conc.dose$dose$data[[col.dose]] <- NA
  }
  if (!(col.time.dose %in% names(conc.dose$dose$data))) {
    col.time.dose <- paste0(max(names(conc.dose$dose$data)), "X")
    conc.dose$dose$data[[col.time.dose]] <- NA
  }
  # Exclude data once at the beginning
  conc_data_all <-
    conc.dose$conc$data[
      # Remove rows to be excluded from all calculations
      is.na(normalize_exclude(conc.dose$conc$data[[conc.dose$conc$exclude]])),,
      drop=FALSE]
  dose_data_all <-
    conc.dose$dose$data[
      # Remove rows to be excluded from all calculations
      is.na(normalize_exclude(conc.dose$dose$data[[conc.dose$dose$exclude]])),,
      drop=FALSE]
  for (i in seq_len(nrow(all.intervals))) {
    ## Subset the data down to the group of current interest, and make 
    ## the first column of each the dependent variable and the second 
    ## column the independent variable.
    conc_data_group <-
      merge(conc_data_all,
        all.intervals[
          i,
          intersect(shared.names, names(all.intervals)),
          drop=FALSE])[,
                       c(col.conc,
                         col.time,
                         col.include_half.life,
                         col.exclude_half.life,
                         col.volume,
                         col.duration.conc)]
    dose_data_group <-
      merge(dose_data_all,
        all.intervals[
          i,
          intersect(shared.names, names(all.intervals)),
          drop=FALSE])[,
                       c(col.dose,
                         col.time.dose,
                         col.duration.dose,
                         col.route)]
    ## Choose only times between the start and end.
    mask.keep.conc <- (all.intervals$start[i] <= conc_data_group[[col.time]] &
                         conc_data_group[[col.time]] <= all.intervals$end[i])
    conc_data_interval <- conc_data_group[mask.keep.conc,]
    mask.keep.dose <- (is.na(dose_data_group[,col.time.dose]) |
                         (all.intervals$start[i] <= dose_data_group[[col.time.dose]] &
                            dose_data_group[[col.time.dose]] < all.intervals$end[i]))
    dose_data_interval <- dose_data_group[mask.keep.dose,]
    ## Sort the data in time order
    conc_data_interval <- conc_data_interval[order(conc_data_interval[[col.time]]),]
    dose_data_interval <- dose_data_interval[order(dose_data_interval[[col.time.dose]]),]
    ## Setup for detailed error reporting in case it's needed
    error.preamble <-
      paste("Error with interval",
            paste(c(shared.names, c("start", "end")),
                  c(unlist(conc_data_all[1,shared.names]),
                    unlist(all.intervals[i,c("start", "end")])),
                  sep="=", collapse=", "))
    if (nrow(conc_data_interval) == 0) {
      warning(paste(error.preamble, "No data for interval", sep=": "))
    } else {
      tryCatch(
        {
          args <- list(
            # Interval-level data
            conc=conc_data_interval[[col.conc]],
            time=conc_data_interval[[col.time]],
            volume=conc_data_interval[[col.volume]],
            duration.conc=conc_data_interval[[col.duration.conc]],
            dose=dose_data_interval[[col.dose]],
            time.dose=dose_data_interval[[col.time.dose]],
            duration.dose=dose_data_interval[[col.duration.dose]],
            route=dose_data_interval[[col.route]],
            # Group-level data
            conc.group=conc_data_group[[col.conc]],
            time.group=conc_data_group[[col.time]],
            volume.group=conc_data_group[[col.volume]],
            duration.conc.group=conc_data_group[[col.duration.conc]],
            dose.group=dose_data_group[[col.dose]],
            time.dose.group=dose_data_group[[col.time.dose]],
            duration.dose.group=dose_data_group[[col.duration.dose]],
            route.group=dose_data_group[[col.route]],
            # Generic data
            interval=all.intervals[i, , drop=FALSE],
            options=options)
          if (!is.null(col.include_half.life)) {
            args$include_half.life <- conc_data_interval[[col.include_half.life]]
          }
          if (!is.null(col.exclude_half.life)) {
            args$exclude_half.life <- conc_data_interval[[col.exclude_half.life]]
          }
          ## Try the calculation
          calculated.interval <- do.call(pk.nca.interval, args)
        },
        error=function(e) {
          e$message <- paste(error.preamble, e$message, sep=": ")
          stop(e)
        })
      ## Add all the new data into the output
      ret <- rbind(ret,
                   cbind(all.intervals[i,c("start", "end")],
                         conc_data_all[1, shared.names, drop=FALSE],
                         calculated.interval,
                         row.names=NULL))
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
          call.args[[arg_formal]] <- interval[1,arg_mapped]
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

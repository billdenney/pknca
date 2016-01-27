#' Compute NCA parameters for each interval for each subject.
#'
#' The computations assume that all calculation options are set in
#' \code{\link{PKNCA.options}}.
#'
#' @param data A PKNCAdata object
#' @return A data frame with a row for each interval for each grouping
#' in the concentration data.
#' @seealso \code{\link{PKNCAdata}}, \code{\link{PKNCA.options}}
#' @export
pk.nca <- function(data) {
  ## Separate the concentration and dosing data by their 
  conc.split <- doBy::splitBy(parseFormula(data$conc)$groupFormula,
                              data=stats::model.frame(data$conc))
  ## If dose is not given, give a false dose column with all NAs to
  ## simplify dose handling in subsequent steps.
  if (identical(all.vars(parseFormula(data$dose)$lhs), character())) {
    col.dose <- paste0(max(names(data$dose$data)), "X")
    data$dose$data[,col.dose] <- NA
    data$dose$formula <-
      stats::update.formula(data$dose$formula, paste0(col.dose, "~."))
  }
  dose.split <- doBy::splitBy(parseFormula(data$dose)$groupFormula,
                              data=stats::model.frame(data$dose))
  conc.dose <- merge(conc=conc.split, dose=dose.split)
  if (nrow(data$intervals) == 0) {
    warning("No intervals given; no calculations done.")
    results <- data.frame()
  } else {
    ## Merge the options into the default options.
    tmp.opt <- PKNCA.options()
    tmp.opt[names(options)] <- data$options
    data$options <- tmp.opt
    ## Calculate the results
    tmp.results <-
      parallel::mclapply(X=conc.dose,
                         FUN=pk.nca.intervals,
                         intervals=data$intervals,
                         options=data$options)
    ## Put the group parameters with the results
    for (i in seq_len(length(tmp.results))) {
      ## If no calculations were performed, the results are NULL.
      ## That will move through rbind.fill without adding rows, so no
      ## action is required.
      if (!is.null(tmp.results[[i]])) {
        ## If calculations were performed, the results are non-NULL,
        ## so add the grouping information to the results, if
        ## applicable.
        keep.group.names <- setdiff(names(attr(conc.split, "groupid")),
                                    names(tmp.results[[i]]))
        if (length(keep.group.names) > 0) {
          tmp.results[[i]] <-
            cbind(attr(conc.split, "groupid")[i,keep.group.names,drop=FALSE],
                  tmp.results[[i]])
        }
      }
    }
    ## Generate the outputs
    results <- do.call(plyr::rbind.fill, tmp.results)
    rownames(results) <- NULL
  }
  PKNCAresults(result=results,
               data=data,
               provenance=list(
                 hash=digest::digest(list(results, data)),
                 sessionInfo=utils::sessionInfo(),
                 datetime=Sys.time(),
                 sysInfo=Sys.info()))
}

## Subset data down to just the times of interest and then pass it
## further to the calculation routines.
##
## This is simply a helper for pk.nca
pk.nca.intervals <- function(conc.dose, intervals, options) {
  if (is.null(conc.dose$conc)) {
    ## No data; potentially placebo data (the warning would have
    ## already been generated from making the PKNCAdata object.
    return(NULL)
  }
  ## Merge the intervals with the group columns from the concentration
  ## data
  if (ncol(conc.dose$conc) >= 3) {
    ## Subset the intervals down to the intervals for the current group.
    shared.names <- names(conc.dose$conc)[3:ncol(conc.dose$conc)]
    all.intervals <- merge(unique(conc.dose$conc[, shared.names, drop=FALSE]),
                           intervals)
  } else {
    ## Likely single-subject data
    all.intervals <- intervals
    shared.names <- c()
  }
  ret <- data.frame()
  ## Column names to use
  col.conc <- names(conc.dose$conc)[1]
  col.time <- names(conc.dose$conc)[2]
  col.dose <- names(conc.dose$dose)[1]
  col.time.dose <- names(conc.dose$dose)[2]
  for (i in seq_len(nrow(all.intervals))) {
    ## Subset the data down to the group of current interest
    tmpconcdata <-
      merge(conc.dose$conc,
            all.intervals[i,shared.names])[,c(col.conc, col.time)]
    tmpdosedata <-
      merge(conc.dose$dose,
            all.intervals[i,])[,c(col.dose, col.time.dose)]
    ## Choose only times between the start and end.
    mask.keep <- (all.intervals$start[i] <= tmpconcdata[,col.time] &
                  tmpconcdata[,col.time] <= all.intervals$end[i])
    tmpconcdata <- tmpconcdata[mask.keep,]
    if (nrow(tmpconcdata) == 0) {
      ## TODO: Improve this error message with additional information
      ## on the specific interval that has no data.
      warning("No data for interval")
    } else {
      calculated.interval <-
        pk.nca.interval(conc=tmpconcdata[,col.conc],
                        time=tmpconcdata[,col.time],
                        dose=tmpdosedata[,col.dose],
                        time.dose=tmpdosedata[,col.time.dose],
                        interval=all.intervals[i,],
                        options=options)
      ## Add all the new data into the output
      ret <- rbind(ret,
                   cbind(all.intervals[i,c("start", "end", shared.names)],
                         calculated.interval,
                         row.names=NULL))
    }
  }
  ret
}

#' Compute all PK parameters for a single concentration-time data set
#'
#' For one subject/time range, compute all available PK parameters.
#' All the internal options should be set by
#' \code{\link{PKNCA.options}} prior to running.  The only part that
#' changes with a call to this function is the \code{conc}entration
#' and \code{time}.
#'
#' @param conc Concentration measured
#' @param time Time of concentration measurement
#' @param dose Dose amount (may be a scalar or vector)
#' @param time.dose Time of the dose (must be the same length as
#'   \code{dose})
#' @param interval One row of an interval definition (see
#'   \code{\link{check.interval.specification}} for how to define the
#'   interval.
#' @param options List of changes to the default
#'   \code{\link{PKNCA.options}} for calculations.
#' @return A data frame with the start and end time along with all PK
#'   parameters for the \code{interval}
#' 
#' @seealso \code{\link{check.interval.specification}}
#' @export
pk.nca.interval <- function(conc, time,
                            dose, time.dose,
                            interval, options=list()) {
  if (!is.data.frame(interval))
    stop("interval must be a data.frame")
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
  ## Make sure that we calculate all of the dependencies.  Do this in
  ## reverse order for dependencies of dependencies.
  for (n in rev(names(all.intervals)))
    if (interval[1,n])
      for (deps in all.intervals[[n]]$depends)
        interval[1,deps] <- TRUE
  ## Do the calculations
  for (n in names(all.intervals))
    if (interval[1,n] & !is.na(all.intervals[[n]]$FUN)) {
      call.args <- list()
      ## Prepare to call the function by setting up its arguments.
      ## Ignore the "..." argument if it exists.
      for (arg in setdiff(names(formals(get(all.intervals[[n]]$FUN))),
                          "...")) {
        if (arg == "conc") {
          call.args[[arg]] <- conc
        } else if (arg == "time") {
          ## Realign the time to be relative to the start of the
          ## interval
          call.args[[arg]] <- time - interval$start[1]
        } else if (arg == "dose") {
          call.args[[arg]] <- dose
        } else if (arg == "time.dose") {
          ## Realign the time to be relative to the start of the
          ## interval
          call.args[[arg]] <- time.dose - interval$start[1]
        } else if (arg == "options") {
          call.args[[arg]] <- options
        } else if (any(mask.arg <- ret$PPTESTCD %in% arg)) {
          call.args[[arg]] <- ret$PPORRES[mask.arg]
        } else {
          ## Give an error if there is not a default argument.
          ## FIXME: checking if the class is a name isn't perfect.  
          if (class(formals(get(all.intervals[[n]]$FUN))[[arg]]) == "name")
            stop(sprintf(
              "Cannot find argument '%s' for NCA function '%s'",
              arg, all.intervals[[n]]$FUN))
        }
      }
      tmp.result <- do.call(all.intervals[[n]]$FUN, call.args)
      ## If the function returns a data frame, save all the returned
      ## values, otherwise, save the value returned.
      if (is.data.frame(tmp.result)) {
        ret <- rbind(ret, data.frame(PPTESTCD=names(tmp.result),
                                     PPORRES=unlist(tmp.result, use.names=FALSE)))
      } else {
        ret <- rbind(ret, data.frame(PPTESTCD=n,
                                     PPORRES=tmp.result))
      }
    }
  ret
}

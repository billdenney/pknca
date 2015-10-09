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
  tmp.data <- doBy::splitBy(parseFormula(data$conc)$groupFormula,
                            data=model.frame(data$conc))
  if (nrow(data$intervals) == 0) {
    warning("No intervals given; no calculations done.")
    results <- data.frame()
  } else {
    tmp.results <-
      parallel::mclapply(X=tmp.data,
                         FUN=pk.nca.intervals,
                         intervals=data$intervals,
                         options=data$options)
    ## Put the group parameters with the results
    for (i in seq_len(length(tmp.results)))
      tmp.results[[i]] <- cbind(attr(tmp.data, "groupid")[i,],
                                tmp.results[[i]])
    ## Generate the outputs
    results <- do.call(plyr::rbind.fill, tmp.results)
    rownames(results) <- NULL
  }
  ## FIXME: Add sessionInfo, username, computer name, date/time
  PKNCAresults(result=results,
               formula=formula(data$conc),
               options=data$options)
}

## Subset data down to just the times of interest and then pass it
## further to the calculation routines.
##
## This is simply a helper for pk.nca
pk.nca.intervals <- function(data, intervals, options) {
  ## Merge the intervals with the group columns from the data
  if (ncol(data) >= 3) {
    ## Subset the intervals down to the intervals for the current group.
    shared.names <- names(data)[3:ncol(data)]
    ret <- merge(unique(data[, shared.names, drop=FALSE]),
                 intervals[,c("start", "end", shared.names)])
  } else {
    ## Likely single-subject data
    ret <- intervals[,c("start", "end")]
    shared.names <- c()
  }
  ## Column names to use
  col.conc <- names(data)[1]
  col.time <- names(data)[2]
  for (i in seq_len(nrow(ret))) {
    ## Subset the data down to the group of current interest
    tmpdata <-
      merge(data, ret[i,shared.names])[,c(col.conc, col.time)]
    ## Choose only times between the start and end.
    mask.keep <- (ret$start[i] <= tmpdata[,col.time] &
                  tmpdata[,col.time] <= ret$end[i])
    tmpdata <- tmpdata[mask.keep,]
    if (nrow(tmpdata) == 0) {
      ## TODO: Improve this error message with additional information
      ## on the specific interval that has no data.
      warning("No data for interval")
    } else {
      calculated.interval <-
        pk.nca.interval(conc=tmpdata[,col.conc],
                        time=tmpdata[,col.time],
                        interval=intervals[i,],
                        options=options)
      ## Add new columns to the return dataset
      for (n in setdiff(names(calculated.interval),
                        names(ret)))
        ret[,n] <- NA
      ## Add all the new data into the output
      ret[i,names(calculated.interval)] <- calculated.interval
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
#' @param interval One row of an interval definition (see
#' \code{\link{check.interval.specification}} for how to define the
#' interval.
#' @param options List of changes to the default
#' \code{\link{PKNCA.options}} for calculations.
#' @return A data frame with the start and end time along with all PK
#' parameters for the \code{interval}
#' 
#' @seealso \code{\link{check.interval.specification}}
#' @export
pk.nca.interval <- function(conc, time, interval, options=list()) {
  if (!is.data.frame(interval))
    stop("interval must be a data.frame")
  if (nrow(interval) != 1)
    stop("interval must be a one-row data.frame")
  ## Prepare the return value
  ret <- data.frame()
  ## Determine exactly what needs to be calculated in what order.
  ## Start with the interval specification and find any dependencies
  ## that are not listed for calculation.  Then loop over the
  ## calculations in order confirming what needs to be passed from a
  ## previous calculation to a later calculation.
  all.intervals <- get.interval.cols()
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
          call.args[[arg]] <- time
        } else if (arg == "options") {
          call.args[[arg]] <- options
        } else if (arg %in% names(ret)) {
          call.args[[arg]] <- ret[1,arg]
        } else {
          ## Give an error if there is not a default argument
          if (formals(get(all.intervals[[n]]$FUN))[[arg]] == "")
            stop(sprintf(
              "Cannot find argument '%s' for NCA function '%s'",
              arg, all.intervals[[n]]$FUN))
        }
      }
      tmp.result <- do.call(all.intervals[[n]]$FUN, call.args)
      ## If the function returns a data frame, save all the returned
      ## values, otherwise, save the value returned.
      if (is.data.frame(tmp.result)) {
        ret <- cbind(ret, tmp.result)
      } else {
        ret[1,n] <- tmp.result
      }
    }
  ret
}

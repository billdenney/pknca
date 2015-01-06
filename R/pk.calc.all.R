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
  tmp.data <- splitBy(parseFormula(data$conc)$groupFormula,
                      data=model.frame(data$conc))
  if (nrow(data$intervals) == 0) {
    warning("No intervals given; no calculations done.")
    results <- data.frame()
  } else {
    tmp.results <- mclapply(X=tmp.data,
                            FUN=pk.nca.intervals,
                            intervals=data$intervals,
                            options=data$options)
    ## Put the group parameters with the results
    for (i in seq_len(length(tmp.results)))
      tmp.results[[i]] <- cbind(attr(tmp.data, "groupid")[i,],
                                tmp.results[[i]])
    ## Generate the outputs
    results <- do.call(rbind.fill, tmp.results)
    rownames(results) <- NULL
  }
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
                 intervals)
  } else {
    ## Likely single-subject data
    ret <- intervals
    shared.names <- c()
  }
  ## Column names to use
  col.conc <- names(data)[1]
  col.time <- names(data)[2]
  ## The half.life column will be filled with the computed half-life.
  ## Change its name
  ret$calculate.half.life <- ret$half.life
  ret$half.life <- NULL
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
                        auc.start=ret$start[i],
                        auc.end=ret$end[i],
                        auc.type=ret$auc.type[i],
                        half.life=ret$calculate.half.life[i])
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
#' @param auc.start The start time for the calculations.
#' @param auc.end The end time for the calculations.
#' @param auc.type The type of AUC to calculate.  See
#' \code{\link{pk.calc.auc}}
#' @param half.life logical. Should the half-life be calculated?
#' @param options List of changes to the default
#' \code{\link{PKNCA.options}} for calculations.
#' @param method The AUC integration method.  See
#' \code{\link{pk.calc.auc}}.
#' @param conc.blq The handling instructions for BLQ concentration
#' values.  See \code{\link{clean.conc.blq}}.
#' @param conc.na The handling instructions for missing concentration
#' values.  See \code{\link{clean.conc.na}}.
#' @param first.tmax The calculation instructions for using the first
#' or last Tmax.  See \code{\link{pk.calc.tmax}}.
#' @return A data frame with all PK parameters
#' 
#' @seealso \code{\link{pk.calc.half.life}}
#' @export
pk.nca.interval <- function(conc, time,
                            auc.start, auc.end, auc.type,
                            half.life,
                            options=list(),
                            method=PKNCA.choose.option("auc.method", options),
                            conc.blq=PKNCA.choose.option("conc.blq", options),
                            conc.na=PKNCA.choose.option("conc.na", options),
                            first.tmax=PKNCA.choose.option("first.tmax", options)) {
  if (half.life | auc.type %in% "AUCinf") {
    ret <- pk.calc.half.life(conc, time - auc.start)
  } else {
    ret <- data.frame(lambda.z=NA)
  }
  if (!is.na(auc.type)) {
    ## Compute the AUCs
    auc.name <- paste("auc", substr(auc.type, 4, 100),
                      0, auc.end - auc.start, sep=".")
    aumc.name <- paste("aumc", substr(auc.type, 4, 100),
                       0, auc.end - auc.start, sep=".")
    lambda.z <- NA
    if ("lambda.z" %in% names(ret))
      lambda.z <- ret$lambda.z
    ## Do the AUC computation
    ret[1,auc.name] <-
      pk.calc.auc(conc, time - auc.start,
                  interval=c(auc.start, auc.end) - auc.start,
                  lambda.z=lambda.z,
                  method=method,
                  conc.blq=conc.blq,
                  conc.na=conc.na)
    ret[1,aumc.name] <-
      pk.calc.aumc(conc, time - auc.start,
                   interval=c(auc.start, auc.end) - auc.start,
                   lambda.z=lambda.z,
                   method=method,
                   conc.blq=conc.blq,
                   conc.na=conc.na)
    ret$cmin <- pk.calc.cmin(conc)
  }
  ret$cmax <- pk.calc.cmax(conc)
  ret$tmax <-
    pk.calc.tmax(conc, time, use.first=first.tmax) - auc.start
  ret$tfirst <- pk.calc.tfirst(conc, time) - auc.start
  ret$tlast <- pk.calc.tlast(conc, time) - auc.start
  ret$clast.obs <- pk.calc.clast.obs(conc, time)
  ret
}

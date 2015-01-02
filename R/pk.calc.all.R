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
                      data=model.frame(data))
  if (nrow(data$intervals) == 0) {
    warning("No intervals given; no calculations done.")
    ret <- data.frame()
  } else {
    tmp.ret <- mclapply(X=tmp.data,
                        FUN=pk.calc.intervals,
                        intervals=data$intervals)
    ## Put the group parameters with the results
    for (i in seq_len(length(tmp.ret)))
      tmp.ret[[i]] <- cbind(attributes(tmp.data, "groupid")[i,],
                            tmp.ret[[i]])
    ## Generate the outputs
    ret <- do.call(rbind, tmp.ret)
  }
  rownames(ret) <- NULL
  class(ret) <- c("PKNCAresults", class(ret))
  ret
}

## Subset data down to just the times of interest and then pass it
## further to the calculation routines.
##
## This is simply a helper for pk.nca
pk.nca.intervals <- function(data, intervals) {
  ## Merge the intervals with the group columns from the data
  if (ncol(data) >= 3) {
    ret <- merge(unique(data[, 3:ncol(data), drop=FALSE]),
                     intervals)
  } else {
    ## Likely single-subject data
    ret <- intervals
  }
  ## Column names to use
  col.conc <- names(data)[1]
  col.time <- names(data)[2]
  shared.names <- intersect(names(interval.data), names(ret))
  ## The half.life column will be filled with the computed half-life.
  ## Change its name
  ret$calculate.half.life <- ret$half.life
  ret$half.life <- NULL
  for (i in seq_len(nrow(ret))) {
    ## Subset the data down to the group of current interest
    tmpdata <-
      merge(data, interval.data[i,shared.names])[,c(col.conc, col.time)]
    mask.keep <- (ret$auc.start[i] <= tmpdata[,col.time] &
                  tmpdata[,col.time] <= ret$auc.end[i])
    tmpdata <- tmpdata[mask.keep,]
    if (nrow(tmpdata) == 0) {
      ## TODO: Improve this error message with additional information
      ## on the specific interval that has no data.
      warning("No data for interval")
    } else {
      calculated.interval <-
        pk.nca.interval(conc=tmpdata[,col.conc],
                        time=tmpdata[,col.time],
                        auc.start=ret$auc.start[i],
                        auc.end=ret$auc.end[i],
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

#' Compute all PK parameters
#'
#' For one subject/time range, compute all available PK parameters.
#' All the internal options should be set by
#' \code{\link{PKNCA.options}} prior to running.  The only part that
#' changes with a call to this function is the \code{conc}entration
#' and \code{time}.
#'
#' @param conc Concentration measured
#' @param time Time of concentration measurement
#' @param auc.intervals A list of intervals to compute AUCs.  See
#' details for a description of the \code{auc.intervals} list.
#' @return A data frame with all PK parameters
#' 
#' @seealso \code{\link{pk.calc.half.life}}
#' @export
pk.nca.interval <- function(conc, time,
                            auc.start, auc.end, auc.type,
                            half.life,
                            method=PKNCA.options("auc.method"),
                            conc.blq=PKNCA.options("conc.blq"),
                            conc.na=PKNCA.options("conc.na"),
                            first.tmax=PKNCA.options("first.tmax")) {
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

#' This differs from pk.nca.interval because it allows many intervals
#' to be calculated for a single individual.
#'
#' @param conc The plasma concentration
#' @param time The time of plasma concentration measurement
#' @param auc.intervals The AUC intervals (see
#' \code{\link{check.auc.specification}})
#' @return A data frame with the NCA parameters calculated.
#' @export
pk.nca.individual <- function(conc, time,
                              auc.intervals=data.frame(),
                              ...) {
  ## Create an empty data frame ready for output data with the right
  ## number of rows.
  ret <- data.frame(start.time=auc.intervals$start,
                    end.time=auc.intervals$end)
  for (i in 1:nrow(auc.intervals)) {
    mask.time <- time >= auc.intervals$start[i]
    if (is.na(auc.intervals$end[i]) |
        is.infinite(auc.intervals$end[i])) {
      ## Do nothing, keep all the other times
    } else {
      mask.time <- (mask.time &
                    time <= auc.intervals$end[i])
    }
    if (!any(mask.time)) {
      warning(sprintf("No times matched AUC interval %d", i))
    } else {
      tmp <- pk.nca.interval(conc[mask.time], time[mask.time],
                             auc.start=auc.intervals$start[i],
                             auc.end=auc.intervals$end[i],
                             auc.type=auc.intervals$auc.type[i],
                             half.life=auc.intervals$half.life[i],
                             ...)
      ## Add the output data into the output data frame
      for (n in names(tmp)) {
        if (!(n %in% names(ret)))
          ret[,n] <- NA
        ret[i,n] <- tmp[1,n]
      }
    }
  }
  ret
}

#' Compute the NCA parameters for many groups (e.g. subjects,
#' treatments, studies, analytes, etc.) simultaneously.
#'
#' @param formula.conc the formula defining which columns in the data
#' correspond to \code{concentration~time|group1+group2+...}
#' @param data.conc the data corresponding with \code{formula.conc}
#' @param formula.dose the formula defining which columns in the data
#' correspond to \code{~time.dose|group1+group2+...}
#' @param data.dose the data corresponding with \code{formula.dose}
#' @param auc.intervals specifications for calculation of AUC and
#' half-life as defined in \code{\link{check.auc.specification}}
#' @param ... Parameters passed to other methods.
#' @return a data.frame with columns for each concentration group
#' (rows change slowest in the first group and fastest in the last
#' group) and for each calculated PK parameter.
pk.nca.groups <- function(...)
  UseMethod("pk.nca.groups")

pk.nca.groups.formula <- function(formula.conc, data.conc,
                                  formula.dose, data.dose,
                                  auc.intervals=data.frame(),
                                  ...) {
  form.parts.conc <- parseGroupFormula(formula.conc,
                                       require.groups=FALSE,
                                       require.two.sided=TRUE)
  form.parts.dose <- parseGroupFormula(formula.dose,
                                       require.groups=FALSE,
                                       require.two.sided=FALSE)
  ## Check that all the formula have the correct setup
  conc.column <- list(conc=all.vars(form.pars.conc$lhs),
                      time=all.vars(form.pars.conc$rhs),
                      groups=all.vars(form.parts.conc$groups))
  dose.column <- list(time=all.vars(form.pars.dose$rhs),
                      groups=all.vars(form.parts.dose$groups))
  if (length(conc.column$conc) != 1)
    stop("There must be exactly one concentration variable")
  if (length(conc.column$time) != 1)
    stop("There must be exactly one time variable for concentration")
  if (!identical(form.pars.dose$lhs, NA))
    stop("The dosing formula must be one-sided")
  if (length(dose.time.column) != 1)
    stop("There must be exactly one time variable for dosing")
  ## Check that all the formula variables are present in the
  ## associated data frame.
  for (n in names(conc.column))
    if (!(all(conc.column[[n]] %in% names(data.conc))))
      stop(sprintf("The %s column '%s' is not in data.conc",
                   n,
                   paste(setdiff(dose.column[[n]], names(data.dose)),
                         collapse=", ")))
  for (n in names(dose.column))
    if (!(all(dose.column[[n]] %in% names(data.dose))))
      stop(sprintf("The %s column '%s' is not in data.dose",
                   n,
                   paste(setdiff(dose.column[[n]], names(data.dose)),
                         collapse=", ")))
  ## Check that the groups are either the same or that dosing is
  ## shorter than concentration.
  if (length(conc.column$groups) == length(dose.column$groups)) {
    if (!all(dose.column$groups %in% conc.column.groups))
      stop("Some groups for the dose are not in the concentration groups")
  } else if (length(dose.column$groups) > length(conc.column$groups)) {
    
  }
}

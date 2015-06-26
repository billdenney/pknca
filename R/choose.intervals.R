#' Choose intervals to compute AUCs from time and dosing information
#'
#' Intervals for AUC are selected by the following metrics:
#' \enumerate{
#'   \item If only one dose is administered, use the
#'         \code{PKNCA.options("single.dose.auc")}
#'   \item If more than one dose is administered, estimate the AUC
#'         between any two doses that have PK taken at both of the
#'         dosing times and at least one time  between the doses.
#'   \item For the final dose of multiple doses, try to determine the
#'         dosing interval (\eqn{\tau}) and estimate the AUC in that
#'         interval if multiple samples are taken in the interval.
#'   \item If there are samples \eqn{> \tau} after the last dose,
#'         calculate the half life after the last dose.
#'  }
#' 
#' @param time.conc Time of concentration measurement
#' @param time.dosing Time of dosing
#' @param options List of changes to the default
#' \code{\link{PKNCA.options}} for calculations.
#' @param single.dose.aucs The AUC specification for single dosing.
#' @return A data frame with columns for \code{start}, \code{end},
#' \code{auc.type}, and \code{half.life}.  See
#' \code{\link{check.interval.specification}} for column definitions.
#' The data frame may have zero rows if no intervals could be found.
#' @seealso \code{\link{pk.calc.auc}}, \code{\link{pk.calc.aumc}},
#' \code{\link{pk.calc.half.life}}, \code{\link{find.tau}},
#' \code{\link{check.interval.specification}}, \code{\link{PKNCA.options}}
#' @export
choose.auc.intervals <- function(time.conc, time.dosing,
                                 options=list(),
                                 single.dose.aucs=PKNCA.choose.option("single.dose.aucs", options)) {
  if (any(is.na(time.conc)))
    stop("time.conc may not have any NA values")
  if (any(is.na(time.dosing)))
    stop("time.dosing may not have any NA values")
  if (length(unique(time.dosing)) == 1) {
    ## If it is single-dose data, use the time of dosing and then
    ## offset it by the dosing time (allowing the case where dosing
    ## time is not 0).
    ret <- check.interval.specification(single.dose.aucs)
    ## If there is an offset from 0, use that offset
    ret$start <- ret$start + unique(time.dosing)
    ret$end <- ret$end + unique(time.dosing)
  } else {
    ## Sort the times so sorting can be assumed farther down in the
    ## algorithm
    time.dosing <- sort(time.dosing)
    time.conc <- sort(time.conc)
    ## Find the doses that have concentration measurements
    mask.dose.conc <- time.dosing %in% time.conc
    ## Find indexes of pairs of doses that both have predose PK associated
    idx.paired.dose <-
      (1:(length(time.dosing)-1))[mask.dose.conc[-1] &
                                  mask.dose.conc[-length(mask.dose.conc)]]
    ## A data frame with all the right columns and classes but no data
    ret <- check.interval.specification(data.frame(start=0, end=1, auclast=TRUE))[-1,]
    ## Find the pairs that have at least one measurement between them
    for (n in idx.paired.dose) {
      if (any(time.dosing[n] < time.conc &
              time.conc < time.dosing[n+1])) {
        ## If there are measurements between the doses, add it to the
        ## output.
        ret <- rbind(
          ret,
          check.interval.specification(
            data.frame(start=time.dosing[n],
                       end=time.dosing[n+1],
                       auclast=TRUE,
                       cmax=TRUE,
                       tmax=TRUE)))
      }
    }
    ## Find the repeating dosing interval if possible and add it to
    ## the last dose if there is PK beyond the last dose to that time.
    tau <- find.tau(time.dosing)
    if (!is.na(tau)) {
      if ((max(time.dosing) + tau) %in% time.conc) {
        ret <- rbind(
          ret,
          check.interval.specification(
            data.frame(start=max(time.dosing),
                       end=max(time.dosing) + tau,
                       cmax=TRUE,
                       tmax=TRUE,
                       auclast=TRUE,
                       stringsAsFactors=FALSE)))
      }
      ## If the maximum concentration measurement time is beyond the
      ## max dosing time + tau, calculate a half-life.
      if ((max(time.dosing) + tau) < max(time.conc)) {
        ret <- rbind(
          ret,
          check.interval.specification(
            data.frame(start=max(time.dosing),
                       end=Inf,
                       half.life=TRUE)))
      }
    }
  }
  ret
}

#' Find the repeating interval within a vector of doses
#'
#' This is intended to find the interval over which x repeats by the
#' rule unique(mod(x, interval)) is minimized.
#'
#' @param x the vector to find the interval within
#' @param na.action What to do with NAs in \code{x}
#' @param tau.choices the intervals to look for if the doses are not
#' all equally spaced.
#' @param options List of changes to the default
#' \code{\link{PKNCA.options}} for calculations.
#' @return A scalar indicating the repeating interval with the most
#' repetition.
#' \enumerate{
#'   \item If all values are \code{NA} then NA is returned.
#'   \item If all values are the same, then 0 is returned.
#'   \item If all values are equally spaced, then that spacing is
#'         returned.
#'   \item If one of the \code{choices} can minimize the number of
#'         unique values, then that is returned.
#'   \item If none of the \code{choices} can minimize the number of
#'         unique values, then -1 is returned.
#' }
#' @export
find.tau <- function(x, na.action=na.omit,
                     options=list(),
                     tau.choices=PKNCA.choose.option("tau.choices", options)) {
  ret <- NA
  x <- na.action(x)
  if (length(unique(x)) == 1) {
    ## Single dose, no more effort needed
    ret <- 0
  } else if (identical(tau.choices, NA)) {
    all.deltas <-
      sort(unique(
        as.vector(sapply(x, FUN=function(x, y) x - y, y=x))))
    tau.choices <- all.deltas[all.deltas > 0]
  }
  if (is.na(ret) &
      length(x) > 1) {
    delta.1 <- x[2] - x[1]
    if (all((x[-1] - x[-length(x)]) == delta.1)) {
      ## Only one interval through the full data set
      ret <- delta.1
    } else {
      ## Drop any tau.choices that are >= the difference in the range
      ## of x because those are uninformative (i.e. if the maximum
      ## time is 12 hours, don't test an interval of 12, 24, or
      ## ... hours because they will match the x - tau < 0 test in a
      ## meaningless way.
      r.x <- range(x)
      tau.choices <- tau.choices[tau.choices < (r.x[2] - r.x[1])]
      ## Ensure that the choices are in order so that we find the
      ## minimum interval.
      tau.choices <- sort(tau.choices)
      ## Test all the tau.choices until we find the first (and thereby
      ## smallest) usable one
      i <- 0
      while (is.na(ret) &
             i < length(tau.choices)) {
        i <- i+1
        tau <- tau.choices[i]
        ## Is the dose either within the first tau or there is a dose
        ## that far before it?
        dose.before <- ((x - tau < 0) | ((x - tau) %in% x))
        ## And
        ## Is the dose either within the last tau or there is a dose
        ## that far after it? 
        dose.after <- ((x + tau > max(x)) | ((x + tau) %in% x))
        if (all(dose.before & dose.after))
          ret <- tau
      }
    }
  }
  ret
}

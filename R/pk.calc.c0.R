#' Estimate the concentration at dosing time for an IV bolus dose.
#'
#' @param conc The observed concentrations
#' @param time The observed times
#' @param time.dose The time when dosing occurred
#' @param method The order of methods to test (see details)
#' @param check Check the \code{conc} and \code{time} inputs
#' @return The estimated concentration at time 0.
#'
#' @details Methods available for interpolation are below, and each
#' has its own specific function.
#' 
#' \describe{
#'   \item{\code{c0}}{If the observed \code{conc} at \code{time.dose} is nonzero, return that.  This method should always be used first.}
#'   \item{\code{logslope}}{Compute the semilog line between the first two measured times, and use that line to extrapolate backward to \code{time.dose}}
#'   \item{\code{c1}}{Use the first point after \code{time.dose}}
#' }
#' @export
pk.calc.c0 <- function(conc, time, time.dose=0,
                       method=c("c0", "logslope", "c1"),
                       check=TRUE) {
  ## Check the inputs
  if (check)
    check.conc.time(conc, time)
  if (length(time.dose) != 1)
    stop("time.dose must be a scalar")
  if (!is.numeric(time.dose) | is.factor(time.dose))
    stop("time.dose must be a number")
  if (time.dose > max(time)) {
    warning("time.dose is after all available data")
    return(NA)
  }
  method <- match.arg(method, several.ok=TRUE)
  ## Find the value
  ret <- NA
  while (is.na(ret) &
         length(method) > 0) {
    current.method <- method[1]
    method <- method[-1]
    ret <- do.call(paste("pk.calc.c0.method", current.method, sep="."),
                   args=list(conc=conc, time=time,
                     time.dose=time.dose, check=FALSE))
  }
  ret
}

#' @describeIn pk.calc.c0 Semilog regress the first and second points
#' after time.dose.  This method will return \code{NA} if the second
#' \code{conc} after \code{time.dose} is 0 or greater than the first.
pk.calc.c0.method.logslope <- function(conc, time, time.dose=0,
                                       check=TRUE) {
  if (check)
    check.conc.time(conc, time)
  mask.positive.time <- (time > time.dose &
                         !(is.na(conc)))
  positive.time <- time[mask.positive.time]
  if (length(positive.time) < 2)
    return(NA)
  ## If there is enough data, proceed to calculate
  mask.1 <- time %in% positive.time[1]
  mask.2 <- time %in% positive.time[2]
  c1 <- conc[mask.1]
  c2 <- conc[mask.2]
  t1 <- time[mask.1]
  t2 <- time[mask.2]
  if (c2 < c1 &
      c2 != 0) {
    exp(log(c1) - (log(c2)-log(c1))/(t2-t1)*(t1 - time.dose))
  } else {
    NA
  }
}

#' @describeIn pk.calc.c0 Use \code{C0} = \code{conc[time %in%
#' time.dose]} if it is nonzero.
pk.calc.c0.method.c0 <- function(conc, time, time.dose=0, check=TRUE) {
  if (check)
    check.conc.time(conc, time)
  ## If there is a non-missing and nonzero concentration measurement
  ## at time.dose, that's our answer.
  mask.dose <- (time %in% time.dose &
                !(conc %in% c(NA, 0)))
  if (any(mask.dose)) {
    conc[mask.dose]
  } else {
    NA
  }
}

#' @describeIn pk.calc.c0 Use \code{C0} = \code{C1}.
pk.calc.c0.method.c1 <- function(conc, time, time.dose=0, check=TRUE) {
  if (check)
    check.conc.time(conc, time)
  mask.post.dose <- (time > time.dose &
                     !is.na(conc))
  if (any(mask.post.dose)) {
    conc[mask.post.dose][1]
  } else {
    NA
  }
}

#' @describeIn pk.calc.c0 Use \code{C0} = 0 (typically used for single
#' dose oral and IV infusion)
pk.calc.c0.method.set0 <- function(conc, time, time.dose=0, check=TRUE)
  0

#' @describeIn pk.calc.c0 Use \code{C0} = Cmin (typically used for
#' multiple dose oral and IV infusion)
pk.calc.c0.method.cmin <- function(conc, time, time.dose=0, check=TRUE)
  pk.calc.cmin(conc, check=check)

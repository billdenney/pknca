#' Compute noncompartmental superposition for repeated dosing
#'
#' @param conc Concentration measured
#' @param time Time of concentration measurement
#' @param tau The dosing interval
#' @param dose.times The time of dosing within the dosing interval.
#' The \code{min(dose.times)} must be >= 0, and the
#' \code{max(dose.times)} must be < \code{tau}.  There may be more
#' than one dose times given as a vector.
#' @param n.tau The number of tau dosing intervals to simulate or
#' \code{Inf} for steady-state.
#' @param lambda.z The elimination rate (from the half life
#' calculation, used to extrapolate beyond the maximum time observed).
#' @param clast.pred To use predicted as opposed to observed Clast,
#' either give the value for clast.pred here or set it to true (for
#' automatic calculation from the half-life).
#' @param additional.times Times to include in the final outputs in
#' addition to the standard times (see details).  All
#' \code{min(additional.times)} must be >= 0, and the
#' \code{max(additional.times)} must be <= \code{tau}.
#' @param steady.state.tol The tolerance for assessing if steady-state
#' has been achieved (between 0 and 1, exclusive).
#' @param interp.method See \code{\link{interp.extrap.conc}}
#' @param extrap.method See \code{\link{interp.extrap.conc}}
#' @param ... Additional arguments passed to the \code{half.life}
#' function if required to compute \code{lambda.z}.
#' @return A data frame with columns named "conc" and "time".
#'
#' @details
#' The returned superposition times will include all of the
#' following times: 0 (zero), \code{dose.times}, \code{time modulo tau}
#' (shifting \code{time} for each dose time as well),
#' \code{additional.times}, and \code{tau}.
#'
#' @seealso \code{\link{interp.extrap.conc}}
#' 
#' @export
superposition <- function(conc, time, tau, dose.times=0, n.tau=Inf,
                          options=list(),
                          lambda.z, clast.pred=FALSE, tlast,
                          additional.times=c(),
                          interp.method=PKNCA.choose.option("auc.method",
                            options),
                          extrap.method="AUCinf",
                          steady.state.tol=1e-3, ...) {
  ## Check the inputs
  check.conc.time(conc, time)
  if (!is.numeric(tau) | is.factor(tau))
    stop("tau must be a number")
  if (length(tau) != 1)
    stop("tau must be a scalar")
  if (tau <= 0)
    stop("tau must be > 0")
  if (!is.numeric(dose.times) | is.factor(dose.times))
    stop("dose.times must be a number")
  if (length(dose.times) < 1)
    stop("There must be at least one dose time")
  if (any(dose.times < 0))
    stop("All dose.times must be non-negative")
  if (any(dose.times >= tau))
    stop("dose.times must be < tau")
  if (!is.numeric(n.tau) | is.factor(n.tau))
    stop("n.tau must be a number")
  if (length(n.tau) != 1)
    stop("n.tau must be a scalar")
  if (n.tau < 1)
    stop("n.tau must be >= 1")
  ## FIXME: this should be something like 10*machine precision
  if (!is.infinite(n.tau) &
      ((n.tau %% 1) > 1e-7))
    stop("n.tau must be an integer or Inf")
  if (!missing(lambda.z)) {
    if (!is.numeric(lambda.z) | is.na(lambda.z) | is.factor(lambda.z))
      stop("lambda.z must be a number or NA")
    if (length(lambda.z) != 1)
      stop("lambda.z must be a scalar")
  }
  if (!is.logical(clast.pred) &
      (!is.numeric(clast.pred) |
       is.factor(clast.pred)))
    stop("clast.pred must either be a logical (TRUE/FALSE) or numeric value")
  if (length(clast.pred) != 1)
    stop("clast.pred must be a scalar")
  if (is.na(clast.pred))
    clast.pred <- FALSE
  if (is.numeric(clast.pred) & clast.pred <= 0)
    stop("clast.pred must be positive (if it is a number)")
  if (!missing(tlast)) {
    if (!is.numeric(tlast) | is.factor(tlast))
      stop("tlast must be a number")
    if (length(tlast) != 1)
      stop("tlast must be a scalar")
  }
  if (missing(lambda.z) & is.numeric(clast.pred))
    warning("clast.pred was given, but lambda.z was not given.  Usually clast.pred would only be given with a lambda.z value.")
  if (length(additional.times) > 0) {
    if (any(is.na(additional.times)))
      stop("No additional.times may be NA (to not include any additional.times, enter c() as the function argument)")
    if (!is.numeric(additional.times) | is.factor(additional.times))
      stop("additional.times must be a number")
    if (any(additional.times < 0))
      stop("All additional.times must be nonnegative")
    if (any(additional.times > tau))
      stop("All additional.times must be <= tau")
  }
  if (!is.numeric(steady.state.tol) | is.factor(steady.state.tol))
    stop("steady.state.tol must be a number")
  if (length(steady.state.tol) != 1)
    stop("steady.state.tol must be a scalar")
  if (steady.state.tol <= 0 |
      steady.state.tol >= 1)
    stop("steady.state.tol must be between 0 and 1, exclusive.")
  if (steady.state.tol > 0.01)
    warning("steady.state.tol is usually <= 0.01")
  ## We get all or none of lambda.z, clast, and tlast
  has.lambda.z <- !missing(lambda.z)
  has.clast.pred <- !is.logical(clast.pred)
  has.tlast <- !missing(tlast)
  if (any(c(has.lambda.z, has.clast.pred, has.tlast)) &
      !all(c(has.lambda.z, has.clast.pred, has.tlast)))
    stop("Either give all or none of the values for these arguments: lambda.z, clast.pred, and tlast")
  ## Determine the output times
  tmp.times <- c(0, tau, dose.times, additional.times)
  ## For the output, ensure that we have each of the input times
  ## shifted for the dosing times.
  for (d in dose.times)
    tmp.times <- c(tmp.times, (d + times) %% tau)
  ret <- data.frame(conc=NA, time=sort(unique(tmp.times)))
  ## Check if all the input concentrations are 0, and if so, give that
  ## simple answer.
  if (all(input.data$conc == 0)) {
    ret$conc <- 0
  } else {
    ## Check if we will need to extrapolate
    if (tau*n.tau > tlast) {
      if (missing(lambda.z)) {
        tmp <- pk.calc.half.life(conc, time, options=options,
                                 ..., check=FALSE)
        tlast <- tmp$tlast
        lambda.z <- tmp$lambda.z
      }
      if (identical(clast.pred, FALSE)) {
        clast <- pk.calc.clast.obs(conc, time, check=FALSE)
      } else if (identical(clast.pred, TRUE)) {
        clast <- tmp$clast.pred
      }
    }
  }
  
  time.now <- as.numeric(times)
  ## This mask allows superposition to continue if a time is NA
  ## (i.e. if the time is missing)
  mask <- !is.na(time.now)
  if (any(time.now[mask] > tau) & is.infinite(n.doses)) {
    stop(sprintf("Cannot interpolate times beyond steady state tau (%g): %s",
                 tau, paste(time.now, collapse=", ")))
  }
  if (missing(lz))
    lz = lambda.z(data, ...)

  current <- data.frame(time=time.now, conc=0)
  current$conc[mask] <- interpolate.conc(data, time.now[mask], lz=lz, ...)$conc
  if (any(mask)) {
    ## If there is at least one time to extrapolate
    delta <- superpos.tol+1
    dose.num <- 1
    while (!(all(is.na(delta))) &
           any(na.omit(delta) > superpos.tol) &
           dose.num < n.doses) {
      previous <- current
      time.now[mask] <- time.now[mask] + tau
      current$conc[mask] <- current$conc[mask] +
        interpolate.conc(data, time.now[mask], lz=lz, ...)$conc
      delta <- current$conc[mask]/previous$conc[mask] - 1
    }
  }
  return(current)
}

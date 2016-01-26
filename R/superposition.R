#' Compute noncompartmental superposition for repeated dosing
#'
#' @param conc Concentration measured
#' @param time Time of concentration measurement
#' @param dose.input The dose given to generate the \code{conc} and
#' \code{time} inputs.  If missing, output doses will be assumed to be
#' equal to the input dose.
#' @param tau The dosing interval
#' @param dose.times The time of dosing within the dosing interval.
#' The \code{min(dose.times)} must be >= 0, and the
#' \code{max(dose.times)} must be < \code{tau}.  There may be more
#' than one dose times given as a vector.
#' @param dose.amount The doses given for the output.  Linear
#' proportionality will be used from the input to output if they are
#' not equal.  The length of dose.amount must be either 1 or matching
#' the length of \code{dose.times}.
#' @param n.tau The number of tau dosing intervals to simulate or
#' \code{Inf} for steady-state.
#' @param options The PKNCA.options to use for the calculation (passed
#' on to subsequent functions like \code{pk.calc.half.life}).
#' @param lambda.z The elimination rate (from the half life
#' calculation, used to extrapolate beyond the maximum time observed).
#' @param clast.pred To use predicted as opposed to observed Clast,
#' either give the value for clast.pred here or set it to true (for
#' automatic calculation from the half-life).
#' @param tlast The time of last observed concentration above the
#' limit of quantificaiton.  This is calculated if not provided.
#' @param additional.times Times to include in the final outputs in
#' addition to the standard times (see details).  All
#' \code{min(additional.times)} must be >= 0, and the
#' \code{max(additional.times)} must be <= \code{tau}.
#' @param check.blq Must the first concentration measurement be below
#' the limit of quantification?
#' @param interp.method See \code{\link{interp.extrap.conc}}
#' @param extrap.method See \code{\link{interp.extrap.conc}}
#' @param steady.state.tol The tolerance for assessing if steady-state
#' has been achieved (between 0 and 1, exclusive).
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
superposition <- function(conc, ...)
  UseMethod("superposition", conc)

#' @rdname superposition
#' @export
superposition.PKNCAconc <- function(conc, ...) {
  conc.col <- as.character(parseFormula(conc)$lhs)
  time.col <- as.character(parseFormula(conc)$rhs)
  ## Split the data by grouping and extract just the concentration and
  ## time columns
  tmp.data <- doBy::splitBy(parseFormula(conc)$groupFormula,
                            getData(conc))
  groupinfo <- attributes(tmp.data)$groupid
  tmp.data <-
    parallel::mclapply(X=tmp.data,
                       FUN=function(x, conc.col, time.col) {
                         renameCol(x, c(conc.col, time.col),
                                   c("conc", "time"))[,c("conc", "time")]
                       },
                       conc.col=conc.col,
                       time.col=time.col)
  tmp.results <-
    parallel::mclapply(X=tmp.data,
                       FUN=function(x, ...) {
                         superposition.numeric(x$conc, x$time, ...)
                       }, ...)
  ret <- data.frame()
  for (i in seq_along(tmp.results))
    ret <- rbind(ret,
                 cbind(groupinfo[i,,drop=FALSE],
                       tmp.results[[i]],
                       row.names=NULL))
  ret
}

#' @rdname superposition
#' @export
superposition.numeric <- function(conc, time, dose.input,
                                  tau, dose.times=0, dose.amount, n.tau=Inf,
                                  options=list(),
                                  lambda.z, clast.pred=FALSE, tlast,
                                  additional.times=c(),
                                  check.blq=TRUE,
                                  interp.method=PKNCA.choose.option("auc.method",
                                    options),
                                  extrap.method="AUCinf",
                                  steady.state.tol=1e-3, ...) {
  ## Check the inputs
  ## Concentration and time
  check.conc.time(conc, time)
  if (check.blq)
    if (!(conc[1] %in% 0))
      stop("The first concentration must be 0 (and not NA).  To change this set check.blq=FALSE.")
  ## dose.input
  if (!missing(dose.input)) {
    if (length(dose.input) != 1)
      stop("dose.input must be a scalar")
    if (!is.numeric(dose.input) | is.factor(dose.input) | is.na(dose.input))
      stop("dose.input must be a number")
    if (dose.input <= 0)
      stop("dose.input must be > 0")
  }
  ## tau
  if (length(tau) != 1)
    stop("tau must be a scalar")
  if (!is.numeric(tau) | is.factor(tau) | is.na(tau))
    stop("tau must be a number")
  if (tau <= 0)
    stop("tau must be > 0")
  ## dose.times
  if (length(dose.times) < 1)
    stop("There must be at least one dose time")
  if (!is.numeric(dose.times) | is.factor(dose.times) | any(is.na(dose.times)))
    stop("dose.times must be a number")
  if (any(dose.times < 0))
    stop("All dose.times must be non-negative")
  if (any(dose.times >= tau))
    stop("dose.times must be < tau")
  ## dose.amount
  if (!missing(dose.amount)) {
    if (missing(dose.input))
      stop("must give dose.input to give dose.amount")
    if (!(length(dose.amount) %in% c(1, length(dose.times))))
      stop("dose.amount must either be a scalar or match the length of dose.times")
    if (!is.numeric(dose.amount) | is.factor(dose.amount))
      stop("dose.amount must be a number")
    if (any(is.infinite(dose.amount)) | any(is.na(dose.amount)))
      stop("dose.amount must be finite and not NA")
    if (any(dose.amount <= 0))
      stop("All dose.amount must be positive")
  }
  ## n.tau
  if (length(n.tau) != 1)
    stop("n.tau must be a scalar")
  if (!is.numeric(n.tau) | is.factor(n.tau) | is.na(n.tau))
    stop("n.tau must be a number")
  if (n.tau < 1)
    stop("n.tau must be >= 1")
  if (!is.infinite(n.tau) &
      ((n.tau %% 1) > 10*.Machine$double.eps))
    stop("n.tau must be an integer or Inf")
  ## lambda.z
  if (!missing(lambda.z)) {
    if (length(lambda.z) != 1)
      stop("lambda.z must be a scalar")
    if (!is.numeric(lambda.z) | is.factor(lambda.z))
      stop("lambda.z must be a number")
  }
  ## clast.pred
  if (length(clast.pred) != 1)
    stop("clast.pred must be a scalar")
  if (is.na(clast.pred))
    clast.pred <- FALSE
  if (!is.logical(clast.pred) &
      (!is.numeric(clast.pred) |
       is.factor(clast.pred)))
    stop("clast.pred must either be a logical (TRUE/FALSE) or numeric value")
  if (is.numeric(clast.pred) & clast.pred <= 0)
    stop("clast.pred must be positive (if it is a number)")
  ## tlast
  if (!missing(tlast)) {
    if (length(tlast) != 1)
      stop("tlast must be a scalar")
    if (!is.numeric(tlast) | is.factor(tlast) | is.na(tlast))
      stop("tlast must be a number")
  }
  ## additional.times
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
  ## steady.state.tol
  if (length(steady.state.tol) != 1)
    stop("steady.state.tol must be a scalar")
  if (!is.numeric(steady.state.tol) | is.factor(steady.state.tol) | is.na(steady.state.tol))
    stop("steady.state.tol must be a number")
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
  ## combine dose.input and dose.amount as applicable to scale the
  ## outputs.
  if (!missing(dose.amount)) {
    dose.scaling <- dose.amount / dose.input
    if (length(dose.scaling) != length(dose.times)) {
      if (length(dose.scaling) != 1)
        stop("bug in dose.amount, dose.times, and dose.input handling")
      ## it is a scalar and there is more than one dose
      dose.scaling <- rep(dose.scaling, length(dose.times))
    }
  } else {
    dose.scaling <- rep(1, length(dose.times))
  }
  ## Determine the output times
  tmp.times <- c(0, tau, dose.times, additional.times)
  ## For the output, ensure that we have each of the input times
  ## shifted for the dosing times.
  for (d in dose.times)
    tmp.times <- c(tmp.times, (d + time) %% tau)
  ret <- data.frame(conc=0, time=sort(unique(tmp.times)))
  ## Check if all the input concentrations are 0, and if so, give that
  ## simple answer.
  if (all(conc == 0)) {
    return(ret)
  } else {
    ## Check if we will need to extrapolate
    tlast <- pk.calc.tlast(conc, time, check=FALSE)
    if ((tau*n.tau) > tlast) {
      if (missing(lambda.z)) {
        tmp <- pk.calc.half.life(conc, time, options=options,
                                 ..., check=FALSE)
        lambda.z <- tmp$lambda.z
      }
      if (identical(clast.pred, FALSE)) {
        clast <- pk.calc.clast.obs(conc, time, check=FALSE)
      } else if (identical(clast.pred, TRUE)) {
        clast <- tmp$clast.pred
      } else {
        ## Use the provided clast.pred
        clast <- clast.pred
      }
    } else {
      ## We don't need lambda.z for calculations, but it is simpler to
      ## define it.
      lambda.z <- NA
    }
  }
  ## cannot continue extrapolating due to missing data (likely due to
  ## half-life not calculable)
  if ((n.tau * tau) > tlast & is.na(lambda.z)) {
    ret$conc <- NA
  } else {
    ## Do the math! (Finally)
    current.tol <- steady.state.tol + 1
    tau.count <- 0
    ## Stop either for reaching steady-state or for reaching the requested number of doses
    while (tau.count < n.tau &
           !is.na(current.tol) & 
           current.tol >= steady.state.tol) {
      prev.conc <- ret$conc
      ## Perform the dosing for a single dosing interval.
      for (i in seq_along(dose.times)) {
        tmp.time <- ret$time - dose.times[i] + (tau * tau.count)
        ## For the first dosing interval, make sure that we do not
        ## assign concentrations before the dose is given.
        mask.time <- tmp.time >= 0
        ## Update the current concentration (previous concentration +
        ## new concentration scaled by the relative dose)
        ret$conc[mask.time] <-
          (ret$conc[mask.time] +
           dose.scaling[i]*
           interp.extrap.conc(conc, time, time.out=tmp.time[mask.time],
                              lambda.z=lambda.z, clast=clast,
                              options=options, check=FALSE))
      }
      tau.count <- tau.count + 1
      if (any(ret$conc %in% 0)) {
        ## prevent division by 0.  Since not all concentrations are 0,
        ## all values will eventually be nonzero.
        current.tol <- steady.state.tol + 1
      } else {
        current.tol <- max(1-(prev.conc/ret$conc))
      }
    }
  }
  ret
}

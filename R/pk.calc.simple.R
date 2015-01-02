## Calculate the simple parameters for PK.

#' Calculate the adjusted r-squared value
#'
#' @param r.sq The r-squared value
#' @param n The number of points
#' @return The numeric adjusted r-squared value
#' @export
adj.r.squared <- function(r.sq, n) {
  if (n <= 2)
    stop("n must be > 2")
  1-(1-r.sq)*(n-1)/(n-2)
}

#' Determine maximum observed PK concentration
#'
#' @param conc Concentration measured
#' @param check Run \code{\link{check.conc.time}}?
#' @return a number for the maximum concentration or NA if all
#' concentrations are missing
#' @export
pk.calc.cmax <- function(conc, check=TRUE) {
  if (check)
    check.conc.time(conc=conc)
  if (length(conc) == 0 | all(is.na(conc))) {
    NA
  } else {
    max(conc, na.rm=TRUE)
  }
}

#' @describeIn pk.calc.cmax Determine the minimum observed PK
#' concentration
pk.calc.cmin <- function(conc, check=TRUE) {
  if (check)
    check.conc.time(conc=conc)
  if (length(conc) == 0 | all(is.na(conc))) {
    NA
  } else {
    min(conc, na.rm=TRUE)
  }
}

#' Determine time of maximum observed PK concentration
#'
#' Input restrictions are:
#' \enumerate{
#'   \item the \code{conc} and \code{time} must be the same length,
#'   \item the \code{time} may have no NAs,
#' }
#' \code{NA} will be returned if:
#' \enumerate{
#'   \item the length of \code{conc} and \code{time} is 0
#'   \item all \code{conc} is 0 or \code{NA}
#' }
#' 
#' @param conc Concentration measured
#' @param time Time of concentration measurement
#' @param use.first If there is more than time that matches the
#' maximum concentration, should the first be considered as Tmax?  If
#' not, then the last is considered Tmax.  (Default: TRUE)
#' @param check Run \code{\link{check.conc.time}}?
#' @return the time of the maximum concentration
#' @export
pk.calc.tmax <- function(conc, time,
                         use.first=PKNCA.options("first.tmax"),
                         check=TRUE) {
  if (missing(conc))
    stop("conc must be given")
  if (missing(time))
    stop("time must be given")
  if (check)
    check.conc.time(conc, time)
  if (length(conc) == 0 | all(conc %in% c(NA, 0))) {
    NA
  } else {
    ret <- time[conc %in% pk.calc.cmax(conc, check=FALSE)]
    if (use.first) {
      ret[1]
    } else {
      ret[length(ret)]
    }
  }
}

#' Determine time of last observed concentration above the limit of
#' quantification.
#'
#' \code{NA} will be returned if all \code{conc} are \code{NA} or 0.
#'
#' @param conc Concentration measured
#' @param time Time of concentration measurement
#' @param check Run \code{\link{check.conc.time}}?
#' @return The time of the last observed concentration measurement
#' @export
pk.calc.tlast <- function(conc, time, check=TRUE) {
  if (missing(conc))
    stop("conc must be given")
  if (missing(time))
    stop("time must be given")
  if (check)
    check.conc.time(conc, time)
  if (all(conc %in% c(NA, 0))) {
    NA
  } else {
    max(time[!(conc %in% c(NA, 0))])
  }
}

#' @describeIn pk.calc.tlast Determine the first concentration above
#' the limit of quantification.
pk.calc.tfirst <- function(conc, time, check=TRUE) {
  if (missing(conc))
    stop("conc must be given")
  if (missing(time))
    stop("time must be given")
  if (check)
    check.conc.time(conc, time)
  if (all(conc %in% c(NA, 0))) {
    NA
  } else {
    min(time[!(conc %in% c(NA, 0))])
  }
}

#' Determine the last observed concentration above the limit of
#' quantification (LOQ).
#'
#' If Tlast is NA (due to no non-missing above LOQ measurements), this
#' will return NA.
#' @param conc Concentration measured
#' @param time Time of concentration measurement
#' @param check Run \code{\link{check.conc.time}}?
#' @return The last observed concentration above the LOQ
#' @export
pk.calc.clast.obs <- function(conc, time, check=TRUE) {
  if (check)
    check.conc.time(conc, time)
  tlast <- pk.calc.tlast(conc, time)
  if (!is.na(tlast)) {
    conc[time %in% tlast]
  } else {
    NA
  }
}

#' Calculate the effective half-life
#'
#' @param mrt the mean residence time
#' @return the numeric value of the effective half-life
#' @export
pk.calc.thalf.eff <- function(mrt)
  log(2)*mrt

#' Calculate the AUC percent extrapolated
#'
#' @param auclast the area under the curve from time 0 to the last
#' measurement above the limit of quantification
#' @param aucinf the area under the curve from time 0 to infinity
#' @return the numeric value of the AUC percent extrapolated
#' @export
pk.calc.aucpext <- function(auclast, aucinf)
  100*(1-auclast/aucinf)

#' Calculate the elimination rate (Kel)
#'
#' @param mrt the mean residence time
#' @return the numeric value of the elimination rate
#' @export
pk.calc.kel <- function(mrt)
  1/mrt

#' Calculate the (observed oral) clearance
#'
#' @param dose the dose administered
#' @param auc the area under the curve from 0 to infinity or 0 to tau
#' (the next dose on a regular schedule at steady-state)
#' @param unitconv the multiplied factor to use for unit conversion
#' (e.g. 1000 for mg \code{dose}, time*ng/mL for \code{auc}, and
#' output in L/time)
#' @return the numeric value of the total (CL) or observed oral
#' clearance (CL/F)
#' @references Gabrielsson J, Weiner D.  "Section 2.5.1 Derivation of
#' clearance."  Pharmacokinetic & Pharmacodynamic Data Analysis:
#' Concepts and Applications, 4th Edition.  Stockholm, Sweden: Swedish
#' Pharmaceutical Press, 2000.  86-7.
#' @export
pk.calc.cl <- function(dose, auc, unitconv)
  if (missing(unitconv)) {
    dose/auc
  } else {
    unitconv*dose/auc
  }

#' Calculate the absolute (or relative) bioavailability
#'
#' @param dose1 The dose administered in route or method 1
#' @param dose2 The dose administered in route or method 2
#' @param auc1 The AUC from 0 to infinity or 0 to tau administered in
#' route or method 1
#' @param auc2 The AUC from 0 to infinity or 0 to tau administered in
#' route or method 2
#' @export
pk.calc.f <- function(dose1, auc1, dose2, auc2)
  (auc2/dose2)/(auc1/dose1)

#' Calcuate the mean residence time (MRT)
#'
#' @param auc the AUC from 0 to infinity or 0 ot tau at steady-state
#' @param aumc the AUMC from 0 to infinity or 0 ot tau at steady-state
#' @return the numeric value of the mean residence time
#' @export
pk.calc.mrt <- function(auc, aumc)
  aumc/auc

#' Calculate the terminal volume of distribution (Vz)
#'
#' @param dose the dose administered
#' @param auc the AUC from 0 to infinity
#' @param kel the elimination rate
#' @param unitconv the factor to use for unit conversion (e.g. 1000
#' for mg \code{dose}, time*ng/mL for \code{auc}, 1/time for
#' \code{kel}, and output in L)
#' @export
pk.calc.vz <- function(dose, auc, kel, unitconv) {
  ## Ensure that dose is either a scalar or the same length as AUC
  ## (more complex repeating patterns while valid for general R are
  ## likely errors here).
  if ((length(dose) != 1) & (length(dose) != length(auc)))
      stop("'dose' and 'auc' must be the same length")
  ## While dose may be the same for many measurements, it is highly
  ## unlikely that kel and auc are the same.  Ensure that kel and auc
  ## are the same length.
  if (length(auc) != length(kel))
    stop("'auc' and 'kel' must be the same length")
  if (missing(unitconv)) {
    dose/(auc*kel)
  } else {
    dose*unitconv/(auc*kel)
  }
}

#' Calculate the steady-state volume of distribution (Vss)
#'
#' @param cl the clearance
#' @param mrt the mean residence time
#' @return the volume of distribution at steady-state
#' @export
pk.calc.vss <- function(cl, mrt)
  cl*mrt

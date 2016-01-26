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
## Add the column to the interval specification
add.interval.col("cmax",
                 FUN="pk.calc.cmax",
                 values=c(FALSE, TRUE),
                 desc="Maximum observed concentration",
                 depends=c())
PKNCA.set.summary("cmax", business.geomean, business.geocv)

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
## Add the column to the interval specification
add.interval.col("cmin",
                 FUN="pk.calc.cmin",
                 values=c(FALSE, TRUE),
                 desc="Minimum observed concentration",
                 depends=c())
PKNCA.set.summary("cmin", business.geomean, business.geocv)

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
#' @param options List of changes to the default
#' \code{\link{PKNCA.options}} for calculations.
#' @param first.tmax If there is more than time that matches the
#' maximum concentration, should the first be considered as Tmax?  If
#' not, then the last is considered Tmax.
#' @param check Run \code{\link{check.conc.time}}?
#' @return the time of the maximum concentration
#' @export
pk.calc.tmax <- function(conc, time,
                         options=list(),
                         first.tmax=PKNCA.choose.option("first.tmax", options),
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
    if (first.tmax) {
      ret[1]
    } else {
      ret[length(ret)]
    }
  }
}
## Add the column to the interval specification
add.interval.col("tmax",
                 FUN="pk.calc.tmax",
                 values=c(FALSE, TRUE),
                 desc="Time of the maximum observed concentration",
                 depends=c())
PKNCA.set.summary("tmax", business.median, business.range)

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
## Add the column to the interval specification
add.interval.col("tlast",
                 FUN="pk.calc.tlast",
                 values=c(FALSE, TRUE),
                 desc="Time of the last concentration observed above the limit of quantification",
                 depends=c())
PKNCA.set.summary("tlast", business.median, business.range)

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
## Add the column to the interval specification
add.interval.col("tfirst",
                 FUN="pk.calc.tfirst",
                 values=c(FALSE, TRUE),
                 desc="Time of the first concentration above the limit of quantification",
                 depends=c())
PKNCA.set.summary("tfirst", business.median, business.range)

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
## Add the column to the interval specification
add.interval.col("clast.obs",
                 FUN="pk.calc.clast.obs",
                 values=c(FALSE, TRUE),
                 desc="The last concentration observed above the limit of quantification",
                 depends=c())
PKNCA.set.summary("clast.obs", business.geomean, business.geocv)

#' Calculate the effective half-life
#'
#' @param mrt the mean residence time
#' @return the numeric value of the effective half-life
#' @export
pk.calc.thalf.eff <- function(mrt)
  log(2)*mrt
## Add the column to the interval specification
add.interval.col("thalf.eff",
                 FUN="pk.calc.thalf.eff",
                 values=c(FALSE, TRUE),
                 desc="The effective half-life (as determined from the MRT)",
                 depends=c("mrt"))
PKNCA.set.summary("thalf.eff", business.geomean, business.geocv)

#' Calculate the AUC percent extrapolated
#'
#' @param auclast the area under the curve from time 0 to the last
#' measurement above the limit of quantification
#' @param aucinf the area under the curve from time 0 to infinity
#' @return the numeric value of the AUC percent extrapolated
#' @export
pk.calc.aucpext <- function(auclast, aucinf) {
  if (auclast >= aucinf)
    warning("auclast should be less than aucinf")
  100*(1-auclast/aucinf)
}
## Add the column to the interval specification
add.interval.col("aucpext",
                 FUN="pk.calc.aucpext",
                 values=c(FALSE, TRUE),
                 desc="Percent of the AUCinf that is extrapolated after Tlast",
                 depends=c("auclast", "aucinf"))
PKNCA.set.summary("aucpext", business.mean, business.sd)

#' Calculate the elimination rate (Kel)
#'
#' @param mrt the mean residence time
#' @return the numeric value of the elimination rate
#' @export
pk.calc.kel <- function(mrt)
  1/mrt
## Add the column to the interval specification
add.interval.col("kel",
                 FUN="pk.calc.kel",
                 values=c(FALSE, TRUE),
                 desc="Elimination rate (as calculated from the MRT)",
                 depends=c("mrt"))
PKNCA.set.summary("kel", business.geomean, business.geocv)

#' Calculate the (observed oral) clearance
#'
#' @param dose the dose administered
#' @param aucinf the area under the curve from 0 to infinity or 0 to
#'   tau (the next dose on a regular schedule at steady-state)
#' @param unitconv the multiplied factor to use for unit conversion
#'   (e.g. 1000 for mg \code{dose}, time*ng/mL for \code{auc}, and
#'   output in L/time)
#' @return the numeric value of the total (CL) or observed oral
#'   clearance (CL/F)
#' @references Gabrielsson J, Weiner D.
#'   "Section 2.5.1 Derivation of clearance."  Pharmacokinetic &
#'   Pharmacodynamic Data Analysis: Concepts and Applications, 4th
#'   Edition.  Stockholm, Sweden: Swedish Pharmaceutical Press, 2000.
#'   86-7.
#' @export
pk.calc.cl <- function(dose, aucinf, unitconv=NA)
  if (is.na(unitconv)) {
    dose/aucinf
  } else {
    unitconv*dose/aucinf
  }
## Add the column to the interval specification
add.interval.col("cl",
                 FUN="pk.calc.cl",
                 values=c(FALSE, TRUE),
                 desc="Clearance or observed oral clearance",
                 depends=list("aucinf", "auclast"))
PKNCA.set.summary("cl", business.geomean, business.geocv)

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
add.interval.col("f",
                 FUN="pk.calc.f",
                 values=c(FALSE, TRUE),
                 desc="Bioavailability or relative bioavailability",
                 depends=c())
PKNCA.set.summary("f", business.geomean, business.geocv)

#' Calcuate the mean residence time (MRT)
#'
#' @param auc the AUC from 0 to infinity or 0 ot tau at steady-state
#' @param aumc the AUMC from 0 to infinity or 0 ot tau at steady-state
#' @return the numeric value of the mean residence time
#' @export
pk.calc.mrt <- function(auc, aumc)
  aumc/auc
## Add the column to the interval specification
add.interval.col("mrt",
                 FUN="pk.calc.mrt",
                 values=c(FALSE, TRUE),
                 desc="The mean residence time",
                 depends=list(c("auclast", "aumclast"),
                              c("aucinf", "aumcinf")))
PKNCA.set.summary("mrt", business.geomean, business.geocv)

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
## Add the column to the interval specification
add.interval.col("vz",
                 FUN="pk.calc.vz",
                 values=c(FALSE, TRUE),
                 desc="The terminal volume of distribution",
                 depends=list(c("aucinf", "kel"),
                              c("auclast", "kel")))
PKNCA.set.summary("vz", business.geomean, business.geocv)

#' Calculate the steady-state volume of distribution (Vss)
#'
#' @param cl the clearance
#' @param mrt the mean residence time
#' @return the volume of distribution at steady-state
#' @export
pk.calc.vss <- function(cl, mrt)
  cl*mrt
## Add the column to the interval specification
add.interval.col("vss",
                 FUN="pk.calc.vss",
                 values=c(FALSE, TRUE),
                 desc="The steady-state volume of distribution",
                 depends=c("cl", "mrt"))
PKNCA.set.summary("vss", business.geomean, business.geocv)

#' Calculate the volume of distribution (Vd) or observed volume of
#' distribution (Vd/F)
#'
#' @param dose Dose given
#' @param aucinf Area under the curve to infinity (either predicted or
#' observed).
#' @param lambda.z Elimination rate constant
#' @return The observed volume of distribution
#' @export
pk.calc.vd <- function(dose, aucinf, lambda.z)
  dose/(aucinf * lambda.z)
add.interval.col("vd",
                 FUN="pk.calc.vd",
                 values=c(FALSE, TRUE),
                 desc="Apparent observed volume of distribution",
                 depends=c("aucinf", "lambda.z"))
PKNCA.set.summary("vd", business.geomean, business.geocv)

#' Calculate the average concentration during an interval.
#'
#' @param auclast The area under the curve during the interval
#' @param start The starting time of the interval
#' @param end The ending time of the interval
#' @return The Cav (average concentration during the interval)
#' @export
pk.calc.cav <- function(auclast, start, end)
  auclast/(end-start)
add.interval.col("cav",
                 FUN="pk.calc.cav",
                 values=c(FALSE, TRUE),
                 desc="The average concentration during an interval",
                 depends="auclast")
PKNCA.set.summary("cav", business.geomean, business.geocv)

#' Determine the trough (predose) concentration
#'
#' @param conc Observed concentrations during the interval
#' @param time Times of \code{conc} observations
#' @param start Starting time of the interval
#' @return The concentration when \code{time == start}.  If none
#' match, then \code{NA}
#' @export
pk.calc.ctrough <- function(conc, time, start) {
  check.conc.time(conc, time)
  mask.start <- time %in% start
  if (sum(mask.start) == 1) {
    conc[mask.start]
  } else if (sum(mask.start) == 0) {
    NA
  } else {
    stop("More than one time matches the starting time")
  }
}
add.interval.col("ctrough",
                 FUN="pk.calc.ctrough",
                 values=c(FALSE, TRUE),
                 desc="The trough (predose) concentration",
                 depends=c())
PKNCA.set.summary("ctrough", business.geomean, business.geocv)

#' Determine the peak-to-trough ratio
#'
#' @param cmax The maximum observed concentration
#' @param cmin The minimum observed concentration
#' @return The ratio of cmax to cmin (if cmin == 0, NA)
#' @export
pk.calc.ptr <- function(cmax, cmin) {
  ret <- cmax/cmin
  ret[cmin %in% 0] <- NA
  ret
}
add.interval.col("ptr",
                 FUN="pk.calc.ptr",
                 values=c(FALSE, TRUE),
                 desc="Peak-to-Trough ratio (fraction)",
                 depends=c("cmax", "cmin"))
PKNCA.set.summary("ptr", business.geomean, business.geocv)

#' Determine the observed lag time (time before the first
#' concentration above the limit of quantification or above the first
#' concentration in the interval)
#'
#' @param conc The observed concentrations
#' @param time The observed times
#' @return The time associated with the first increasing concentration
#' @export
pk.calc.tlag <- function(conc, time) {
  check.conc.time(conc, time)
  mask.increase <- conc[-1] > conc[-length(conc)]
  if (any(mask.increase)) {
    time[c(mask.increase, FALSE)][1]
  } else {
    NA
  }
}
add.interval.col("tlag",
                 FUN="pk.calc.tlag",
                 values=c(FALSE, TRUE),
                 desc="Lag time",
                 depends=c())
PKNCA.set.summary("tlag", business.median, business.range)

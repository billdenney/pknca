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
#'   concentration
#' @export
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
#'   the limit of quantification.
#' @export
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
#' @param mrt the mean residence time to infinity
#' @return the numeric value of the effective half-life
#' @export
pk.calc.thalf.eff <- function(mrt)
  log(2)*mrt
pk.calc.thalf.eff.obs <- function(mrt.obs)
  pk.calc.thalf.eff(mrt.obs)
pk.calc.thalf.eff.pred <- function(mrt.pred)
  pk.calc.thalf.eff(mrt.pred)
## Add the columns to the interval specification
add.interval.col("thalf.eff.obs",
                 FUN="pk.calc.thalf.eff.obs",
                 values=c(FALSE, TRUE),
                 desc="The effective half-life (as determined from the MRTobs)",
                 depends=c("mrt.obs"))
PKNCA.set.summary("thalf.eff.obs", business.geomean, business.geocv)
add.interval.col("thalf.eff.pred",
                 FUN="pk.calc.thalf.eff.pred",
                 values=c(FALSE, TRUE),
                 desc="The effective half-life (as determined from the MRTpred)",
                 depends=c("mrt.pred"))
PKNCA.set.summary("thalf.eff.pred", business.geomean, business.geocv)

#' Calculate the AUC percent extrapolated
#' 
#' @param auclast the area under the curve from time 0 to the last measurement
#'   above the limit of quantification
#' @param aucinf,aucinf.obs,aucinf.pred the area under the curve from time 0 to
#'   infinity
#' @return the numeric value of the AUC percent extrapolated
#' @export
pk.calc.aucpext <- function(auclast, aucinf) {
  if (auclast >= aucinf)
    warning("auclast should be less than aucinf")
  100*(1-auclast/aucinf)
}

#' @describeIn pk.calc.aucpext Compute the percent extrapolated AUCinf from the 
#'   observed Clast
#' @export
pk.calc.aucpext.obs <- function(auclast, aucinf.obs)
  pk.calc.aucpext(auclast, aucinf.obs)
#' @describeIn pk.calc.aucpext Compute the percent extrapolated AUCinf from the 
#'   observed Clast
#' @export
pk.calc.aucpext.pred <- function(auclast, aucinf.pred)
  pk.calc.aucpext(auclast, aucinf.pred)

## Add the columns to the interval specification
add.interval.col("aucpext.obs",
                 FUN="pk.calc.aucpext.obs",
                 values=c(FALSE, TRUE),
                 desc="Percent of the AUCinf that is extrapolated after Tlast calculated from the observed Clast",
                 depends=c("auclast", "aucinf.obs"))
PKNCA.set.summary("aucpext.obs", business.mean, business.sd)
add.interval.col("aucpext.pred",
                 FUN="pk.calc.aucpext.pred",
                 values=c(FALSE, TRUE),
                 desc="Percent of the AUCinf that is extrapolated after Tlast calculated from the predicted Clast",
                 depends=c("auclast", "aucinf.pred"))
PKNCA.set.summary("aucpext.pred", business.mean, business.sd)

#' Calculate the elimination rate (Kel)
#'
#' @param mrt,mrt.obs,mrt.pred the mean residence time
#' @return the numeric value of the elimination rate
#' @export
pk.calc.kel <- function(mrt)
  1/mrt
#' @describeIn pk.calc.kel Calculate Kel with observed Clast
#' @export
pk.calc.kel.obs <- function(mrt.obs)
  pk.calc.kel(mrt.obs)
#' @describeIn pk.calc.kel Calculate Kel with predicted Clast
#' @export
pk.calc.kel.pred <- function(mrt.pred)
  pk.calc.kel(mrt.pred)
## Add the columns to the interval specification
add.interval.col("kel.obs",
                 FUN="pk.calc.kel.obs",
                 values=c(FALSE, TRUE),
                 desc="Elimination rate (as calculated from the MRT with observed Clast)",
                 depends=c("mrt.obs"))
PKNCA.set.summary("kel.obs", business.geomean, business.geocv)
add.interval.col("kel.pred",
                 FUN="pk.calc.kel.pred",
                 values=c(FALSE, TRUE),
                 desc="Elimination rate (as calculated from the MRT with predicted Clast)",
                 depends=c("mrt.pred"))
PKNCA.set.summary("kel.pred", business.geomean, business.geocv)

#' Calculate the (observed oral) clearance
#' 
#' @param dose the dose administered
#' @param aucall The area under the concentration-time curve from 0 to the last 
#'   measurement above the limit of quantifiation (LOQ) plus the triangle to the
#'   first concentration below the LOQ.
#' @param auclast The area under the concentration-time curve from 0 to the last
#'   measurement above the LOQ.
#' @param aucinf,aucinf.obs,aucinf.pred The area under the concentration-time 
#'   curve from 0 to infinity (the next dose on a regular schedule at 
#'   steady-state)
#' @return the numeric value of the total (CL) or observed oral clearance (CL/F)
#' @details If \code{dose} is the same length as the other inputs, then the 
#'   output will be the same length as all of the inputs; the function assumes 
#'   that you are calculating for multiple intervals simultaneously.  If the 
#'   inputs other than \code{dose} are scalars and \code{dose} is a vector, then
#'   the function assumes multiple doses were given in a single interval, and 
#'   the sum of the \code{dose}s will be used for the calculation.
#' @references Gabrielsson J, Weiner D. "Section 2.5.1 Derivation of clearance."
#'   Pharmacokinetic & Pharmacodynamic Data Analysis: Concepts and Applications,
#'   4th Edition.  Stockholm, Sweden: Swedish Pharmaceutical Press, 2000. 86-7.
#' @export
pk.calc.cl <- function(dose, aucinf) {
  if (length(aucinf) == 1)
    dose <- sum(dose)
  dose/aucinf
}
#' @describeIn pk.calc.cl Compute the clearance from AUClast
#' @export
pk.calc.cl.last <- function(dose, auclast)
  pk.calc.cl(dose, auclast)
#' @describeIn pk.calc.cl Compute the clearance from AUCall
#' @export
pk.calc.cl.all <- function(dose, aucall)
  pk.calc.cl(dose, aucall)
#' @describeIn pk.calc.cl Compute the clearance from AUCinf (calculated from 
#'   observed Clast)
#' @export
pk.calc.cl.obs <- function(dose, aucinf.obs)
  pk.calc.cl(dose, aucinf.obs)
#' @describeIn pk.calc.cl Compute the clearance from AUCinf (calculated from
#'   predicted Clast)
#' @export
pk.calc.cl.pred <- function(dose, aucinf.pred)
  pk.calc.cl(dose, aucinf.pred)

## Add the columns to the interval specification
add.interval.col("cl.last",
                 FUN="pk.calc.cl.last",
                 values=c(FALSE, TRUE),
                 desc="Clearance or observed oral clearance calculated to Clast",
                 depends=c("auclast"))
PKNCA.set.summary("cl.last", business.geomean, business.geocv)
add.interval.col("cl.all",
                 FUN="pk.calc.cl.all",
                 values=c(FALSE, TRUE),
                 desc="Clearance or observed oral clearance calculated with AUCall",
                 depends=c("aucall"))
PKNCA.set.summary("cl.all", business.geomean, business.geocv)
add.interval.col("cl.obs",
                 FUN="pk.calc.cl.obs",
                 values=c(FALSE, TRUE),
                 desc="Clearance or observed oral clearance calculated with observed Clast",
                 depends=c("aucinf.obs"))
PKNCA.set.summary("cl.obs", business.geomean, business.geocv)
add.interval.col("cl.pred",
                 FUN="pk.calc.cl.pred",
                 values=c(FALSE, TRUE),
                 desc="Clearance or observed oral clearance calculated with predicted Clast",
                 depends=list("aucinf.pred"))
PKNCA.set.summary("cl.pred", business.geomean, business.geocv)

#' Calculate the absolute (or relative) bioavailability
#'
#' @param dose1 The dose administered in route or method 1
#' @param dose2 The dose administered in route or method 2
#' @param auc1 The AUC from 0 to infinity or 0 to tau administered in
#'   route or method 1
#' @param auc2 The AUC from 0 to infinity or 0 to tau administered in
#'   route or method 2
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
#' @param aucinf,aucinf.obs,aucinf.pred the AUC from 0 to infinity or
#'   0 to tau at steady-state
#' @param auclast the AUC from 0 to the last concentration above the
#'   limit of quantification (LOQ)
#' @param aumcinf,aumcinf.obs,aumcinf.pred the AUMC from 0 to infinity
#'   or 0 to tau at steady-state
#' @param aumclast the AUMC from 0 to the last concentration above the
#'   LOQ
#' @return the numeric value of the mean residence time
#' @export
pk.calc.mrt <- function(aucinf, aumcinf)
  aumcinf/aucinf
#' @describeIn pk.calc.mrt Calcuate the mean residence time (MRT)
#'   using observed Clast
#' @export
pk.calc.mrt.obs <- function(aucinf.obs, aumcinf.obs)
  pk.calc.mrt(aucinf.obs, aumcinf.obs)
#' @describeIn pk.calc.mrt Calcuate the mean residence time (MRT)
#'   using predicted Clast
#' @export
pk.calc.mrt.pred <- function(aucinf.pred, aumcinf.pred)
  pk.calc.mrt(aucinf.pred, aumcinf.pred)
## Add the columns to the interval specification
add.interval.col("mrt.obs",
                 FUN="pk.calc.mrt.obs",
                 values=c(FALSE, TRUE),
                 desc="The mean residence time to infinity using observed Clast",
                 depends=c("aucinf.obs", "aumcinf.obs"))
PKNCA.set.summary("mrt.obs", business.geomean, business.geocv)
add.interval.col("mrt.pred",
                 FUN="pk.calc.mrt.pred",
                 values=c(FALSE, TRUE),
                 desc="The mean residence time to infinity using predicted Clast",
                 depends=c("aucinf.pred", "aumcinf.pred"))
PKNCA.set.summary("mrt.pred", business.geomean, business.geocv)

#' @describeIn pk.calc.mrt Calculate the mean residence time (MRT) to the last
#'   concentration above the limit of quantification
#' @export
pk.calc.mrt.last <- function(auclast, aumclast)
  pk.calc.mrt(auclast, aumclast)
## Add the column to the interval specification
add.interval.col("mrt.last",
                 FUN="pk.calc.mrt.last",
                 values=c(FALSE, TRUE),
                 desc="The mean residence time to the last observed concentration above the LOQ",
                 depends=list("auclast", "aumclast"))
PKNCA.set.summary("mrt.last", business.geomean, business.geocv)

#' Calculate the terminal volume of distribution (Vz)
#'
#' @param cl,cl.obs,cl.pred the clearance (or apparent observed clearance)
#' @param lambda.z the elimination rate
#' @export
pk.calc.vz <- function(cl, lambda.z) {
  ## Ensure that cl is either a scalar or the same length as AUC
  ## (more complex repeating patterns while valid for general R are
  ## likely errors here).
  if (!(length(cl) %in% c(1, length(lambda.z))) |
      !(length(lambda.z) %in% c(1, length(cl))))
    stop("'cl' and 'lambda.z' must be the same length")
  cl/lambda.z
}
#' @describeIn pk.calc.vz Calculate Vz using observed Clast
#' @export
pk.calc.vz.obs <- function(cl.obs, lambda.z)
  pk.calc.vz(cl.obs, lambda.z)
#' @describeIn pk.calc.vz Calculate Vz using predicted Clast
#' @export
pk.calc.vz.pred <- function(cl.pred, lambda.z)
  pk.calc.vz(cl.pred, lambda.z)
## Add the columns to the interval specification
add.interval.col("vz.obs",
                 FUN="pk.calc.vz.obs",
                 values=c(FALSE, TRUE),
                 desc="The terminal volume of distribution using observed Clast",
                 depends=c("cl.obs", "lambda.z"))
PKNCA.set.summary("vz.obs", business.geomean, business.geocv)
add.interval.col("vz.pred",
                 FUN="pk.calc.vz.pred",
                 values=c(FALSE, TRUE),
                 desc="The terminal volume of distribution using predicted Clast",
                 depends=c("cl.pred", "lambda.z"))
PKNCA.set.summary("vz.pred", business.geomean, business.geocv)

#' Calculate the steady-state volume of distribution (Vss)
#'
#' @param cl,cl.obs,cl.pred the clearance
#' @param mrt,mrt.obs,mrt.pred the mean residence time
#' @return the volume of distribution at steady-state
#' @export
pk.calc.vss <- function(cl, mrt)
  cl*mrt
#' @describeIn pk.calc.vss Vss calculation using observed Clast
#' @export
pk.calc.vss.obs <- function(cl.obs, mrt.obs)
  pk.calc.vss(cl.obs, mrt.obs)
#' @describeIn pk.calc.vss Vss calculation using predicted Clast
#' @export
pk.calc.vss.pred <- function(cl.pred, mrt.pred)
  pk.calc.vss(cl.pred, mrt.pred)
# Add the columns to the interval specification
add.interval.col("vss.obs",
                 FUN="pk.calc.vss.obs",
                 values=c(FALSE, TRUE),
                 desc="The steady-state volume of distribution using observed Clast",
                 depends=c("cl.obs", "mrt.obs"))
PKNCA.set.summary("vss.obs", business.geomean, business.geocv)
add.interval.col("vss.pred",
                 FUN="pk.calc.vss.pred",
                 values=c(FALSE, TRUE),
                 desc="The steady-state volume of distribution using predicted Clast",
                 depends=c("cl.pred", "mrt.pred"))
PKNCA.set.summary("vss.pred", business.geomean, business.geocv)

#' Calculate the volume of distribution (Vd) or observed volume of distribution 
#' (Vd/F)
#' 
#' @param dose One or more doses given during an interval
#' @param aucinf,aucinf.obs,aucinf.pred Area under the curve to
#'   infinity (either predicted or observed).
#' @param lambda.z Elimination rate constant
#' @details If \code{dose} is the same length as the other inputs,
#'   then the output will be the same length as all of the inputs; the
#'   function assumes that you are calculating for multiple intervals
#'   simultaneously.  If the inputs other than \code{dose} are scalars
#'   and \code{dose} is a vector, then the function assumes multiple
#'   doses were given in a single interval, and the sum of the
#'   \code{dose}s will be used for the calculation.
#' @return The observed volume of distribution
#' @export
pk.calc.vd <- function(dose, aucinf, lambda.z) {
  if (length(aucinf) == 1 &
      length(lambda.z) == 1) {
    dose <- sum(dose)
  }
  dose/(aucinf * lambda.z)
}
#' @describeIn pk.calc.vd Compute the Vd with observed Clast
#' @export
pk.calc.vd.obs <- function(dose, aucinf.obs, lambda.z)
  pk.calc.vd(dose, aucinf.obs, lambda.z)
#' @describeIn pk.calc.vd Compute the Vd with predicted Clast
#' @export
pk.calc.vd.pred <- function(dose, aucinf.pred, lambda.z)
  pk.calc.vd(dose, aucinf.pred, lambda.z)
## Add the columns to the interval specification
add.interval.col("vd.obs",
                 FUN="pk.calc.vd.obs",
                 values=c(FALSE, TRUE),
                 desc="Apparent observed volume of distribution calculated with observed Clast",
                 depends=c("aucinf.obs", "lambda.z"))
PKNCA.set.summary("vd.obs", business.geomean, business.geocv)
add.interval.col("vd.pred",
                 FUN="pk.calc.vd.pred",
                 values=c(FALSE, TRUE),
                 desc="Apparent observed volume of distribution calculated with predicted Clast",
                 depends=c("aucinf.pred", "lambda.z"))
PKNCA.set.summary("vd.pred", business.geomean, business.geocv)

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
#'   match, then \code{NA}
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
  mask.increase <- c(conc[-1] > conc[-length(conc)], FALSE)
  if (any(mask.increase)) {
    time[mask.increase][1]
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

#' Determine the degree of fluctuation
#' 
#' @param cmax The maximum observed concentration
#' @param cmin The minimum observed concentration
#' @param cav The average concentration in the interval
#' @return The degree of fluctuation around the average concentration.
#' @export
pk.calc.deg.fluc <- function(cmax, cmin, cav) {
  100*(cmax - cmin)/cav
}
add.interval.col("deg.fluc",
                 FUN="pk.calc.deg.fluc",
                 desc="Degree of fluctuation",
                 depends=c("cmax", "cmin", "cav"))
PKNCA.set.summary("deg.fluc", business.mean, business.sd)

#' Determine the PK swing
#' 
#' @param cmax The maximum observed concentration
#' @param cmin The minimum observed concentration
#' @return The swing above the minimum concentration.  If \code{cmin} is zero,
#'   then the result is infinity.
#' @export
pk.calc.swing <- function(cmax, cmin) {
  if (cmin > 0) {
    100*(cmax - cmin)/cmin
  } else {
    Inf
  }
}
add.interval.col("swing",
                 FUN="pk.calc.swing",
                 desc="Swing relative to Cmin",
                 depends=c("cmax", "cmin"))
PKNCA.set.summary("swing", business.mean, business.sd)

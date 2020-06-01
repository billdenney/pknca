## Calculate the simple parameters for PK.

#' Calculate the adjusted r-squared value
#'
#' @param r.sq The r-squared value
#' @param n The number of points
#' @return The numeric adjusted r-squared value
#' @export
adj.r.squared <- function(r.sq, n) {
  if (n <= 2) {
    warning("n must be > 2 for adj.r.squared")
    structure(NA_real_, exclude="n must be > 2")
  } else {
    1-(1-r.sq)*(n-1)/(n-2)
  }
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
PKNCA.set.summary(
  name="cmax",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

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
PKNCA.set.summary(
  name="cmin",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

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
                         first.tmax=NULL,
                         check=TRUE) {
  first.tmax <- PKNCA.choose.option(name="first.tmax", value=first.tmax, options=options)
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
PKNCA.set.summary(
  name="tmax",
  description="median and range",
  point=business.median,
  spread=business.range
)

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
PKNCA.set.summary(
  name="tlast",
  description="median and range",
  point=business.median,
  spread=business.range
)

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
PKNCA.set.summary(
  name="tfirst",
  description="median and range",
  point=business.median,
  spread=business.range
)

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
PKNCA.set.summary(
  name="clast.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Calculate the effective half-life
#' 
#' @details thalf.eff is \code{log(2)*mrt}.
#'
#' @param mrt the mean residence time to infinity
#' @return the numeric value of the effective half-life
#' @export
pk.calc.thalf.eff <- function(mrt)
  log(2)*mrt
## Add the columns to the interval specification
add.interval.col("thalf.eff.obs",
                 FUN="pk.calc.thalf.eff",
                 values=c(FALSE, TRUE),
                 desc="The effective half-life (as determined from the MRTobs)",
                 formalsmap=list(mrt="mrt.obs"),
                 depends=c("mrt.obs"))
PKNCA.set.summary(
  name="thalf.eff.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("thalf.eff.pred",
                 FUN="pk.calc.thalf.eff",
                 values=c(FALSE, TRUE),
                 desc="The effective half-life (as determined from the MRTpred)",
                 formalsmap=list(mrt="mrt.pred"),
                 depends=c("mrt.pred"))
PKNCA.set.summary(
  name="thalf.eff.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("thalf.eff.last",
                 FUN="pk.calc.thalf.eff",
                 values=c(FALSE, TRUE),
                 desc="The effective half-life (as determined from the MRTlast)",
                 formalsmap=list(mrt="mrt.last"),
                 depends=c("mrt.last"))
PKNCA.set.summary(
  name="thalf.eff.last",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("thalf.eff.iv.obs",
                 FUN="pk.calc.thalf.eff",
                 values=c(FALSE, TRUE),
                 desc="The effective half-life (as determined from the intravenous MRTobs)",
                 formalsmap=list(mrt="mrt.iv.obs"),
                 depends=c("mrt.iv.obs"))
PKNCA.set.summary(
  name="thalf.eff.iv.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("thalf.eff.iv.pred",
                 FUN="pk.calc.thalf.eff",
                 values=c(FALSE, TRUE),
                 desc="The effective half-life (as determined from the intravenous MRTpred)",
                 formalsmap=list(mrt="mrt.iv.pred"),
                 depends=c("mrt.iv.pred"))
PKNCA.set.summary(
  name="thalf.eff.iv.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("thalf.eff.iv.last",
                 FUN="pk.calc.thalf.eff",
                 values=c(FALSE, TRUE),
                 desc="The effective half-life (as determined from the intravenous MRTlast)",
                 formalsmap=list(mrt="mrt.iv.last"),
                 depends=c("mrt.iv.last"))
PKNCA.set.summary(
  name="thalf.eff.iv.last",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Calculate the AUC percent extrapolated
#' 
#' @details aucpext is \code{100*(1-auclast/aucinf)}.
#' 
#' @param auclast the area under the curve from time 0 to the last 
#'   measurement above the limit of quantification
#' @param aucinf the area under the curve from time 0 to infinity
#' @return The numeric value of the AUC percent extrapolated or 
#'   \code{NA_real_} if any of the following are true 
#'   \code{is.na(aucinf)}, \code{is.na(auclast)}, \code{aucinf <= 0},
#'   or \code{auclast <= 0}.
#' @export
pk.calc.aucpext <- function(auclast, aucinf) {
  scalar_auclast <- length(auclast) == 1
  scalar_aucinf <- length(aucinf) == 1
  if (scalar_auclast | scalar_aucinf) {
    # no length checking needs to occur
  } else if ((!scalar_auclast & !scalar_aucinf) &
             length(auclast) != length(aucinf)) {
    stop("auclast and aucinf must either be a scalar or the same length.")
  }
  ret <- rep(NA_real_, max(c(length(auclast), length(aucinf))))
  mask_na <-
    is.na(auclast) |
    is.na(aucinf)
  mask_negative <-
    !mask_na &
    (aucinf <= 0 |
       auclast <= 0)
  mask_greater <-
    !mask_na &
    (auclast >= aucinf)
  mask_calc <- !mask_na
  if (any(mask_greater))
    warning("aucpext is typically only calculated when aucinf is greater than auclast.")
  if (any(mask_negative))
    warning("aucpext is typically only calculated when both aucinf and auclast are positive.")
  ret[mask_calc] <-
    100*(1-auclast[mask_calc]/aucinf[mask_calc])
  ret
}

## Add the columns to the interval specification
add.interval.col("aucpext.obs",
                 FUN="pk.calc.aucpext",
                 values=c(FALSE, TRUE),
                 desc="Percent of the AUCinf that is extrapolated after Tlast calculated from the observed Clast",
                 formalsmap=list(aucinf="aucinf.obs"),
                 depends=c("auclast", "aucinf.obs"))
PKNCA.set.summary(
  name="aucpext.obs",
  description="arithmetic mean and standard deviation",
  point=business.mean,
  spread=business.sd
)
add.interval.col("aucpext.pred",
                 FUN="pk.calc.aucpext",
                 values=c(FALSE, TRUE),
                 desc="Percent of the AUCinf that is extrapolated after Tlast calculated from the predicted Clast",
                 formalsmap=list(aucinf="aucinf.pred"),
                 depends=c("auclast", "aucinf.pred"))
PKNCA.set.summary(
  name="aucpext.pred",
  description="arithmetic mean and standard deviation",
  point=business.mean,
  spread=business.sd
)

#' Calculate the elimination rate (Kel)
#'
#' @param kel is \code{1/mrt}, not to be confused with lambda.z.
#'
#' @param mrt the mean residence time
#' @return the numeric value of the elimination rate
#' @export
pk.calc.kel <- function(mrt)
  1/mrt
## Add the columns to the interval specification
add.interval.col("kel.obs",
                 FUN="pk.calc.kel",
                 values=c(FALSE, TRUE),
                 desc="Elimination rate (as calculated from the MRT with observed Clast)",
                 formalsmap=list(mrt="mrt.obs"),
                 depends=c("mrt.obs"))
PKNCA.set.summary(
  name="kel.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("kel.pred",
                 FUN="pk.calc.kel",
                 values=c(FALSE, TRUE),
                 desc="Elimination rate (as calculated from the MRT with predicted Clast)",
                 formalsmap=list(mrt="mrt.pred"),
                 depends=c("mrt.pred"))
PKNCA.set.summary(
  name="kel.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("kel.last",
                 FUN="pk.calc.kel",
                 values=c(FALSE, TRUE),
                 desc="Elimination rate (as calculated from the MRT using AUClast)",
                 formalsmap=list(mrt="mrt.last"),
                 depends=c("mrt.last"))
PKNCA.set.summary(
  name="kel.last",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("kel.iv.obs",
                 FUN="pk.calc.kel",
                 values=c(FALSE, TRUE),
                 desc="Elimination rate (as calculated from the intravenous MRTobs)",
                 formalsmap=list(mrt="mrt.iv.obs"),
                 depends=c("mrt.iv.obs"))
PKNCA.set.summary(
  name="kel.iv.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("kel.iv.pred",
                 FUN="pk.calc.kel",
                 values=c(FALSE, TRUE),
                 desc="Elimination rate (as calculated from the intravenous MRTpred)",
                 formalsmap=list(mrt="mrt.iv.pred"),
                 depends=c("mrt.iv.pred"))
PKNCA.set.summary(
  name="kel.iv.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("kel.iv.last",
                 FUN="pk.calc.kel",
                 values=c(FALSE, TRUE),
                 desc="Elimination rate (as calculated from the intravenous MRTlast)",
                 formalsmap=list(mrt="mrt.iv.last"),
                 depends=c("mrt.iv.last"))
PKNCA.set.summary(
  name="kel.iv.last",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Calculate the (observed oral) clearance
#' 
#' @details cl is \code{dose/auc}.
#' 
#' @param dose the dose administered
#' @param auc The area under the concentration-time curve.
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
pk.calc.cl <- function(dose, auc) {
  if (length(auc) == 1) {
    dose <- sum(dose)
  }
  ret <- dose/auc
  mask_zero <- !is.na(auc) & (auc <= 0)
  if (any(mask_zero)) {
    ret[mask_zero] <- NA_real_
  }
  ret
}

## Add the columns to the interval specification
add.interval.col("cl.last",
                 FUN="pk.calc.cl",
                 values=c(FALSE, TRUE),
                 desc="Clearance or observed oral clearance calculated to Clast",
                 formalsmap=list(auc="auclast"),
                 depends=c("auclast"))
PKNCA.set.summary(
  name="cl.last",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("cl.all",
                 FUN="pk.calc.cl",
                 values=c(FALSE, TRUE),
                 desc="Clearance or observed oral clearance calculated with AUCall",
                 formalsmap=list(auc="aucall"),
                 depends=c("aucall"))
PKNCA.set.summary(
  name="cl.all",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("cl.obs",
                 FUN="pk.calc.cl",
                 values=c(FALSE, TRUE),
                 desc="Clearance or observed oral clearance calculated with observed Clast",
                 formalsmap=list(auc="aucinf.obs"),
                 depends=c("aucinf.obs"))
PKNCA.set.summary(
  name="cl.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("cl.pred",
                 FUN="pk.calc.cl",
                 values=c(FALSE, TRUE),
                 desc="Clearance or observed oral clearance calculated with predicted Clast",
                 formalsmap=list(auc="aucinf.pred"),
                 depends=list("aucinf.pred"))
PKNCA.set.summary(
  name="cl.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Calculate the absolute (or relative) bioavailability
#' 
#' @details f is \code{(auc2/dose2)/(auc1/dose1)}.
#'
#' @param dose1 The dose administered in route or method 1
#' @param dose2 The dose administered in route or method 2
#' @param auc1 The AUC from 0 to infinity or 0 to tau administered in
#'   route or method 1
#' @param auc2 The AUC from 0 to infinity or 0 to tau administered in
#'   route or method 2
#' @export
pk.calc.f <- function(dose1, auc1, dose2, auc2) {
  ret <- (auc2/dose2)/(auc1/dose1)
  mask_zero <-
    is.na(auc1)  | (auc1 <= 0) |
    is.na(dose2) | (dose2 <= 0) |
    is.na(dose1) | (dose1 <= 0)
  if (any(mask_zero)) {
    ret[mask_zero] <- NA_real_
  }
  ret
}
add.interval.col("f",
                 FUN="pk.calc.f",
                 values=c(FALSE, TRUE),
                 desc="Bioavailability or relative bioavailability",
                 depends=c())
PKNCA.set.summary(
  name="f",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Calculate the mean residence time (MRT) for single-dose data or linear
#' multiple-dose data.
#'
#' @details mrt is \code{aumc/auc - duration.dose/2} where \code{duration.dose =
#'   0} for oral administration.
#' 
#' @param auc the AUC from 0 to infinity or 0 to tau
#' @param aumc the AUMC from 0 to infinity or 0 to tau
#' @param duration.dose The duration of the dose (usually an infusion
#'   duration for an IV infusion)
#' @return the numeric value of the mean residence time
#' @seealso \code{\link{pk.calc.mrt.md}}
#' @export
pk.calc.mrt <- function(auc, aumc) {
  pk.calc.mrt.iv(auc, aumc, duration.dose=0)
}
## Add the columns to the interval specification
add.interval.col("mrt.obs",
                 FUN="pk.calc.mrt",
                 values=c(FALSE, TRUE),
                 desc="The mean residence time to infinity using observed Clast",
                 formalsmap=list(auc="aucinf.obs", aumc="aumcinf.obs"),
                 depends=c("aucinf.obs", "aumcinf.obs"))
PKNCA.set.summary(
  name="mrt.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("mrt.pred",
                 FUN="pk.calc.mrt",
                 values=c(FALSE, TRUE),
                 desc="The mean residence time to infinity using predicted Clast",
                 formalsmap=list(auc="aucinf.pred", aumc="aumcinf.pred"),
                 depends=c("aucinf.pred", "aumcinf.pred"))
PKNCA.set.summary(
  name="mrt.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("mrt.last",
                 FUN="pk.calc.mrt",
                 values=c(FALSE, TRUE),
                 desc="The mean residence time to the last observed concentration above the LOQ",
                 formalsmap=list(auc="auclast", aumc="aumclast"),
                 depends=list("auclast", "aumclast"))
PKNCA.set.summary(
  name="mrt.last",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' @describeIn pk.calc.mrt MRT for an IV infusion
#' @export
pk.calc.mrt.iv <- function(auc, aumc, duration.dose) {
  ret <- aumc/auc - duration.dose/2
  mask_zero <- is.na(auc) | auc <= 0
  if (any(mask_zero)) {
    ret[mask_zero] <- NA_real_
  }
  ret
}
## Add the columns to the interval specification
add.interval.col("mrt.iv.obs",
                 FUN="pk.calc.mrt.iv",
                 values=c(FALSE, TRUE),
                 desc="The mean residence time to infinity using observed Clast correcting for dosing duration",
                 formalsmap=list(auc="aucinf.obs", aumc="aumcinf.obs"),
                 depends=c("aucinf.obs", "aumcinf.obs"))
PKNCA.set.summary(
  name="mrt.iv.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("mrt.iv.pred",
                 FUN="pk.calc.mrt.iv",
                 values=c(FALSE, TRUE),
                 desc="The mean residence time to infinity using predicted Clast correcting for dosing duration",
                 formalsmap=list(auc="aucinf.pred", aumc="aumcinf.pred"),
                 depends=c("aucinf.pred", "aumcinf.pred"))
PKNCA.set.summary(
  name="mrt.iv.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("mrt.iv.last",
                 FUN="pk.calc.mrt.iv",
                 values=c(FALSE, TRUE),
                 desc="The mean residence time to the last observed concentration above the LOQ correcting for dosing duration",
                 formalsmap=list(auc="auclast", aumc="aumclast"),
                 depends=list("auclast", "aumclast"))
PKNCA.set.summary(
  name="mrt.iv.last",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Calculate the mean residence time (MRT) for multiple-dose data with nonlinear
#' kinetics.
#'
#' @details mrt.md is \code{aumctau/auctau + tau*(aucinf-auctau)/auctau} and
#' should only be used for multiple dosing with equal intervals between doses.
#' 
#' @param auctau the AUC from time 0 to the end of the dosing interval 
#'   (tau).
#' @param aumctau the AUMC from time 0 to the end of the dosing interval
#'   (tau).
#' @param aucinf the AUC from time 0 to infinity (typically using 
#'   single-dose data)
#' @param tau the dosing interval
#' @details Note that if \code{aucinf == auctau} (as would be the 
#'   assumption with linear kinetics), the equation becomes the same as 
#'   the single-dose MRT.
#' @seealso \code{\link{pk.calc.mrt}}
#' @export
pk.calc.mrt.md <- function(auctau, aumctau, aucinf, tau) {
  ret <- aumctau/auctau + tau*(aucinf-auctau)/auctau
  mask_zero <- is.na(auctau) | auctau <= 0
  if (any(mask_zero)) {
    ret[mask_zero] <- NA_real_
  }
  ret
}
add.interval.col("mrt.md.obs",
                 FUN="pk.calc.mrt.md",
                 values=c(FALSE, TRUE),
                 desc="The mean residence time with multiple dosing and nonlinear kinetics using observed Clast",
                 formalsmap=list(auctau="auclast", aumctau="aumclast", aucinf="aucinf.obs"),
                 depends=c("auclast", "aumclast", "aucinf.obs"))
PKNCA.set.summary(
  name="mrt.md.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("mrt.md.pred",
                 FUN="pk.calc.mrt.md",
                 values=c(FALSE, TRUE),
                 desc="The mean residence time with multiple dosing and nonlinear kinetics using predicted Clast",
                 formalsmap=list(auctau="auclast", aumctau="aumclast", aucinf="aucinf.pred"),
                 depends=c("auclast", "aumclast", "aucinf.pred"))
PKNCA.set.summary(
  name="mrt.md.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Calculate the terminal volume of distribution (Vz)
#'
#' @details vz is \code{cl/lambda.z}.
#'
#' @param cl the clearance (or apparent observed clearance)
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
## Add the columns to the interval specification
add.interval.col("vz.obs",
                 FUN="pk.calc.vz",
                 values=c(FALSE, TRUE),
                 desc="The terminal volume of distribution using observed Clast",
                 formalsmap=list(cl="cl.obs"),
                 depends=c("cl.obs", "lambda.z"))
PKNCA.set.summary(
  name="vz.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("vz.pred",
                 FUN="pk.calc.vz",
                 values=c(FALSE, TRUE),
                 desc="The terminal volume of distribution using predicted Clast",
                 formalsmap=list(cl="cl.pred"),
                 depends=c("cl.pred", "lambda.z"))
PKNCA.set.summary(
  name="vz.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Calculate the steady-state volume of distribution (Vss)
#'
#' @details vss is \code{cl*mrt}.
#' @param cl the clearance
#' @param mrt the mean residence time
#' @return the volume of distribution at steady-state
#' @export
pk.calc.vss <- function(cl, mrt)
  cl*mrt
# Add the columns to the interval specification
add.interval.col("vss.obs",
                 FUN="pk.calc.vss",
                 values=c(FALSE, TRUE),
                 desc="The steady-state volume of distribution using observed Clast",
                 formalsmap=list(cl="cl.obs", mrt="mrt.obs"),
                 depends=c("cl.obs", "mrt.obs"))
PKNCA.set.summary(
  name="vss.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("vss.pred",
                 FUN="pk.calc.vss",
                 values=c(FALSE, TRUE),
                 desc="The steady-state volume of distribution using predicted Clast",
                 formalsmap=list(cl="cl.pred", mrt="mrt.pred"),
                 depends=c("cl.pred", "mrt.pred"))
PKNCA.set.summary(
  name="vss.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("vss.last",
                 FUN="pk.calc.vss",
                 values=c(FALSE, TRUE),
                 desc="The steady-state volume of distribution calculating through Tlast",
                 formalsmap=list(cl="cl.last", mrt="mrt.last"),
                 depends=c("cl.last", "mrt.last"))
PKNCA.set.summary(
  name="vss.last",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("vss.iv.obs",
                 FUN="pk.calc.vss",
                 values=c(FALSE, TRUE),
                 desc="The steady-state volume of distribution with intravenous infusion using observed Clast",
                 formalsmap=list(cl="cl.obs", mrt="mrt.iv.obs"),
                 depends=c("cl.obs", "mrt.iv.obs"))
PKNCA.set.summary(
  name="vss.iv.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("vss.iv.pred",
                 FUN="pk.calc.vss",
                 values=c(FALSE, TRUE),
                 desc="The steady-state volume of distribution with intravenous infusion using predicted Clast",
                 formalsmap=list(cl="cl.pred", mrt="mrt.iv.pred"),
                 depends=c("cl.pred", "mrt.iv.pred"))
PKNCA.set.summary(
  name="vss.iv.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("vss.iv.last",
                 FUN="pk.calc.vss",
                 values=c(FALSE, TRUE),
                 desc="The steady-state volume of distribution with intravenous infusion calculating through Tlast",
                 formalsmap=list(cl="cl.last", mrt="mrt.iv.last"),
                 depends=c("cl.last", "mrt.iv.last"))
PKNCA.set.summary(
  name="vss.iv.last",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

add.interval.col("vss.md.obs",
                 FUN="pk.calc.vss",
                 values=c(FALSE, TRUE),
                 desc="The steady-state volume of distribution for nonlinear multiple-dose data using observed Clast",
                 formalsmap=list(cl="cl.last", mrt="mrt.md.obs"),
                 depends=c("cl.last", "mrt.md.obs"))
PKNCA.set.summary(
  name="vss.md.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("vss.md.pred",
                 FUN="pk.calc.vss",
                 values=c(FALSE, TRUE),
                 desc="The steady-state volume of distribution for nonlinear multiple-dose data using predicted Clast",
                 formalsmap=list(cl="cl.last", mrt="mrt.md.pred"),
                 depends=c("cl.last", "mrt.md.pred"))
PKNCA.set.summary(
  name="vss.md.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Calculate the volume of distribution (Vd) or observed volume of
#' distribution (Vd/F)
#' 
#' @details vd is \code{dose/(aucinf * lambda.z)}.
#' 
#' @param dose One or more doses given during an interval
#' @param aucinf Area under the curve to infinity (either predicted or
#'   observed).
#' @param lambda.z Elimination rate constant
#' @details If \code{dose} is the same length as the other inputs, then
#'   the output will be the same length as all of the inputs; the 
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
  ret <- dose/(aucinf * lambda.z)
  mask_zero <-
    is.na(aucinf) | aucinf <= 0 |
    is.na(lambda.z) | lambda.z <= 0
  if (any(mask_zero)) {
    ret[mask_zero] <- NA_real_
  }
  ret
}
## Add the columns to the interval specification
add.interval.col("vd.obs",
                 FUN="pk.calc.vd",
                 values=c(FALSE, TRUE),
                 desc="Apparent observed volume of distribution calculated with observed Clast",
                 formalsmap=list(aucinf="aucinf.obs"),
                 depends=c("aucinf.obs", "lambda.z"))
PKNCA.set.summary(
  name="vd.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("vd.pred",
                 FUN="pk.calc.vd",
                 values=c(FALSE, TRUE),
                 desc="Apparent observed volume of distribution calculated with predicted Clast",
                 formalsmap=list(aucinf="aucinf.pred"),
                 depends=c("aucinf.pred", "lambda.z"))
PKNCA.set.summary(
  name="vd.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Calculate the average concentration during an interval.
#' 
#' @details cav is \code{auclast/(end-start)}.
#'
#' @param auclast The area under the curve during the interval
#' @param start The starting time of the interval
#' @param end The ending time of the interval
#' @return The Cav (average concentration during the interval)
#' @export
pk.calc.cav <- function(auclast, start, end) {
  ret <- auclast/(end-start)
  mask_zero <- is.na(end) | is.na(start) | end == start
  if (any(mask_zero)) {
    ret[mask_zero] <- NA_real_
  }
  ret
}
add.interval.col("cav",
                 FUN="pk.calc.cav",
                 values=c(FALSE, TRUE),
                 desc="The average concentration during an interval",
                 depends="auclast")
PKNCA.set.summary(
  name="cav",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Determine the trough (predose) concentration
#'
#' @param conc Observed concentrations during the interval
#' @param time Times of \code{conc} observations
#' @param end End time of the interval
#' @return The concentration when \code{time == end}.  If none
#'   match, then \code{NA}
#' @export
pk.calc.ctrough <- function(conc, time, end) {
  check.conc.time(conc, time)
  mask_end <- time %in% end
  if (sum(mask_end) == 1) {
    conc[mask_end]
  } else if (sum(mask_end) == 0) {
    NA_real_
  } else {
    # This should be impossible as check.conc.time should catch
    # duplicates.
    stop("More than one time matches the starting time.  Please report this as a bug with a reproducible example.") # nocov
  }
}
add.interval.col("ctrough",
                 FUN="pk.calc.ctrough",
                 values=c(FALSE, TRUE),
                 desc="The trough (predose) concentration",
                 depends=c())
PKNCA.set.summary(
  name="ctrough",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Determine the peak-to-trough ratio
#'
#' @details ptr is \code{cmax/ctrough}.
#'
#' @param cmax The maximum observed concentration
#' @param ctrough The last concentration in an interval
#' @return The ratio of cmax to ctrough (if ctrough == 0, NA)
#' @export
pk.calc.ptr <- function(cmax, ctrough) {
  ret <- cmax/ctrough
  ret[ctrough %in% 0] <- NA_real_
  ret
}
add.interval.col("ptr",
                 FUN="pk.calc.ptr",
                 values=c(FALSE, TRUE),
                 desc="Peak-to-Trough ratio (fraction)",
                 depends=c("cmax", "ctrough"))
PKNCA.set.summary(
  name="ptr",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

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
    NA_real_
  }
}
add.interval.col("tlag",
                 FUN="pk.calc.tlag",
                 values=c(FALSE, TRUE),
                 desc="Lag time",
                 depends=c())
PKNCA.set.summary(
  name="tlag",
  description="median and range",
  point=business.median,
  spread=business.range
)

#' Determine the degree of fluctuation
#' 
#' @details deg.fluc is \code{100*(cmax - cmin)/cav}.
#'
#' @param cmax The maximum observed concentration
#' @param cmin The minimum observed concentration
#' @param cav The average concentration in the interval
#' @return The degree of fluctuation around the average concentration.
#' @export
pk.calc.deg.fluc <- function(cmax, cmin, cav) {
  ret <- 100*(cmax - cmin)/cav
  mask_zero <- is.na(cav) | cav == 0
  if (any(mask_zero)) {
    ret[mask_zero] <- NA_real_
  }
  ret
}
add.interval.col("deg.fluc",
                 FUN="pk.calc.deg.fluc",
                 desc="Degree of fluctuation",
                 depends=c("cmax", "cmin", "cav"))
PKNCA.set.summary(
  name="deg.fluc",
  description="arithmetic mean and standard deviation",
  point=business.mean,
  spread=business.sd
)

#' Determine the PK swing
#' 
#' @details swing is \code{100*(cmax - cmin)/cmin}.
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
PKNCA.set.summary(
  name="swing",
  description="arithmetic mean and standard deviation",
  point=business.mean,
  spread=business.sd
)

#' Determine the concentration at the end of infusion
#' 
#' @param conc Concentration measured
#' @param time Time of concentration measurement
#' @param duration.dose The duration for the dosing administration 
#'   (typically from IV infusion)
#' @param check Run \code{\link{check.conc.time}}?
#' @return The concentration at the end of the infusion, \code{NA} if
#'   duration.dose is \code{NA}, or \code{NA} if all \code{time != duration.dose}
#' @export
pk.calc.ceoi <- function(conc, time, duration.dose=NA, check=TRUE) {
  if (check)
    check.conc.time(conc=conc, time=time)
  if (is.na(duration.dose)) {
    NA_real_
  } else if (all(time != duration.dose)) {
    NA_real_
  } else {
    conc[time == duration.dose][1]
  }
}
add.interval.col("ceoi",
                 FUN="pk.calc.ceoi",
                 desc="Concentration at the end of infusion",
                 depends=c())
PKNCA.set.summary(
  name="ceoi",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

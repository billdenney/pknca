# Calculate the simple parameters for PK.

#' Calculate the adjusted r-squared value
#'
#' @param r.sq The r-squared value
#' @param n The number of points
#' @return The numeric adjusted r-squared value
#' @export
adj.r.squared <- function(r.sq, n) {
  if (n <= 2) {
    rlang::warn(
      message = "n must be > 2 for adj.r.squared",
      class = "pknca_adjr2_2points"
    )
    structure(NA_real_, exclude="n must be > 2")
  } else {
    1-(1-r.sq)*(n-1)/(n-2)
  }
}

#' Determine maximum observed PK concentration
#'
#' @inheritParams assert_conc_time
#' @param check Run [assert_conc()]?
#' @return a number for the maximum concentration or NA if all concentrations
#'   are missing
#' @family NCA parameters for concentrations during the intervals
#' @export
pk.calc.cmax <- function(conc, check=TRUE) {
  if (check) {
    assert_conc(conc = conc)
  }
  if (length(conc) == 0 | all(is.na(conc))) {
    NA
  } else {
    max(conc, na.rm=TRUE)
  }
}
# Add the column to the interval specification
add.interval.col("cmax",
                 FUN="pk.calc.cmax",
                 values=c(FALSE, TRUE),
                 unit_type="conc",
                 pretty_name="Cmax",
                 desc="Maximum observed concentration",
                 depends=NULL)
PKNCA.set.summary(
  name="cmax",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' @describeIn pk.calc.cmax Determine the minimum observed PK
#'   concentration
#' @family NCA parameters for concentrations during the intervals
#' @examples
#' conc_data <- Theoph[Theoph$Subject == 1,]
#' pk.calc.cmin(conc_data$conc)
#' @export
pk.calc.cmin <- function(conc, check=TRUE) {
  if (check) {
    assert_conc(conc=conc)
  }
  if (length(conc) == 0 | all(is.na(conc))) {
    NA
  } else {
    min(conc, na.rm=TRUE)
  }
}
# Add the column to the interval specification
add.interval.col("cmin",
                 FUN="pk.calc.cmin",
                 values=c(FALSE, TRUE),
                 unit_type="conc",
                 pretty_name="Cmin",
                 desc="Minimum observed concentration",
                 depends=NULL)
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
#'   \item the `conc` and `time` must be the same length,
#'   \item the `time` may have no NAs,
#' }
#' `NA` will be returned if:
#' \enumerate{
#'   \item the length of `conc` and `time` is 0
#'   \item all `conc` is 0 or `NA`
#' }
#'
#' @inheritParams assert_conc_time
#' @inheritParams PKNCA.choose.option
#' @param first.tmax If there is more than time that matches the maximum
#'   concentration, should the first be considered as Tmax?  If not, then the
#'   last is considered Tmax.
#' @param check Run [assert_conc_time()]?
#' @returns The time of the maximum concentration
#' @examples
#' conc_data <- Theoph[Theoph$Subject == 1,]
#' pk.calc.tmax(conc = conc_data$conc, time = conc_data$Time)
#' @export
pk.calc.tmax <- function(conc, time,
                         options=list(),
                         first.tmax=NULL,
                         check=TRUE) {
  first.tmax <- PKNCA.choose.option(name="first.tmax", value=first.tmax, options=options)
  if (check) {
    assert_conc_time(conc = conc, time = time)
  }
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
# Add the column to the interval specification
add.interval.col("tmax",
                 FUN="pk.calc.tmax",
                 values=c(FALSE, TRUE),
                 unit_type="time",
                 pretty_name="Tmax",
                 desc="Time of the maximum observed concentration",
                 depends=NULL)
PKNCA.set.summary(
  name="tmax",
  description="median and range",
  point=business.median,
  spread=business.range
)

#' Determine time of last observed concentration above the limit of
#' quantification.
#'
#' `NA` will be returned if all `conc` are `NA` or 0.
#'
#' @inheritParams assert_conc_time
#' @param check Run [assert_conc_time()]?
#' @returns The time of the last observed concentration measurement
#' @export
pk.calc.tlast <- function(conc, time, check=TRUE) {
  if (check) {
    assert_conc_time(conc = conc, time = time)
  }
  if (all(conc %in% c(NA, 0))) {
    NA
  } else {
    max(time[!(conc %in% c(NA, 0))])
  }
}
# Add the column to the interval specification
add.interval.col("tlast",
                 FUN="pk.calc.tlast",
                 values=c(FALSE, TRUE),
                 unit_type="time",
                 pretty_name="Tlast",
                 desc="Time of the last concentration observed above the limit of quantification",
                 depends=NULL)
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
  if (check) {
    assert_conc_time(conc, time)
  }
  if (all(conc %in% c(NA, 0))) {
    NA
  } else {
    min(time[!(conc %in% c(NA, 0))])
  }
}
# Add the column to the interval specification
add.interval.col("tfirst",
                 FUN="pk.calc.tfirst",
                 values=c(FALSE, TRUE),
                 unit_type="time",
                 pretty_name="Tfirst",
                 desc="Time of the first concentration above the limit of quantification",
                 depends=NULL)
PKNCA.set.summary(
  name="tfirst",
  description="median and range",
  point=business.median,
  spread=business.range
)

#' Determine the last observed concentration above the limit of quantification
#' (LOQ).
#'
#' If all concentrations are missing, `NA_real_` is returned.  If all
#' concentrations are zero (below the limit of quantification) or missing, zero
#' is returned.  If Tlast is NA (due to no non-missing above LOQ measurements),
#' this will return `NA_real_`.
#'
#' @inheritParams assert_conc_time
#' @param check Run [assert_conc_time()]?
#' @returns The last observed concentration above the LOQ
#' @family NCA parameters for concentrations during the intervals
#' @export
pk.calc.clast.obs <- function(conc, time, check=TRUE) {
  if (check) {
    assert_conc_time(conc = conc, time = time)
  }
  if (all(is.na(conc))) {
    NA_real_
  } else if (all(conc %in% c(0, NA))) {
    0
  } else {
    tlast <- pk.calc.tlast(conc, time, check = FALSE)
    if (!is.na(tlast)) {
      conc[time %in% tlast]
    } else {
      NA_real_
    }
  }
}
# Add the column to the interval specification
add.interval.col("clast.obs",
                 FUN="pk.calc.clast.obs",
                 values=c(FALSE, TRUE),
                 unit_type="conc",
                 pretty_name="Clast",
                 desc="The last concentration observed above the limit of quantification",
                 depends=NULL)
PKNCA.set.summary(
  name="clast.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Calculate the effective half-life
#'
#' @details thalf.eff is `log(2)*mrt`.
#'
#' @param mrt the mean residence time to infinity
#' @return the numeric value of the effective half-life
#' @export
pk.calc.thalf.eff <- function(mrt) {
  log(2)*mrt
}
# Add the columns to the interval specification
add.interval.col("thalf.eff.obs",
                 FUN="pk.calc.thalf.eff",
                 values=c(FALSE, TRUE),
                 desc="The effective half-life (as determined from the MRTobs)",
                 unit_type="time",
                 pretty_name="Effective half-life (based on MRT,obs)",
                 formalsmap=list(mrt="mrt.obs"),
                 depends="mrt.obs")
PKNCA.set.summary(
  name="thalf.eff.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("thalf.eff.pred",
                 FUN="pk.calc.thalf.eff",
                 values=c(FALSE, TRUE),
                 unit_type="time",
                 pretty_name="Effective half-life (based on MRT,pred)",
                 desc="The effective half-life (as determined from the MRTpred)",
                 formalsmap=list(mrt="mrt.pred"),
                 depends="mrt.pred")
PKNCA.set.summary(
  name="thalf.eff.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("thalf.eff.last",
                 FUN="pk.calc.thalf.eff",
                 values=c(FALSE, TRUE),
                 unit_type="time",
                 pretty_name="Effective half-life (based on MRT,last)",
                 desc="The effective half-life (as determined from the MRTlast)",
                 formalsmap=list(mrt="mrt.last"),
                 depends="mrt.last")
PKNCA.set.summary(
  name="thalf.eff.last",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("thalf.eff.iv.obs",
                 FUN="pk.calc.thalf.eff",
                 values=c(FALSE, TRUE),
                 unit_type="time",
                 pretty_name="Effective half-life (for IV dosing, based on MRT,obs)",
                 desc="The effective half-life (as determined from the intravenous MRTobs)",
                 formalsmap=list(mrt="mrt.iv.obs"),
                 depends="mrt.iv.obs")
PKNCA.set.summary(
  name="thalf.eff.iv.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("thalf.eff.iv.pred",
                 FUN="pk.calc.thalf.eff",
                 values=c(FALSE, TRUE),
                 unit_type="time",
                 pretty_name="Effective half-life (for IV dosing, based on MRT,pred)",
                 desc="The effective half-life (as determined from the intravenous MRTpred)",
                 formalsmap=list(mrt="mrt.iv.pred"),
                 depends="mrt.iv.pred")
PKNCA.set.summary(
  name="thalf.eff.iv.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("thalf.eff.iv.last",
                 FUN="pk.calc.thalf.eff",
                 values=c(FALSE, TRUE),
                 unit_type="time",
                 pretty_name="Effective half-life (for IV dosing, based on MRTlast)",
                 desc="The effective half-life (as determined from the intravenous MRTlast)",
                 formalsmap=list(mrt="mrt.iv.last"),
                 depends="mrt.iv.last")
PKNCA.set.summary(
  name="thalf.eff.iv.last",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Calculate the AUC percent extrapolated
#'
#' @details aucpext is `100*(1-auclast/aucinf)`.
#'
#' @param auclast the area under the curve from time 0 to the last measurement
#'   above the limit of quantification
#' @param aucinf the area under the curve from time 0 to infinity
#' @returns The numeric value of the AUC percent extrapolated or `NA_real_` if
#'   any of the following are true `is.na(aucinf)`, `is.na(auclast)`,
#'   `aucinf <= 0`, or `auclast <= 0`.
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
  mask_calc <- !mask_na & !(aucinf %in% 0)
  if (any(mask_greater))
    rlang::warn(
      message = "aucpext is typically only calculated when aucinf is greater than auclast.",
      class = "pknca_aucpext_aucinf_le_auclast"
    )
  if (any(mask_negative))
    rlang::warn(
      message = "aucpext is typically only calculated when both aucinf and auclast are positive.",
      class = "pknca_aucpext_aucinf_auclast_positive"
    )
  ret[mask_calc] <-
    100*(1-auclast[mask_calc]/aucinf[mask_calc])
  ret
}

# Add the columns to the interval specification
add.interval.col("aucpext.obs",
                 FUN="pk.calc.aucpext",
                 values=c(FALSE, TRUE),
                 unit_type="%",
                 pretty_name="AUCpext (based on AUCinf,obs)",
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
                 unit_type="%",
                 pretty_name="AUCpext (based on AUCinf,pred)",
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
#' @param kel is `1/mrt`, not to be confused with lambda.z.
#'
#' @param mrt the mean residence time
#' @returns the numeric value of the elimination rate
#' @export
pk.calc.kel <- function(mrt) {
  1/mrt
}
# Add the columns to the interval specification
add.interval.col("kel.obs",
                 FUN="pk.calc.kel",
                 values=c(FALSE, TRUE),
                 unit_type="inverse_time",
                 pretty_name="Kel (based on AUCinf,obs)",
                 desc="Elimination rate (as calculated from the MRT with observed Clast)",
                 formalsmap=list(mrt="mrt.obs"),
                 depends="mrt.obs")
PKNCA.set.summary(
  name="kel.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("kel.pred",
                 FUN="pk.calc.kel",
                 values=c(FALSE, TRUE),
                 unit_type="inverse_time",
                 pretty_name="Kel (based on AUCinf,pred)",
                 desc="Elimination rate (as calculated from the MRT with predicted Clast)",
                 formalsmap=list(mrt="mrt.pred"),
                 depends="mrt.pred")
PKNCA.set.summary(
  name="kel.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("kel.last",
                 FUN="pk.calc.kel",
                 values=c(FALSE, TRUE),
                 unit_type="inverse_time",
                 pretty_name="Kel (based on AUClast)",
                 desc="Elimination rate (as calculated from the MRT using AUClast)",
                 formalsmap=list(mrt="mrt.last"),
                 depends="mrt.last")
PKNCA.set.summary(
  name="kel.last",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("kel.iv.obs",
                 FUN="pk.calc.kel",
                 values=c(FALSE, TRUE),
                 unit_type="inverse_time",
                 pretty_name="Kel (for IV dosing, based on AUCinf,obs)",
                 desc="Elimination rate (as calculated from the intravenous MRTobs)",
                 formalsmap=list(mrt="mrt.iv.obs"),
                 depends="mrt.iv.obs")
PKNCA.set.summary(
  name="kel.iv.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("kel.iv.pred",
                 FUN="pk.calc.kel",
                 values=c(FALSE, TRUE),
                 unit_type="inverse_time",
                 pretty_name="Kel (for IV dosing, based on AUCinf,pred)",
                 desc="Elimination rate (as calculated from the intravenous MRTpred)",
                 formalsmap=list(mrt="mrt.iv.pred"),
                 depends="mrt.iv.pred")
PKNCA.set.summary(
  name="kel.iv.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("kel.iv.last",
                 FUN="pk.calc.kel",
                 values=c(FALSE, TRUE),
                 unit_type="inverse_time",
                 pretty_name="Kel (for IV dosing, based on AUClast)",
                 desc="Elimination rate (as calculated from the intravenous MRTlast)",
                 formalsmap=list(mrt="mrt.iv.last"),
                 depends="mrt.iv.last")
PKNCA.set.summary(
  name="kel.iv.last",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Calculate the (observed oral) clearance
#'
#' @details cl is `dose/auc`.
#'
#' @param dose the dose administered
#' @param auc The area under the concentration-time curve.
#' @returns the numeric value of the total (CL) or observed oral clearance
#'   (CL/F)
#' @details If `dose` is the same length as the other inputs, then the output
#'   will be the same length as all of the inputs; the function assumes that you
#'   are calculating for multiple intervals simultaneously.  If the inputs other
#'   than `dose` are scalars and `dose` is a vector, then the function assumes
#'   multiple doses were given in a single interval, and the sum of the `dose`s
#'   will be used for the calculation.
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

# Add the columns to the interval specification
add.interval.col("cl.last",
                 FUN="pk.calc.cl",
                 values=c(FALSE, TRUE),
                 unit_type="clearance",
                 pretty_name="CL (based on AUClast)",
                 desc="Clearance or observed oral clearance calculated to Clast",
                 formalsmap=list(auc="auclast"),
                 depends="auclast")
PKNCA.set.summary(
  name="cl.last",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("cl.all",
                 FUN="pk.calc.cl",
                 values=c(FALSE, TRUE),
                 unit_type="clearance",
                 pretty_name="CL (based on AUCall)",
                 desc="Clearance or observed oral clearance calculated with AUCall",
                 formalsmap=list(auc="aucall"),
                 depends="aucall")
PKNCA.set.summary(
  name="cl.all",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("cl.obs",
                 FUN="pk.calc.cl",
                 values=c(FALSE, TRUE),
                 unit_type="clearance",
                 pretty_name="CL (based on AUCinf,obs)",
                 desc="Clearance or observed oral clearance calculated with observed Clast",
                 formalsmap=list(auc="aucinf.obs"),
                 depends="aucinf.obs")
PKNCA.set.summary(
  name="cl.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("cl.pred",
                 FUN="pk.calc.cl",
                 values=c(FALSE, TRUE),
                 unit_type="clearance",
                 pretty_name="CL (based on AUCinf,pred)",
                 desc="Clearance or observed oral clearance calculated with predicted Clast",
                 formalsmap=list(auc="aucinf.pred"),
                 depends="aucinf.pred")
PKNCA.set.summary(
  name="cl.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Calculate the absolute (or relative) bioavailability
#'
#' @details f is `(auc2/dose2)/(auc1/dose1)`.
#'
#' @param dose1 The dose administered in route or method 1
#' @param dose2 The dose administered in route or method 2
#' @param auc1 The AUC from 0 to infinity or 0 to tau administered in route or
#'   method 1
#' @param auc2 The AUC from 0 to infinity or 0 to tau administered in route or
#'   method 2
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
                 unit_type="fraction",
                 pretty_name="Bioavailability",
                 desc="Bioavailability or relative bioavailability",
                 depends=NULL)
PKNCA.set.summary(
  name="f",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Calculate the mean residence time (MRT) for single-dose data or linear
#' multiple-dose data.
#'
#' @details mrt is `aumc/auc - duration.dose/2` where `duration.dose =
#'   0` for oral administration.
#'
#' @param auc the AUC from 0 to infinity or 0 to tau
#' @param aumc the AUMC from 0 to infinity or 0 to tau
#' @param duration.dose The duration of the dose (usually an infusion duration
#'   for an IV infusion)
#' @returns the numeric value of the mean residence time
#' @seealso [pk.calc.mrt.md()]
#' @export
pk.calc.mrt <- function(auc, aumc) {
  pk.calc.mrt.iv(auc, aumc, duration.dose=0)
}
# Add the columns to the interval specification
add.interval.col("mrt.obs",
                 FUN="pk.calc.mrt",
                 values=c(FALSE, TRUE),
                 unit_type="time",
                 pretty_name="MRT (based on AUCinf,obs)",
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
                 unit_type="time",
                 pretty_name="MRT (based on AUCinf,pred)",
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
                 unit_type="time",
                 pretty_name="MRT (based on AUClast)",
                 desc="The mean residence time to the last observed concentration above the LOQ",
                 formalsmap=list(auc="auclast", aumc="aumclast"),
                 depends=c("auclast", "aumclast"))
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
# Add the columns to the interval specification
add.interval.col("mrt.iv.obs",
                 FUN="pk.calc.mrt.iv",
                 values=c(FALSE, TRUE),
                 unit_type="time",
                 pretty_name="MRT (for IV dosing, based on AUCinf,obs)",
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
                 unit_type="time",
                 pretty_name="MRT (for IV dosing, based on AUCinf,pred)",
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
                 unit_type="time",
                 pretty_name="MRT (for IV dosing, based on AUClast)",
                 desc="The mean residence time to the last observed concentration above the LOQ correcting for dosing duration",
                 formalsmap=list(auc="auclast", aumc="aumclast"),
                 depends=c("auclast", "aumclast"))
PKNCA.set.summary(
  name="mrt.iv.last",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Calculate the mean residence time (MRT) for multiple-dose data with nonlinear
#' kinetics.
#'
#' @details mrt.md is `aumctau/auctau + tau*(aucinf-auctau)/auctau` and should
#'   only be used for multiple dosing with equal intervals between doses.
#'
#' @param auctau the AUC from time 0 to the end of the dosing interval (tau).
#' @param aumctau the AUMC from time 0 to the end of the dosing interval (tau).
#' @param aucinf the AUC from time 0 to infinity (typically using single-dose
#'   data)
#' @inheritParams assert_dosetau
#' @details Note that if `aucinf == auctau` (as would be the assumption with
#'   linear kinetics), the equation becomes the same as the single-dose MRT.
#' @seealso [pk.calc.mrt()]
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
                 unit_type="time",
                 pretty_name="MRT (for multiple dosing, based on AUCinf,obs)",
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
                 unit_type="time",
                 pretty_name="MRT (for multiple dosing, based on AUCinf,pred)",
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
#' @details vz is `cl/lambda.z`.
#'
#' @inheritParams assert_lambdaz
#' @param cl the clearance (or apparent observed clearance)
#' @export
pk.calc.vz <- function(cl, lambda.z) {
  assert_lambdaz(lambda.z)
  # Ensure that cl is either a scalar or the same length as AUC
  # (more complex repeating patterns while valid for general R are
  # likely errors here).
  if (!(length(cl) %in% c(1, length(lambda.z))) |
      !(length(lambda.z) %in% c(1, length(cl))))
    stop("'cl' and 'lambda.z' must be the same length")
  cl/lambda.z
}
# Add the columns to the interval specification
add.interval.col("vz.obs",
                 FUN="pk.calc.vz",
                 values=c(FALSE, TRUE),
                 unit_type="volume",
                 pretty_name="Vz (based on AUCinf,obs)",
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
                 unit_type="volume",
                 pretty_name="Vz (based on AUCinf,pred)",
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
#' @details vss is `cl*mrt`.
#' @param cl the clearance
#' @param mrt the mean residence time
#' @return the volume of distribution at steady-state
#' @export
pk.calc.vss <- function(cl, mrt) {
  cl*mrt
}
# Add the columns to the interval specification
add.interval.col("vss.obs",
                 FUN="pk.calc.vss",
                 values=c(FALSE, TRUE),
                 unit_type="volume",
                 pretty_name="Vss (based on AUCinf,obs)",
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
                 unit_type="volume",
                 pretty_name="Vss (based on AUCinf,pred)",
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
                 unit_type="volume",
                 pretty_name="Vss (based on AUClast)",
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
                 unit_type="volume",
                 pretty_name="Vss (for IV dosing, based on AUCinf,obs)",
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
                 unit_type="volume",
                 pretty_name="Vss (for IV dosing, based on AUCinf,pred)",
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
                 unit_type="volume",
                 pretty_name="Vss (for IV dosing, based on AUClast)",
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
                 unit_type="volume",
                 pretty_name="Vss (for multiple-dose, based on Clast,obs)",
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
                 unit_type="volume",
                 pretty_name="Vss (for multiple-dose, based on Clast,pred)",
                 desc="The steady-state volume of distribution for nonlinear multiple-dose data using predicted Clast",
                 formalsmap=list(cl="cl.last", mrt="mrt.md.pred"),
                 depends=c("cl.last", "mrt.md.pred"))
PKNCA.set.summary(
  name="vss.md.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Calculate the average concentration during an interval.
#'
#' @details cav is `auc/(end-start)`.
#'
#' @param auc The area under the curve during the interval
#' @inheritParams assert_intervaltime_single
#' @returns The Cav (average concentration during the interval)
#' @export
pk.calc.cav <- function(auc, start, end) {
  ret <- auc/(end-start)
  mask_zero <- is.na(end) | is.na(start) | end == start
  if (any(mask_zero)) {
    ret[mask_zero] <- NA_real_
  }
  ret
}
add.interval.col(
  "cav",
  FUN = "pk.calc.cav",
  values = c(FALSE, TRUE),
  unit_type = "conc",
  pretty_name = "Cav",
  desc = "The average concentration during an interval (calculated with AUClast)",
  depends = "auclast",
  formalsmap = list(auc = "auclast")
)
add.interval.col(
  "cav.int.last",
  FUN = "pk.calc.cav",
  values = c(FALSE, TRUE),
  unit_type = "conc",
  pretty_name = "Cav",
  desc = "The average concentration during an interval (calculated with AUCint.last)",
  depends = "aucint.last",
  formalsmap = list(auc = "aucint.last"),
)
add.interval.col(
  "cav.int.all",
  FUN = "pk.calc.cav",
  values = c(FALSE, TRUE),
  unit_type = "conc",
  pretty_name = "Cav",
  desc = "The average concentration during an interval (calculated with AUCint.all)",
  depends = "aucint.all",
  formalsmap = list(auc = "aucint.all"),
)
add.interval.col(
  "cav.int.inf.obs",
  FUN = "pk.calc.cav",
  values = c(FALSE, TRUE),
  unit_type = "conc",
  pretty_name = "Cav",
  desc = "The average concentration during an interval (calculated with AUCint.inf.obs)",
  depends = "aucint.inf.obs",
  formalsmap = list(auc = "aucint.inf.obs"),
)
add.interval.col(
  "cav.int.inf.pred",
  FUN = "pk.calc.cav",
  values = c(FALSE, TRUE),
  unit_type = "conc",
  pretty_name = "Cav",
  desc = "The average concentration during an interval (calculated with AUCint.inf.pred)",
  depends = "aucint.inf.pred",
  formalsmap = list(auc = "aucint.inf.pred"),
)

PKNCA.set.summary(
  name = c("cav", "cav.int.last", "cav.int.all", "cav.int.inf.obs", "cav.int.inf.pred"),
  description = "geometric mean and geometric coefficient of variation",
  point = business.geomean,
  spread = business.geocv
)

#' Determine the trough (end of interval) concentration
#'
#' @inheritParams assert_conc_time
#' @inheritParams assert_intervaltime_single
#' @returns The concentration when `time == end`.  If none match, then `NA`
#' @family NCA parameters for concentrations during the intervals
#' @export
pk.calc.ctrough <- function(conc, time, end) {
  assert_conc_time(conc = conc, time = time)
  mask_end <- time %in% end
  if (sum(mask_end) == 1) {
    conc[mask_end]
  } else if (sum(mask_end) == 0) {
    NA_real_
  } else {
    # This should be impossible as assert_conc_time should catch
    # duplicates.
    stop("More than one time matches the starting time.  Please report this as a bug with a reproducible example.") # nocov
  }
}
add.interval.col("ctrough",
                 FUN="pk.calc.ctrough",
                 values=c(FALSE, TRUE),
                 unit_type="conc",
                 pretty_name="Ctrough",
                 desc="The trough (end of interval) concentration",
                 depends=NULL)
PKNCA.set.summary(
  name="ctrough",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Determine the concentration at the beginning of the interval
#'
#' @inheritParams assert_conc_time
#' @inheritParams assert_intervaltime_single
#' @return The concentration when `time == end`.  If none match, then `NA`
#' @family NCA parameters for concentrations during the intervals
#' @export
pk.calc.cstart <- function(conc, time, start) {
  assert_conc_time(conc = conc, time = time)
  mask_start <- time %in% start
  if (sum(mask_start) == 1) {
    conc[mask_start]
  } else if (sum(mask_start) == 0) {
    NA_real_
  } else {
    # This should be impossible as assert_conc_time should catch
    # duplicates.
    stop("More than one time matches the starting time.  Please report this as a bug with a reproducible example.") # nocov
  }
}
add.interval.col("cstart",
                 FUN="pk.calc.cstart",
                 values=c(FALSE, TRUE),
                 unit_type="conc",
                 pretty_name="Cstart",
                 desc="The predose concentration",
                 depends=NULL)
PKNCA.set.summary(
  name="cstart",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Determine the peak-to-trough ratio
#'
#' @details ptr is `cmax/ctrough`.
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
                 unit_type="fraction",
                 pretty_name="Peak-to-trough ratio",
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
#' @inheritParams assert_conc_time
#' @returns The time associated with the first increasing concentration
#' @export
pk.calc.tlag <- function(conc, time) {
  assert_conc_time(conc = conc, time = time)
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
                 unit_type="time",
                 pretty_name="Tlag",
                 desc="Lag time",
                 depends=NULL)
PKNCA.set.summary(
  name="tlag",
  description="median and range",
  point=business.median,
  spread=business.range
)

#' Determine the degree of fluctuation
#'
#' @details deg.fluc is `100*(cmax - cmin)/cav`.
#'
#' @param cmax The maximum observed concentration
#' @param cmin The minimum observed concentration
#' @param cav The average concentration in the interval
#' @returns The degree of fluctuation around the average concentration.
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
                 unit_type="%",
                 pretty_name="Degree of fluctuation",
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
#' @details swing is `100*(cmax - cmin)/cmin`.
#'
#' @param cmax The maximum observed concentration
#' @param cmin The minimum observed concentration
#' @returns The swing above the minimum concentration.  If `cmin` is zero, then
#'   the result is infinity.
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
                 unit_type="%",
                 pretty_name="Swing",
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
#' @inheritParams assert_conc_time
#' @param duration.dose The duration for the dosing administration (typically
#'   from IV infusion)
#' @param check Run [assert_conc_time()]?
#' @returns The concentration at the end of the infusion, `NA` if
#'   `duration.dose` is `NA`, or `NA` if all `time != duration.dose`
#' @export
pk.calc.ceoi <- function(conc, time, duration.dose=NA, check=TRUE) {
  if (check) {
    assert_conc_time(conc = conc, time = time)
  }
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
                 unit_type="conc",
                 pretty_name="Ceoi",
                 desc="Concentration at the end of infusion",
                 depends=NULL)
PKNCA.set.summary(
  name="ceoi",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Calculate the AUC above a given concentration
#'
#' Concentrations below the given concentration (`conc_above`) will be set
#' to zero.
#' @inheritParams pk.calc.time_above
#' @returns The AUC of the concentration above the limit
#' @export
pk.calc.aucabove <- function(conc, time, conc_above = NA_real_, ..., options=list()) {
  stopifnot(length(conc_above) == 1)
  stopifnot(is.numeric(conc_above))
  if (is.na(conc_above)) {
    ret <- structure(NA_real_, exclude = "Missing concentration to be above")
  } else {
    ret <-
      pk.calc.auc(
        conc=pmax(conc - conc_above, 0), time=time, ..., options=options,
        auc.type="AUCall",
        lambda.z=NA
      )
  }
  ret
}
add.interval.col(
  "aucabove.predose.all",
  FUN="pk.calc.aucabove",
  unit_type="auc",
  pretty_name="AUC,above",
  desc="The area under the concentration time the beginning of the interval to the last concentration above the limit of quantification plus the triangle from that last concentration to 0 at the first concentration below the limit of quantification, with a concentration subtracted from all concentrations and values below zero after subtraction set to zero",
  depends="cstart",
  formalsmap = list(conc_above = "cstart")
)
PKNCA.set.summary(
  name="aucabove.predose.all",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col(
  "aucabove.trough.all",
  FUN="pk.calc.aucabove",
  unit_type="auc",
  pretty_name="AUC,above",
  desc="The area under the concentration time the beginning of the interval to the last concentration above the limit of quantification plus the triangle from that last concentration to 0 at the first concentration below the limit of quantification, with a concentration subtracted from all concentrations and values below zero after subtraction set to zero",
  depends="ctrough",
  formalsmap = list(conc_above = "ctrough")
)
PKNCA.set.summary(
  name="aucabove.trough.all",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Count the number of concentration measurements in an interval
#'
#' `count_conc` and `count_conc_measured` are typically used for quality control
#' on the data to ensure that there are a sufficient number of non-missing
#' samples for a calculation and to ensure that data are consistent between
#' individuals.
#'
#' @inheritParams pk.calc.cmax
#' @returns a count of the non-missing concentrations (0 if all concentrations
#'   are missing)
#' @family NCA parameters for concentrations during the intervals
#' @export
pk.calc.count_conc <- function(conc, check=TRUE) {
  if (check) {
    assert_conc(conc)
  }
  sum(!is.na(conc))
}
# Add the column to the interval specification
add.interval.col(
  "count_conc",
  FUN = "pk.calc.count_conc",
  values = c(FALSE, TRUE),
  unit_type = "count",
  pretty_name = "Concentration count",
  desc = "Number of non-missing concentrations for an interval",
  depends = NULL
)

#' @describeIn pk.calc.count_conc Count the number of concentration measurements
#'   that are not missing, above, or below the limit of quantification in an
#'   interval
#'
#' @returns a count of the non-missing, measured (not below or above the limit
#'   of quantification) concentrations (0 if all concentrations are missing)
#' @export
pk.calc.count_conc_measured <- function(conc, check=TRUE) {
  if (check) {
    assert_conc(conc)
  }
  sum(!is.na(conc) & is.finite(conc) & conc > 0)
}
# Add the column to the interval specification
add.interval.col(
  "count_conc_measured",
  FUN="pk.calc.count_conc_measured",
  values=c(FALSE, TRUE),
  unit_type="count",
  pretty_name="Measured concentration count",
  desc="Number of measured and non BLQ/ALQ concentrations for an interval",
  depends=NULL
)
PKNCA.set.summary(
  name=c("count_conc", "count_conc_measured"),
  description="median and range",
  point=business.median,
  spread=business.range
)

#' Extract the dose used for calculations
#'
#' @inheritParams pk.calc.cl
#' @returns The total dose for an interval
#' @export
pk.calc.totdose <- function(dose) {
  sum(dose)
}
add.interval.col(
  "totdose",
  FUN="pk.calc.totdose",
  values=c(FALSE, TRUE),
  unit_type="dose",
  pretty_name="Total dose",
  desc="Total dose administered during an interval"
)
PKNCA.set.summary(
  name="totdose",
  description="median and range",
  point=business.median,
  spread=business.range
)

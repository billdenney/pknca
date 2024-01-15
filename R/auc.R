#' A compute the Area Under the (Moment) Curve
#'
#' Compute the area under the curve (AUC) and the area under the moment curve
#' (AUMC) for pharmacokinetic (PK) data.  AUC and AUMC are used for many
#' purposes when analyzing PK in drug development.
#'
#' `pk.calc.auc.last` is simply a shortcut setting the `interval`
#' parameter to `c(0, "last")`.
#'
#' Extrapolation beyond Clast occurs using the half-life and Clast,obs;
#' Clast,pred is not yet supported.
#'
#' If all conc input are zero, then the AU(M)C is zero.
#'
#' @inheritParams assert_conc_time
#' @inheritParams assert_intervaltime_single
#' @inheritParams choose_interval_method
#' @inheritParams assert_lambdaz
#' @inheritParams PKNCA.choose.option
#' @param clast,clast.obs,clast.pred The last concentration above the limit of
#'   quantification; this is used for AUCinf calculations.  If provided as
#'   clast.obs (observed clast value, default), AUCinf is AUCinf,obs. If
#'   provided as clast.pred, AUCinf is AUCinf,pred.
#' @param conc.blq How to handle BLQ values in between the first and last above
#'   LOQ concentrations. (See [clean.conc.blq()] for usage instructions.)
#' @param conc.na How to handle missing concentration values.  (See
#'   [clean.conc.na()] for usage instructions.)
#' @param check Run [assert_conc_time()], [clean.conc.blq()], and
#'   [clean.conc.na()]?
#' @param fun_linear The function to use for integration of the linear part of
#'   the curve (not required for AUC or AUMC functions)
#' @param fun_log The function to use for integration of the logarithmic part of
#'   the curve (if log integration is used; not required for AUC or AUMC
#'   functions)
#' @param fun_inf The function to use for extrapolation from the final
#'   measurement to infinite time (not required for AUC or AUMC functions.
#' @param ... For functions other than `pk.calc.auxc`, these values are passed
#'   to `pk.calc.auxc`
#' @returns A numeric value for the AU(M)C.
#' @aliases pk.calc.auc pk.calc.aumc pk.calc.auc.last
#' @family AUC calculations
#' @seealso [clean.conc.blq()]
#' @references
#'
#' Gabrielsson J, Weiner D.  "Section 2.8.1 Computation methods - Linear
#' trapezoidal rule."  Pharmacokinetic & Pharmacodynamic Data Analysis: Concepts
#' and Applications, 4th Edition.  Stockholm, Sweden: Swedish Pharmaceutical
#' Press, 2000.  162-4.
#'
#' Gabrielsson J, Weiner D.  "Section 2.8.3 Computation methods - Log-linear
#' trapezoidal rule."  Pharmacokinetic & Pharmacodynamic Data Analysis: Concepts
#' and Applications, 4th Edition.  Stockholm, Sweden: Swedish Pharmaceutical
#' Press, 2000.  164-7.
#' @examples
#' myconc <- c(0, 1, 2, 1, 0.5, 0.25, 0)
#' mytime <- c(0, 1, 2, 3, 4,   5,    6)
#' pk.calc.auc(myconc, mytime, interval=c(0, 6))
#' pk.calc.auc(myconc, mytime, interval=c(0, Inf))
#' @export
pk.calc.auxc <- function(conc, time, interval=c(0, Inf),
                         clast=pk.calc.clast.obs(conc, time, check=FALSE),
                         lambda.z=NA,
                         auc.type=c("AUClast", "AUCinf", "AUCall"),
                         options=list(),
                         method=NULL,
                         conc.blq=NULL,
                         conc.na=NULL,
                         check=TRUE,
                         fun_linear, fun_log, fun_inf) {
  # Check the inputs
  method <- PKNCA.choose.option(name="auc.method", value=method, options=options)
  conc.blq <- PKNCA.choose.option(name="conc.blq", value=conc.blq, options=options)
  conc.na <- PKNCA.choose.option(name="conc.na", value=conc.na, options=options)
  if (check) {
    assert_conc_time(conc = conc, time = time)
    data <-
      clean.conc.blq(
        conc, time,
        conc.blq=conc.blq,
        conc.na=conc.na,
        check=FALSE
      )
  } else {
    data <- data.frame(conc = conc, time = time)
  }
  if (nrow(data) == 0) {
    # All the data were missing
    return(
      structure(
        NA_real_,
        exclude="No data for AUC calculation"
      )
    )
  } else if (nrow(data) == 1) {
    return(
      structure(
        NA_real_,
        exclude="AUC cannot be calculated with only one measured concentration"
      )
    )
  } else if (all(data$conc %in% c(0, NA))) {
    # All the data were missing or 0 before excluding points
    return(structure(0, exclude="DO NOT EXCLUDE"))
  }
  auc.type <- match.arg(auc.type)
  interval <- assert_intervaltime_single(interval = interval)
  if (auc.type %in% "AUCinf" & is.finite(interval[2])) {
    warning("Requesting AUCinf when the end of the interval is not Inf")
  }

  # Subset the data to the range of interest ####
  interval_start <- interval[1]
  interval_end <- interval[2]
  # Find the first time point
  if (interval_start < min(data$time)) {
    warn_message <-
      sprintf(
        "Requesting an AUC range starting (%g) before the first measurement (%g) is not allowed",
        interval_start, min(data$time)
      )
    rlang::warn(message = warn_message, class = "pknca_warn_auc_before_first")
    return(structure(NA_real_, exclude=warn_message))
  } else if (interval_start > max(data$time)) {
    # Give this as a warning, but allow it to continue
    warning(sprintf("AUC start time (%g) is after the maximum observed time (%g)",
                    interval_start, max(data$time)))
  }
  # Ensure that we have clean concentration and time data.  This means that we
  # need to make sure that we have our starting point. Interpolation ensures
  # that (and will give the same answer if it already exists in the right form).
  conc_start <-
    interp.extrap.conc(
      data$conc, data$time,
      time.out = interval_start,
      lambda.z = lambda.z,
      method = method,
      auc.type = auc.type,
      check = FALSE
    )
  # Add that concentration and time to the vectors removing the
  # original version if it was there.
  data <- rbind(data.frame(conc=conc_start,
                           time=interval_start),
                data[!(data$time %in% interval_start), ])
  # * either have our ending point or the ending point is Inf
  if (is.finite(interval_end)) {
    conc_end <-
      interp.extrap.conc(
        data$conc, data$time,
        interval_end,
        method = method,
        auc.type = auc.type,
        check = FALSE
      )
    # !mask.end because we may be replacing an entry that is a 0.
    if (!is.na(conc_end)) {
      data <-
        rbind(
          data[!(data$time %in% interval_end), ],
          data.frame(conc=conc_end, time=interval_end)
        )
    }
  }
  # Subset the conc and time data to the interval of interest
  data <- data[interval_start <= data$time & data$time <= interval_end, , drop=FALSE]
  # Set the overall tlast
  tlast <- pk.calc.tlast(data$conc, data$time, check=FALSE)
  if (all(data$conc %in% 0)) {
    # All the data were missing or 0 after excluding points
    ret <- structure(0, exclude="DO NOT EXCLUDE")
  } else if (is.na(tlast)) {
    # All concentrations are BLQ (note that this has to be checked
    # after full subsetting and interpolation to ensure that it is
    # still true)
    stop("Unknown error with NA tlast but non-BLQ concentrations") # nocov
  } else {
    interval_method <- choose_interval_method(conc = data$conc, time = data$time, tlast = tlast, method = method, auc.type = auc.type, options = options)
    ret <-
      auc_integrate(
        conc = data$conc, time = data$time,
        clast = clast, tlast = tlast, lambda.z = lambda.z,
        interval_method = interval_method,
        fun_linear = fun_linear, fun_log = fun_log, fun_inf = fun_inf
      )
  }
  ret
}

#' @describeIn pk.calc.auxc Compute the area under the curve
#' @export
pk.calc.auc <- function(conc, time, ..., options=list()) {
  pk.calc.auxc(
    conc=conc, time=time, ...,
    options=options,
    fun_linear=aucintegrate_linear,
    fun_log=aucintegrate_log,
    fun_inf=aucintegrate_inf
  )
}

# Note that lambda.z is set to NA for both auc.last and auc.all because all
# interpolation should happen within given points. lambda.z should not be used,
# and if it is used, that should be caught as an error.
#' @describeIn pk.calc.auxc Compute the AUClast.
#' @export
pk.calc.auc.last <- function(conc, time, ..., options=list()) {
  if ("auc.type" %in% names(list(...)))
    stop("auc.type cannot be changed when calling pk.calc.auc.last, please use pk.calc.auc")
  pk.calc.auc(conc=conc, time=time, ...,
              options=options,
              auc.type="AUClast",
              lambda.z=NA)
}

#' @describeIn pk.calc.auxc Compute the AUCinf
#' @export
pk.calc.auc.inf <- function(conc, time, ..., options=list(), lambda.z) {
  if ("auc.type" %in% names(list(...)))
    stop("auc.type cannot be changed when calling pk.calc.auc.inf, please use pk.calc.auc")
  pk.calc.auc(conc=conc, time=time, ...,
              options=options,
              auc.type="AUCinf",
              lambda.z=lambda.z)
}

#' @describeIn pk.calc.auxc Compute the AUCinf with the observed Clast.
#' @export
pk.calc.auc.inf.obs <- function(conc, time, clast.obs, ..., options=list(),
                                lambda.z) {
  pk.calc.auc.inf(conc=conc, time=time, clast=clast.obs, ...,
                  options=options,
                  lambda.z=lambda.z)
}

#' @describeIn pk.calc.auxc Compute the AUCinf with the predicted Clast.
#' @export
pk.calc.auc.inf.pred <- function(conc, time, clast.pred, ..., options=list(),
                                 lambda.z) {
  pk.calc.auc.inf(conc=conc, time=time, clast=clast.pred, ...,
                  options=options,
                  lambda.z=lambda.z)
}

#' @describeIn pk.calc.auxc Compute the AUCall.
#' @export
pk.calc.auc.all <- function(conc, time, ..., options=list()) {
  if ("auc.type" %in% names(list(...)))
    stop("auc.type cannot be changed when calling pk.calc.auc.all, please use pk.calc.auc")
  pk.calc.auc(conc=conc, time=time, ..., options=options,
              auc.type="AUCall",
              lambda.z=NA)
}


#' @describeIn pk.calc.auxc Compute the area under the moment curve
#' @export
pk.calc.aumc <- function(conc, time, ..., options=list()) {
  pk.calc.auxc(conc=conc, time=time, ..., options=options,
    fun_linear = aumcintegrate_linear,
    fun_log = aumcintegrate_log,
    fun_inf = aumcintegrate_inf
  )
}

#' @describeIn pk.calc.auxc Compute the AUMClast.
#' @export
pk.calc.aumc.last <- function(conc, time, ..., options=list()) {
  if ("auc.type" %in% names(list(...)))
    stop("auc.type cannot be changed when calling pk.calc.aumc.last, please use pk.calc.aumc")
  pk.calc.aumc(conc=conc, time=time, ..., options=options,
               auc.type="AUClast",
               lambda.z=NA)
}

#' @describeIn pk.calc.auxc Compute the AUMCinf
#' @export
pk.calc.aumc.inf <- function(conc, time, ..., options=list(),
                             lambda.z) {
  if ("auc.type" %in% names(list(...))) {
    stop("auc.type cannot be changed when calling pk.calc.aumc.inf, please use pk.calc.aumc")
  }
  pk.calc.aumc(conc=conc, time=time, ..., options=options,
               auc.type="AUCinf",
               lambda.z=lambda.z)
}

#' @describeIn pk.calc.auxc Compute the AUMCinf with the observed Clast.
#' @export
pk.calc.aumc.inf.obs <- function(conc, time, clast.obs, ..., options=list(),
                                 lambda.z) {
  pk.calc.aumc.inf(conc=conc, time=time, clast=clast.obs, ..., options=options,
                   lambda.z=lambda.z)
}

#' @describeIn pk.calc.auxc Compute the AUMCinf with the predicted Clast.
#' @export
pk.calc.aumc.inf.pred <- function(conc, time, clast.pred, ..., options=list(),
                             lambda.z) {
  pk.calc.aumc.inf(conc=conc, time=time, clast=clast.pred, ..., options=options,
                   lambda.z=lambda.z)
}

#' @describeIn pk.calc.auxc Compute the AUMCall.
#' @export
pk.calc.aumc.all <- function(conc, time, ..., options=list()) {
  if ("auc.type" %in% names(list(...)))
    stop("auc.type cannot be changed when calling pk.calc.aumc.all, please use pk.calc.aumc")
  pk.calc.aumc(conc=conc, time=time, ..., options=options,
               auc.type="AUCall",
               lambda.z=NA)
}

# Add the columns to the interval specification
add.interval.col("aucinf.obs",
                 FUN="pk.calc.auc.inf.obs",
                 values=c(FALSE, TRUE),
                 unit_type="auc",
                 pretty_name="AUCinf,obs",
                 desc="The area under the concentration time curve from the beginning of the interval to infinity with extrapolation to infinity from the observed Clast",
                 depends=c("lambda.z", "clast.obs"))

add.interval.col("aucinf.pred",
                 FUN="pk.calc.auc.inf.pred",
                 values=c(FALSE, TRUE),
                 unit_type="auc",
                 pretty_name="AUCinf,pred",
                 desc="The area under the concentration time curve from the beginning of the interval to infinity with extrapolation to infinity from the predicted Clast",
                 depends=c("lambda.z", "clast.pred"))

add.interval.col("auclast",
                 FUN="pk.calc.auc.last",
                 values=c(FALSE, TRUE),
                 unit_type="auc",
                 pretty_name="AUClast",
                 desc="The area under the concentration time curve from the beginning of the interval to the last concentration above the limit of quantification")

add.interval.col("aucall",
                 FUN="pk.calc.auc.all",
                 values=c(FALSE, TRUE),
                 unit_type="auc",
                 pretty_name="AUCall",
                 desc="The area under the concentration time curve from the beginning of the interval to the last concentration above the limit of quantification plus the triangle from that last concentration to 0 at the first concentration below the limit of quantification"
)

add.interval.col("aumcinf.obs",
                 FUN="pk.calc.aumc.inf.obs",
                 values=c(FALSE, TRUE),
                 unit_type="aumc",
                 pretty_name="AUMC,inf,obs",
                 desc="The area under the concentration time moment curve from the beginning of the interval to infinity with extrapolation to infinity from the observed Clast",
                 depends=c("lambda.z", "clast.obs"))

add.interval.col("aumcinf.pred",
                 FUN="pk.calc.aumc.inf.pred",
                 values=c(FALSE, TRUE),
                 unit_type="aumc",
                 pretty_name="AUMC,inf,pred",
                 desc="The area under the concentration time moment curve from the beginning of the interval to infinity with extrapolation to infinity from the predicted Clast",
                 depends=c("lambda.z", "clast.pred"))

add.interval.col("aumclast",
                 FUN="pk.calc.aumc.last",
                 values=c(FALSE, TRUE),
                 unit_type="aumc",
                 pretty_name="AUMC,last",
                 desc="The area under the concentration time moment curve from the beginning of the interval to the last concentration above the limit of quantification")

add.interval.col("aumcall",
                 FUN="pk.calc.aumc.all",
                 values=c(FALSE, TRUE),
                 unit_type="aumc",
                 pretty_name="AUMC,all",
                 desc="The area under the concentration time moment curve from the beginning of the interval to the last concentration above the limit of quantification plus the moment of the triangle from that last concentration to 0 at the first concentration below the limit of quantification")

PKNCA.set.summary(
  name=
    c(
      "aucinf.obs", "aucinf.pred", "auclast", "aucall",
      "aumcinf.obs", "aumcinf.pred", "aumclast", "aumcall"
    ),
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

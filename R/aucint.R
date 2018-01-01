#' Calculate the AUC over an interval with interpolation and/or
#' extrapolation of concentrations for the beginning and end of the
#' interval.
#'
#' @param conc Concentration measured
#' @param time Time of concentration measurement (must be monotonically
#'   increasing and the same length as the concentration data)
#' @param interval Numeric vector of two numbers for the start and end
#'   time of integration
#' @param start,end The start and end of the interval (cannot be given
#'   if \code{interval} is given)
#' @param clast,clast.obs,clast.pred The last concentration above the
#'   limit of quantification; this is used for AUCinf calculations.  If
#'   provided as clast.obs (observed clast value, default), AUCinf is
#'   AUCinf,obs. If provided as clast.pred, AUCinf is AUCinf,pred.
#' @param lambda.z The elimination rate (in units of inverse time) for
#'   extrapolation
#' @param time.dose,route,duration.dose The time of doses, route of
#'   administration, and duration of dose used with interpolation and
#'   extrapolation of concentration data (see
#'   \code{\link{interp.extrap.conc.dose}}).  If \code{NULL},
#'   \code{\link{interp.extrap.conc}} will be used instead (assuming
#'   that no doses affecting concentrations are in the interval).
#' @param method The method for integration (either 'lin up/log down' or
#'   'linear')
#' @param auc.type The type of AUC to compute.  Choices are 'AUCinf',
#'   'AUClast', and 'AUCall'.
#' @param ... Additional arguments passed to \code{pk.calc.auxc} and
#'   \code{interp.extrap.conc}
#' @param options List of changes to the default
#'   \code{\link{PKNCA.options}} for calculations.
#' @seealso \code{\link{pk.calc.auxc}}, \code{\link{PKNCA.options}},
#'   \code{\link{interp.extrap.conc.dose}}
#' @export
pk.calc.aucint <- function(conc, time, interval, start, end,
                           clast=pk.calc.clast.obs(conc, time), lambda.z=NA,
                           time.dose=NULL, route="extravascular", duration.dose=0,
                           method=PKNCA.choose.option("auc.method", options),
                           auc.type="AUClast", ...,
                           options=list()) {
  if (missing(interval)) {
    if (missing(start) | missing(end)) {
      stop("If interval is not given, start and end must be given.")
    } else if (length(start) != 1) {
      stop("start must be a scalar")
    } else if (length(end) != 1) {
      stop("end must be a scalar")
    } else if (!is.numeric(start) | is.factor(start) | is.infinite(start)) {
      stop("start must be a finite number")
    } else if (!is.numeric(end) | is.factor(end) | is.infinite(end)) {
      stop("end must be a finite number")
    }
    interval <- c(start, end)
  } else if (!missing(start) | !missing(end)) {
    stop("start and end cannot be given if interval is given")
  } else if (length(interval) != 2) {
    stop("interval must be a vector with 2 elements")
  } else if (!is.numeric(interval) | is.factor(interval) | any(is.infinite(interval))) {
    stop("interval must be numeric and finite")
  }
  if (interval[1] >= interval[2]) {
    stop("interval start must be before interval end.")
  }
  missing_times <- setdiff(c(interval, time.dose), time)
  # Handle the potential double-calculation (before/after tlast) with AUCinf
  conc_clast <- NULL
  time_clast <- NULL
  if (auc.type %in% "AUCinf") {
    tlast <- pk.calc.tlast(conc=conc, time=time)
    if (clast != pk.calc.clast.obs(conc=conc, time=time) &
        interval[2] > tlast) {
      # If using clast.pred, we need to doubly calculate at tlast.
      conc_clast <- clast
      time_clast <- tlast
    }
  }
  if (length(missing_times)) {
    if (is.null(time.dose)) {
      missing_conc <-
        interp.extrap.conc(
          conc=conc, time=time,
          time.out=missing_times,
          interp.method=method,
          extrap.method=auc.type,
          clast=clast, lambda.z=lambda.z,
          options=options,
          ...)
    } else {
      missing_conc <-
        interp.extrap.conc.dose(
          conc=conc, time=time,
          time.out=missing_times,
          interp.method=method,
          extrap.method=auc.type,
          clast=clast, lambda.z=lambda.z,
          options=options,
          # arguments specific to interp.extrap.conc.dose
          time.dose=time.dose,
          route.dose=route,
          duration.dose=duration.dose,
          out.after=FALSE,
          ...)
    }
    new_data <- data.frame(conc=c(conc, conc_clast, missing_conc),
                           time=c(time, time_clast, missing_times))
    new_data <- new_data[new_data$time >= interval[1] &
                           new_data$time <= interval[2],]
    new_data <- new_data[order(new_data$time),]
    conc_interp <- new_data$conc
    time_interp <- new_data$time
    if (any(mask_na_conc <- is.na(conc_interp))) {
      missing_times <- time_interp[mask_na_conc]
      warning_message <-
        if (any(is.na(lambda.z))) {
          paste("Some interpolated/extrapolated concentration values are missing",
                "(may be due to interpolating or extrapolating over a dose with lambda.z=NA).",
                "Time points with missing data are: ",
                paste(missing_times, collapse=", "))
        } else {
          paste("Some interpolated/extrapolated concentration values are missing",
                "Time points with missing data are: ",
                paste(missing_times, collapse=", "))
        }
      warning(warning_message)
      return(NA_real_)
    }
  } else {
    conc_interp <- conc
    time_interp <- time
  }
  # AUCinf traces an AUClast curve (because the interval doesn't go to
  # infinity) while AUCall and AUClast trace their own curves.
  auc.type_map <-
    list(AUCinf="AUClast",
         AUCall="AUCall",
         AUCinf="AUClast")[[auc.type]]
  pk.calc.auc(conc=conc_interp, time=time_interp,
              interval=interval,
              clast=clast, lambda.z=lambda.z,
              auc.type=auc.type_map,
              options=options,
              method=method,
              ...,
              check=FALSE)
}

#' @describeIn pk.calc.aucint Interpolate or extrapolate concentrations
#'   for AUClast
#' @export
pk.calc.aucint.last <- function(conc, time, start, end, ..., options=list()) {
  pk.calc.aucint(conc=conc, time=time, interval=c(start, end), options=options, ...,
                 auc.type="AUClast")
}
#' @describeIn pk.calc.aucint Interpolate or extrapolate concentrations
#'   for AUCall
#' @export
pk.calc.aucint.all <- function(conc, time, start, end, ..., options=list()) {
  pk.calc.aucint(conc=conc, time=time, interval=c(start, end), options=options, ...,
                 auc.type="AUCall")
}
#' @describeIn pk.calc.aucint Interpolate or extrapolate concentrations
#'   for AUCinf.obs
#' @export
pk.calc.aucint.inf.obs <- function(conc, time, start, end, lambda.z, clast.obs, ..., options=list()) {
  pk.calc.aucint(conc=conc, time=time, interval=c(start, end),
                 lambda.z=lambda.z, clast=clast.obs,
                 options=options, ...,
                 auc.type="AUCinf")
}
#' @describeIn pk.calc.aucint Interpolate or extrapolate concentrations
#'   for AUCinf.pred
#' @export
pk.calc.aucint.inf.pred <- function(conc, time, start, end, lambda.z, clast.pred, ..., options=list()) {
  pk.calc.aucint(conc=conc, time=time, interval=c(start, end),
                 lambda.z=lambda.z, clast=clast.pred,
                 options=options, ...,
                 auc.type="AUCinf")
}

add.interval.col("aucint.last",
                 FUN="pk.calc.aucint.last",
                 values=c(FALSE, TRUE),
                 desc="The area under the concentration time curve in the interval extrapolating from Tlast to infinity with zeros (matching AUClast)",
                 formalsmap=list(conc="conc.group", time="time.group"),
                 depends=c())
PKNCA.set.summary("aucint.last", business.geomean, business.geocv)

add.interval.col("aucint.all",
                 FUN="pk.calc.aucint.all",
                 values=c(FALSE, TRUE),
                 desc="The area under the concentration time curve in the interval extrapolating from Tlast to infinity with the triangle from Tlast to the next point and zero thereafter (matching AUCall)",
                 formalsmap=list(conc="conc.group", time="time.group"),
                 depends=c())
PKNCA.set.summary("aucint.all", business.geomean, business.geocv)

add.interval.col("aucint.inf.obs",
                 FUN="pk.calc.aucint.inf.obs",
                 values=c(FALSE, TRUE),
                 desc="The area under the concentration time curve in the interval extrapolating from Tlast to infinity with zeros (matching AUClast)",
                 formalsmap=list(conc="conc.group", time="time.group"),
                 depends=c("lambda.z", "clast.obs"))
PKNCA.set.summary("aucint.last", business.geomean, business.geocv)

add.interval.col("aucint.inf.pred",
                 FUN="pk.calc.aucint.inf.pred",
                 values=c(FALSE, TRUE),
                 desc="The area under the concentration time curve in the interval extrapolating from Tlast to infinity with the triangle from Tlast to the next point and zero thereafter (matching AUCall)",
                 formalsmap=list(conc="conc.group", time="time.group"),
                 depends=c("lambda.z", "clast.pred"))
PKNCA.set.summary("aucint.all", business.geomean, business.geocv)

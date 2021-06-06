#' A compute the Area Under the (Moment) Curve
#' 
#' Compute the area under the curve (AUC) and the area under the moment curve
#' (AUMC) for pharmacokinetic (PK) data.  AUC and AUMC are used for many
#' purposes when analyzing PK in drug development.
#' 
#' \code{pk.calc.auc.last} is simply a shortcut setting the \code{interval}
#' parameter to \code{c(0, "last")}.
#' 
#' Extrapolation beyond Clast occurs using the half-life and Clast,obs;
#' Clast,pred is not yet supported.
#' 
#' If all conc input are zero, then the AU(M)C is zero.
#' 
#' @param conc Concentration measured
#' @param time Time of concentration measurement (must be monotonically
#'   increasing and the same length as the concentration data)
#' @param interval Numeric vector of two numbers for the start and end time of
#'   integration
#' @param auc.type The type of AUC to compute.  Choices are 'AUCinf', 'AUClast',
#'   and 'AUCall'.
#' @param clast,clast.obs,clast.pred The last concentration above the limit of 
#'   quantification; this is used for AUCinf calculations.  If provided as
#'   clast.obs (observed clast value, default), AUCinf is AUCinf,obs. If
#'   provided as clast.pred, AUCinf is AUCinf,pred.
#' @param lambda.z The elimination rate (in units of inverse time) for 
#'   extrapolation
#' @param options List of changes to the default \code{\link{PKNCA.options}} for
#'   calculations.
#' @param method The method for integration (either 'lin up/log down' or
#'   'linear')
#' @param conc.blq How to handle BLQ values in between the first and last above
#'   LOQ concentrations. (See \code{\link{clean.conc.blq}} for usage
#'   instructions.)
#' @param conc.na How to handle missing concentration values.  (See 
#'   \code{\link{clean.conc.na}} for usage instructions.)
#' @param check Run \code{\link{check.conc.time}}, \code{\link{clean.conc.blq}},
#'   and \code{\link{clean.conc.na}}?
#' @param fun.linear The function to use for integration of the linear part of
#'   the curve (not required for AUC or AUMC functions)
#' @param fun.log The function to use for integration of the logarithmic part of
#'   the curve (if log integration is used; not required for AUC or AUMC
#'   functions)
#' @param fun.inf The function to use for extrapolation from the final 
#'   measurement to infinite time (not required for AUC or AUMC functions.
#' @param ... For functions other than \code{pk.calc.auxc}, these values are
#'   passed to \code{pk.calc.auxc}
#' @return A numeric value for the AU(M)C.
#' @aliases pk.calc.auc pk.calc.aumc pk.calc.auc.last
#' @family AUC calculations
#' @seealso \code{\link{clean.conc.blq}}
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
                         fun.linear, fun.log, fun.inf) {
  # Check the inputs
  method <- PKNCA.choose.option(name="auc.method", value=method, options=options)
  conc.blq <- PKNCA.choose.option(name="conc.blq", value=conc.blq, options=options)
  conc.na <- PKNCA.choose.option(name="conc.na", value=conc.na, options=options)
  if (check) {
    check.conc.time(conc, time)
    data <-
      clean.conc.blq(
        conc, time,
        conc.blq=conc.blq,
        conc.na=conc.na,
        check=FALSE
      )
  } else {
    data <- data.frame(conc, time)
  }
  if (nrow(data) == 0) {
    # All the data were missing
    return(NA)
  } else if (all(data$conc %in% c(0, NA))) {
    # All the data were missing or 0 before excluding points
    return(structure(0, exclude="DO NOT EXCLUDE"))
  }
  auc.type <- match.arg(auc.type)
  if (interval[1] >= interval[2])
    stop("The AUC interval must be increasing")
  if (auc.type %in% "AUCinf" &
        is.finite(interval[2]))
    warning("Requesting AUCinf when the end of the interval is not Inf")
  ##############################
  # Subset the data to the range of interest
  interval_start <- interval[1]
  interval_end <- interval[2]
  # Find the first time point
  if (interval_start < min(data$time)) {
    warning(sprintf(
      "Requesting an AUC range starting (%g) before the first measurement (%g) is not allowed",
      interval_start, min(data$time)))
    return(NA)
  } else if (interval_start > max(data$time)) {
    # Give this as a warning, but allow it to continue
    warning(sprintf("AUC start time (%g) is after the maximum observed time (%g)",
                    interval_start, max(data$time)))
  }
  # Ensure that we have clean concentration and time data.  This means that we
  # need to make sure that we have our starting point. Interpolation ensures
  # that (and will give the same answer if it already exists in the right form).
  conc_start <-
    interp.extrap.conc(data$conc, data$time,
                       time.out=interval_start,
                       lambda.z=lambda.z,
                       interp.method=method,
                       extrap.method=auc.type,
                       check=FALSE)
  # Add that concentration and time to the vectors removing the
  # original version if it was there.
  data <- rbind(data.frame(conc=conc_start,
                           time=interval_start),
                data[!(data$time %in% interval_start),])
  # * either have our ending point or the ending point is Inf
  if (is.finite(interval_end)) {
    conc_end <-
      interp.extrap.conc(data$conc, data$time,
                         interval_end,
                         interp.method=method,
                         extrap.method=auc.type,
                         check=FALSE)
    # !mask.end because we may be replacing an entry that is a 0.
    data <- rbind(data[!(data$time %in% interval_end),],
                  data.frame(conc=conc_end,
                             time=interval_end))
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
    # ############################
    # Compute the AUxC
    # Compute it in linear space from the start to Tlast
    if (auc.type %in% "AUCall" &
        tlast != max(data$time)) {
      # Include the first point after tlast if it exists and we are
      # looking for AUCall
      idx_1 <- seq_len(sum(data$time <= tlast))
    } else {
      idx_1 <- seq_len(sum(data$time <= tlast) - 1)
    }
    idx_2 <- idx_1 + 1
    if (method %in% "linear") {
      ret <- fun.linear(data$conc[idx_1], data$conc[idx_2],
                        data$time[idx_1], data$time[idx_2])
    } else if (method %in% "lin up/log down") {
      # Compute log down if required (and if the later point is not 0)
      mask_down <- (data$conc[idx_2] < data$conc[idx_1] &
                      data$conc[idx_2] != 0)
      mask_up <- !mask_down
      idx_1_down <- idx_1[mask_down]
      idx_2_down <- idx_2[mask_down]
      idx_1_up <- idx_1[mask_up]
      idx_2_up <- idx_2[mask_up]
      ret <- rep(NA_real_, length(idx_1))
      ret[mask_up] <-
        fun.linear(data$conc[idx_1_up], data$conc[idx_2_up],
                   data$time[idx_1_up], data$time[idx_2_up])
      ret[mask_down] <-
        fun.log(data$conc[idx_1_down], data$conc[idx_2_down],
                data$time[idx_1_down], data$time[idx_2_down])
    } else if (!(method %in% "linear")) {
      # This should have already been caught, but the test exists to double-check
      stop("Invalid AUC integration method") # nocov
    }
    if (auc.type %in% "AUCinf") {
      # Whether AUCinf,obs or AUCinf,pred is calculated depends on if clast,obs
      # or clast,pred is passed in.
      ret[length(ret)+1] <- fun.inf(clast, tlast, lambda.z)
    }
    ret <- sum(ret)
  }
  ret
}

fun.auc.linear <- function(conc.1, conc.2, time.1, time.2)
  (time.2-time.1)*(conc.2+conc.1)/2
fun.auc.log <- function(conc.1, conc.2, time.1, time.2)
  (time.2-time.1)*(conc.2-conc.1)/log(conc.2/conc.1)
fun.auc.inf <- function(clast, tlast, lambda.z)
  clast/lambda.z

#' @describeIn pk.calc.auxc Compute the area under the curve
#' @export
pk.calc.auc <- function(conc, time, ..., options=list())
  pk.calc.auxc(conc=conc, time=time, ...,
               options=options,
               fun.linear=fun.auc.linear,
               fun.log=fun.auc.log,
               fun.inf=fun.auc.inf)

## Note that lambda.z is set to NA for both auc.last and auc.all
## because all interpolation should happen within given points.
## lambda.z should not be used, and if it is used, that should be
## caught as an error.
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
pk.calc.auc.inf <- function(conc, time, ..., options=list(),
                            lambda.z) {
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
pk.calc.aumc <- function(conc, time, ..., options=list())
  pk.calc.auxc(conc=conc, time=time, ..., options=options,
    fun.linear=function(conc.1, conc.2, time.1, time.2) {
      (time.2-time.1)*(conc.2*time.2+conc.1*time.1)/2
    },
    fun.log=function(conc.1, conc.2, time.1, time.2) {
      ((time.2-time.1)*(conc.2*time.2-conc.1*time.1)/log(conc.2/conc.1)-
       (time.2-time.1)^2*(conc.2-conc.1)/(log(conc.2/conc.1)^2))
    },
    fun.inf=function(conc.last, time.last, lambda.z) {
      (conc.last*time.last/lambda.z) + conc.last/(lambda.z^2)
    })

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
  if ("auc.type" %in% names(list(...)))
    stop("auc.type cannot be changed when calling pk.calc.aumc.inf, please use pk.calc.aumc")
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

## Add the columns to the interval specification
add.interval.col("aucinf.obs",
                 FUN="pk.calc.auc.inf.obs",
                 values=c(FALSE, TRUE),
                 desc="The area under the concentration time curve from the beginning of the interval to infinity with extrapolation to infinity from the observed Clast",
                 depends=c("lambda.z", "clast.obs"))
PKNCA.set.summary(
  name="aucinf.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

add.interval.col("aucinf.pred",
                 FUN="pk.calc.auc.inf.pred",
                 values=c(FALSE, TRUE),
                 desc="The area under the concentration time curve from the beginning of the interval to infinity with extrapolation to infinity from the predicted Clast",
                 depends=c("lambda.z", "clast.pred"))
PKNCA.set.summary(
  name="aucinf.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

add.interval.col("auclast",
                 FUN="pk.calc.auc.last",
                 values=c(FALSE, TRUE),
                 desc="The area under the concentration time curve from the beginning of the interval to the last concentration above the limit of quantification",
                 depends=c())
PKNCA.set.summary(
  name="auclast",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

add.interval.col("aucall",
                 FUN="pk.calc.auc.all",
                 values=c(FALSE, TRUE),
                 desc="The area under the concentration time curve from the beginning of the interval to the last concentration above the limit of quantification plus the triangle from that last concentration to 0 at the first concentration below the limit of quantification",
                 depends=c())
PKNCA.set.summary(
  name="aucall",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

add.interval.col("aumcinf.obs",
                 FUN="pk.calc.aumc.inf.obs",
                 values=c(FALSE, TRUE),
                 desc="The area under the concentration time moment curve from the beginning of the interval to infinity with extrapolation to infinity from the observed Clast",
                 depends=c("lambda.z", "clast.obs"))
PKNCA.set.summary(
  name="aumcinf.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

add.interval.col("aumcinf.pred",
                 FUN="pk.calc.aumc.inf.pred",
                 values=c(FALSE, TRUE),
                 desc="The area under the concentration time moment curve from the beginning of the interval to infinity with extrapolation to infinity from the predicted Clast",
                 depends=c("lambda.z", "clast.pred"))
PKNCA.set.summary(
  name="aumcinf.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

add.interval.col("aumclast",
                 FUN="pk.calc.aumc.last",
                 values=c(FALSE, TRUE),
                 desc="The area under the concentration time moment curve from the beginning of the interval to the last concentration above the limit of quantification",
                 depends=c())
PKNCA.set.summary(
  name="aumclast",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

add.interval.col("aumcall",
                 FUN="pk.calc.aumc.all",
                 values=c(FALSE, TRUE),
                 desc="The area under the concentration time moment curve from the beginning of the interval to the last concentration above the limit of quantification plus the moment of the triangle from that last concentration to 0 at the first concentration below the limit of quantification",
                 depends=c())
PKNCA.set.summary(
  name="aumcall",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

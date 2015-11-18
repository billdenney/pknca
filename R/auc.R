#' A compute the Area Under the (Moment) Curve
#'
#' Compute the area under the curve (AUC) and the area under the
#' moment curve (AUMC) for pharmacokinetic (PK) data.  AUC and AUMC
#' are used for many purposes when analyzing PK in drug development.
#'
#' \code{pk.calc.auc.last} is simply a shortcut setting the
#' \code{interval} parameter to \code{c(0, "last")}.
#'
#' Extrapolation beyond Clast occurs using the half-life and
#' Clast,obs; Clast,pred is not yet supported.
#'
#' @param conc Concentration measured
#' @param time Time of concentration measurement (must be
#' monotonically increasing and the same length as the concentration
#' data)
#' @param interval Numeric vector of two numbers for the start and
#' end time of integration
#' @param auc.type The type of AUC to compute.  Choices are 'AUCinf',
#' 'AUClast', and 'AUCall'.
#' @param clast The last concentration above the limit of
#' quantification; this is used for AUCinf calculations.  If provided
#' as clast.obs (observed clast value, default), AUCinf is AUCinf,obs.
#' If provided as clast.pred, AUCinf is AUCinf,pred.
#' @param lambda.z The elimination rate (in units of inverse time) for
#' extrapolation
#' @param options List of changes to the default
#' \code{\link{PKNCA.options}} for calculations.
#' @param method The method for integration (either 'lin up/log down'
#' or 'linear')
#' @param conc.blq How to handle BLQ values in between the first and
#' last above LOQ concentrations. (See \code{\link{clean.conc.blq}}
#' for usage instructions.)
#' @param conc.na How to handle missing concentration values.  (See
#' \code{\link{clean.conc.na}} for usage instructions.)
#' @param check Run \code{\link{check.conc.time}},
#' \code{\link{clean.conc.blq}}, and \code{\link{clean.conc.na}}?
#' @param fun.linear The function to use for integration of the linear
#' part of the curve (not required for AUC or AUMC functions)
#' @param fun.log The function to use for integration of the
#' logarithmic part of the curve (if log integration is used; not
#' required for AUC or AUMC functions)
#' @param fun.inf The function to use for extrapolation from the final
#' measurement to infinite time (not required for AUC or AUMC
#' functions.
#' @param ... For functions other than \code{pk.calc.auxc}, these
#' values are passed to \code{pk.calc.auxc}
#' @return A numeric value for the AU(M)C
#' @aliases pk.calc.auc pk.calc.aumc pk.calc.auc.last
#' @seealso \code{\link{pk.calc.auc.all}},
#' \code{\link{pk.calc.auc.last}}, \code{\link{clean.conc.blq}}
#' @references
#'
#' Gabrielsson J, Weiner D.  "Section 2.8.1 Computation methods -
#' Linear trapezoidal rule."  Pharmacokinetic & Pharmacodynamic Data
#' Analysis: Concepts and Applications, 4th Edition.  Stockholm,
#' Sweden: Swedish Pharmaceutical Press, 2000.  162-4.
#'
#' Gabrielsson J, Weiner D.  "Section 2.8.3 Computation methods -
#' Log-linear trapezoidal rule."  Pharmacokinetic & Pharmacodynamic
#' Data Analysis: Concepts and Applications, 4th Edition.  Stockholm,
#' Sweden: Swedish Pharmaceutical Press, 2000.  164-7.
#' @examples
#' myconc <- c(0, 1, 2, 1, 0.5, 0.25, 0)
#' mytime <- c(0, 1, 2, 3, 4,   5,    6)
#' pk.calc.auc(myconc, mytime, interval=c(0, 6))
#' pk.calc.auc(myconc, mytime, interval=c(0, Inf))
#' @export
pk.calc.auxc <- function(conc, time, interval=c(0, Inf),
                         clast=pk.calc.clast.obs(conc, time, check=FALSE),
                         lambda.z=NA,
                         auc.type="AUClast",
                         options=list(),
                         method=PKNCA.choose.option("auc.method", options),
                         conc.blq=PKNCA.choose.option("conc.blq", options),
                         conc.na=PKNCA.choose.option("conc.na", options),
                         check=TRUE,
                         fun.linear, fun.log, fun.inf) {
  ## Check the inputs
  method <- PKNCA.options(auc.method=method, check=TRUE)
  if (check) {
    check.conc.time(conc, time)
    data <- clean.conc.blq(conc, time,
                           conc.blq=conc.blq,
                           conc.na=conc.na,
                           check=FALSE)
  } else {
    data <- data.frame(conc, time)
  }
  if (nrow(data) == 0) {
    ## All the data were missing
    return(NA)
  } else if (all(data$conc %in% c(0, NA))) {
    ## All the data were missing or 0
    return(0)
  }
  auc.type <- match.arg(auc.type, c("AUClast", "AUCinf", "AUCall"))
  if (interval[1] >= interval[2])
    stop("The AUC interval must be increasing")
  if (auc.type %in% "AUCinf" &
        is.finite(interval[2]))
    warning("Requesting AUCinf when the end of the interval is not Inf")
  ##############################
  ## Subset the data to the range of interest
  interval.start <- interval[1]
  interval.end <- interval[2]
  ## Find the first time point
  if (interval.start < min(data$time)) {
    warning(sprintf(
      "Requesting an AUC range starting (%g) before the first measurement (%g) is not allowed",
      interval.start, min(data$time)))
    return(NA)
  } else if (interval.start > max(data$time)) {
    ## Give this as a warning, but allow it to continue
    warning(sprintf("AUC start time (%g) is after the maximum observed time (%g)",
                    interval.start, max(data$time)))
  }
  ## Ensure that we have clean concentration and time data.  This
  ## means that we need to make sure that we have our starting point.
  ## Interpolation ensures that (and will give the same answer if it
  ## already exists in the right form).
  conc.start <-
    interp.extrap.conc(data$conc, data$time,
                       time.out=interval.start,
                       lambda.z=lambda.z,
                       interp.method=method,
                       extrap.method=auc.type,
                       check=FALSE)
  ## Add that concentration and time to the vectors removing the
  ## original version if it was there.
  data <- rbind(data.frame(conc=conc.start,
                           time=interval.start),
                data[!(data$time %in% interval.start),])
  ## * either have our ending point or the ending point is Inf
  if (is.finite(interval.end)) {
    conc.end <-
      interp.extrap.conc(data$conc, data$time,
                         interval.end,
                         interp.method=method,
                         extrap.method=auc.type,
                         check=FALSE)
    ## !mask.end because we may be replacing an entry that is a 0.
    data <- rbind(data[!(data$time %in% interval.end),],
                  data.frame(conc=conc.end,
                             time=interval.end))
  }
  ## Subset the conc and time data to the interval of interest
  data <- subset(data, (interval.start <= time &
                          time <= interval.end))
  ## Set the overall tlast
  tlast <- pk.calc.tlast(data$conc, data$time, check=FALSE)
  if (is.na(tlast)) {
    ## All concentrations are BLQ (note that this has to be checked
    ## after full subsetting and interpolation to ensure that it is
    ## still true)
    if (all(data$conc %in% 0)) {
      ret <- 0
    } else {
      stop("Unknown error with NA tlast but non-BLQ concentrations")
    }
  } else {
    ## ############################
    ## Compute the AUxC
    ## Compute it in linear space from the start to Tlast
    if (auc.type %in% "AUCall" &
        tlast != max(data$time)) {
      ## Include the first point after tlast if it exists and we are
      ## looking for AUCall
      idx.1 <- 1:sum(data$time <= tlast)
    } else {
      idx.1 <- 1:(sum(data$time <= tlast) - 1)
    }
    idx.2 <- idx.1 + 1
    ret <- fun.linear(data$conc[idx.1], data$conc[idx.2],
                      data$time[idx.1], data$time[idx.2])
    if (method %in% "lin up/log down") {
      ## Compute log down if required (and if the later point is not 0)
      mask.down <- (data$conc[idx.2] < data$conc[idx.1] &
                      data$conc[idx.2] != 0)
      idx.1 <- idx.1[mask.down]
      idx.2 <- idx.2[mask.down]
      ret[mask.down] <- fun.log(data$conc[idx.1], data$conc[idx.2],
                                data$time[idx.1], data$time[idx.2])
    } else if (!(method %in% "linear")) {
      stop("Invalid AUC integration method")
    }
    if (auc.type %in% "AUCinf") {
      ## Whether AUCinf,obs or AUCinf,pred is calculated depends on if
      ## clast,obs or clast,pred is passed in.
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

#' @describeIn pk.calc.auxc Compute the AUCinf.
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

#' @describeIn pk.calc.auxc Compute the AUMCinf.
#' @export
pk.calc.aumc.inf <- function(conc, time, ..., options=list(),
                             lambda.z) {
  if ("auc.type" %in% names(list(...)))
    stop("auc.type cannot be changed when calling pk.calc.aumc.inf, please use pk.calc.aumc")
  pk.calc.aumc(conc=conc, time=time, ..., options=options,
               auc.type="AUCinf",
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
add.interval.col("aucinf",
                 FUN="pk.calc.auc.inf",
                 values=c(FALSE, TRUE),
                 desc="The area under the concentration time curve from the beginning of the interval to infinity",
                 depends="half.life")
PKNCA.set.summary("aucinf", business.geomean, business.geocv)

add.interval.col("auclast",
                 FUN="pk.calc.auc.last",
                 values=c(FALSE, TRUE),
                 desc="The area under the concentration time curve from the beginning of the interval to the last concentration above the limit of quantification",
                 depends=c())
PKNCA.set.summary("auclast", business.geomean, business.geocv)

add.interval.col("aucall",
                 FUN="pk.calc.auc.all",
                 values=c(FALSE, TRUE),
                 desc="The area under the concentration time curve from the beginning of the interval to the last concentration above the limit of quantification plus the triangle from that last concentration to 0 at the first concentration below the limit of quantification",
                 depends=c())
PKNCA.set.summary("aucall", business.geomean, business.geocv)

add.interval.col("aumcinf",
                 FUN="pk.calc.aumc.inf",
                 values=c(FALSE, TRUE),
                 desc="The area under the concentration time moment curve from the beginning of the interval to infinity",
                 depends="half.life")
PKNCA.set.summary("aumcinf", business.geomean, business.geocv)

add.interval.col("aumclast",
                 FUN="pk.calc.aumc.last",
                 values=c(FALSE, TRUE),
                 desc="The area under the concentration time moment curve from the beginning of the interval to the last concentration above the limit of quantification",
                 depends=c())
PKNCA.set.summary("aumclast", business.geomean, business.geocv)

add.interval.col("aumcall",
                 FUN="pk.calc.aumc.all",
                 values=c(FALSE, TRUE),
                 desc="The area under the concentration time moment curve from the beginning of the interval to the last concentration above the limit of quantification plus the moment of the triangle from that last concentration to 0 at the first concentration below the limit of quantification",
                 depends=c())
PKNCA.set.summary("aumcall", business.geomean, business.geocv)

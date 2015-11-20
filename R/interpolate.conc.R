#' Interpolate concentrations between measurements or extrapolate
#' concentrations after the last measurement.
#'
#' \code{interpolate.conc} and \code{extrapolate.conc} returns an
#' interpolated (or extrapolated) concentration.
#' \code{interp.extrap.conc} will choose whether interpolation or
#' extrapolation is required and will also operate on many
#' concentrations.  These will typically be used to estimate the
#' concentration between two measured concentrations or after the last
#' measured concentration.  Of note, these functions will not
#' extrapolate prior to the first point.
#'
#' @param conc Measured concentrations
#' @param time Time of the concentration measurement
#' @param time.out Time when interpolation is requested (vector for
#' \code{interp.extrap.conc}, scalar otherwise)
#' @param lambda.z The elimination rate constant.  \code{NA} will
#' prevent extrapolation.
#' @param clast The last observed concentration above the limit of
#' quantification.  If not given, \code{clast} is calculated from \code{\link{pk.calc.clast.obs}}
#' @param options List of changes to the default
#' \code{\link{PKNCA.options}} for calculations.
#' @param interp.method The method for interpolation (either
#' 'lin up/log down' or 'linear')
#' @param extrap.method The method for extrapolation: "AUCinf",
#' "AUClast", or "AUCall".  See details for usage.
#' @param conc.blq How to handle BLQ values. (See
#' \code{\link{clean.conc.blq}} for usage instructions.)
#' @param conc.na How to handle NA concentrations.  (See
#' \code{\link{clean.conc.na}})
#' @param check Run \code{\link{check.conc.time}},
#' \code{\link{clean.conc.blq}}, and \code{\link{clean.conc.na}}?
#' @return The interpolated or extrapolated concentration value as a
#' scalar float.
#'
#' @details
#' \describe{
#'   \item{extrap.method}{
#'     \describe{
#'       \item{'AUCinf'}{Use lambda.z to extrapolate beyond the last
#'         point with the half-life.}
#'       \item{'AUCall'}{If the last point is above the limit of
#'         quantification or missing, this is identical to 'AUCinf'.
#'         If the last point is below the limit of quantification,
#'         then linear interpolation between the Clast and the next
#'         BLQ is used for that interval and all additional points
#'         are extrapolated as 0.}
#'       \item{'AUClast'}{Extrapolates all points after the last above
#'         the limit of quantification as 0.}
#'     }
#'   }
#' }
#' @seealso \code{\link{pk.calc.clast.obs}}
#' \code{\link{pk.calc.half.life}}
#' @export
interp.extrap.conc <- function(conc, time, time.out,
                               lambda.z=NA,
                               clast=pk.calc.clast.obs(conc, time),
                               options=list(),
                               interp.method=PKNCA.choose.option("auc.method", options),
                               extrap.method="AUCinf",
                               conc.blq=PKNCA.choose.option("conc.blq", options),
                               conc.na=PKNCA.choose.option("conc.na", options),
                               check=TRUE) {
  if (check) {
    check.conc.time(conc, time)
    data <- clean.conc.blq(conc, time,
                           conc.blq=conc.blq, conc.na=conc.na,
                           check=FALSE)
  } else {
    data <- data.frame(conc, time)
  }
  tlast <- pk.calc.tlast(data$conc, data$time, check=FALSE)
  if (length(time.out) < 1)
    stop("time.out must be a vector with at least one element")
  ret <- rep(NA, length(time.out))
  for (i in seq_len(length(time.out)))
    if (is.na(time.out[i])) {
      warning("An interpolation/extrapolation time is NA")
    } else if (time.out[i] <= tlast) {
      ret[i] <- interpolate.conc(data$conc, data$time,
                                 time.out[i],
                                 interp.method=interp.method,
                                 conc.blq=conc.blq,
                                 conc.na=conc.na,
                                 check=FALSE)
    } else {
      ret[i] <- extrapolate.conc(data$conc, data$time,
                                 time.out[i],
                                 lambda.z=lambda.z,
                                 clast=clast,
                                 extrap.method=extrap.method,
                                 check=FALSE)
    }
  ret
}

#' @describeIn interp.extrap.conc Interpolate concentrations through
#' Tlast (inclusive)
#' @export
interpolate.conc <- function(conc, time, time.out,
                             options=list(),
                             interp.method=PKNCA.choose.option("auc.method", options),
                             conc.blq=PKNCA.choose.option("conc.blq", options),
                             conc.na=PKNCA.choose.option("conc.na", options),
                             check=TRUE) {
  ## Check the inputs
  if (check) {
    check.conc.time(conc, time)
    data <- clean.conc.blq(conc, time,
                           conc.blq=conc.blq, conc.na=conc.na,
                           check=FALSE)
  } else {
    data <- data.frame(conc, time)
  }
  ## Verify that we are interpolating between the first concentration
  ## and the last above LOQ concentration
  if (length(time.out) != 1)
    stop("Can only interpolate for one time point per function call")
  tlast <- pk.calc.tlast(data$conc, data$time, check=FALSE)
  if (time.out < min(data$time))
    stop("Cannot interpolate backward in time")
  if (time.out > tlast)
    stop("interpolate.conc can only works through Tlast, please use interp.extrap.conc to combine both interpolation and extrapolation.")
  if (!(tolower(interp.method) %in% c("lin up/log down", "linear")))
    stop("interp.method must be one of 'linear' or 'lin up/log down'")
  keep.trying <- TRUE
  if (time.out %in% data$time) {
    ## See if there is an exact time match and return that if it
    ## exists.
    ret <- conc[time.out == data$time]
  } else {
    ## Find the last time before and the first time after the output
    ## time.
    time.1 <- max(data$time[data$time <= time.out])
    time.2 <- min(data$time[time.out <= data$time])
    conc.1 <- data$conc[data$time == time.1]
    conc.2 <- data$conc[data$time == time.2]
    interp.method <- tolower(interp.method)
    if ((interp.method == "linear") |
        (interp.method == "lin up/log down" &
         ((conc.1 <= 0 | conc.2 <= 0) |
          (conc.1 <= conc.2)))) {
      ## Do linear interpolation if:
      ##   linear interpolation is selected or
      ##   lin up/log down interpolation is selected and
      ##     one concentration is 0 or
      ##     the concentrations are equal
      ret <- conc.1+(time.out-time.1)/(time.2-time.1)*(conc.2-conc.1)
    } else if (interp.method == "lin up/log down") {
      ret <- exp(log(conc.1)+
                 (time.out-time.1)/(time.2-time.1)*(log(conc.2)-log(conc.1)))
    } else {
      stop("You should never see this error.")
    }
  }
  ret
}

#' @describeIn interp.extrap.conc Extrapolate concentrations after Tlast
#' @export
extrapolate.conc <- function(conc, time, time.out,
                             lambda.z=NA, clast=pk.calc.clast.obs(conc, time),
                             extrap.method="AUCinf",
                             options=list(),
                             conc.na=PKNCA.choose.option("conc.na", options),
                             conc.blq=PKNCA.choose.option("conc.blq", options),
                             check=TRUE) {
  if (check) {
    check.conc.time(conc, time)
    data <- clean.conc.blq(conc, time, conc.na=conc.na, check=FALSE)
  } else {
    data <- data.frame(conc, time)
  }
  extrap.method <- tolower(extrap.method)
  if (!(extrap.method %in%
        c("aucinf", "aucall", "auclast")))
    stop("extrap.method must be one of 'AUCinf', 'AUClast', or 'AUCall'")
  if (length(time.out) != 1)
    stop("Only one time.out value may be estimated at once.")
  tlast <- pk.calc.tlast(data$conc, data$time, check=FALSE)
  if (is.na(tlast)) {
    ## If there are no observed concentrations, return NA
    ret <- NA
  } else if (time.out <= tlast) {
    stop("extrapolate.conc can only work beyond Tlast, please use interp.extrap.conc to combine both interpolation and extrapolation.")
  } else {
    ## Start the interpolation
    if (extrap.method %in% "aucinf") {
      ## If AUCinf is requested, extrapolate using the half-life
      ret <- clast*exp(-lambda.z*(time.out - tlast))
    } else if (extrap.method %in% "auclast" |
                 (extrap.method %in% "aucall" &
                    tlast == max(data$time))) {
      ## If AUClast is requested or AUCall is requested and there are
      ## no BLQ at the end, we are already certain that we are after
      ## Tlast, so the answer is 0.
      ret <- 0
    } else if (extrap.method %in% "aucall") {
      ## If the last non-missing concentration is below the limit of
      ## quantification, extrapolate with the triangle method of
      ## AUCall.
      time.prev <- max(data$time[data$time <= time.out])
      conc.prev <- data$conc[data$time %in% time.prev]
      if (conc.prev %in% 0) {
        ## If we are already BLQ, then we have confirmed that there
        ## are no more ALQ measurements (because we are beyond
        ## Tlast) and therefore we can extrapolate as 0.
        ret <- 0
      } else {
        if (time.prev != max(data$time)) {
          time.next <- min(data$time[data$time >= time.out])
          conc.next <- data$conc[data$time %in% time.next]
        }
        ## If we are not already BLQ, then we have confirmed that we
        ## are in the triangle extrapolation region and need to draw
        ## a line.
        ret <- (time.out - time.prev)/(time.next - time.prev)*conc.prev
      }
    } else {
      stop("Invalid extrap.method caught too late (seeing this error indicates a software bug)")
    }
  }
  ret
}

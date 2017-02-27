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
#' @param time.dose Time of the dose
#' @param time.out Time when interpolation is requested (vector for 
#'   \code{interp.extrap.conc}, scalar otherwise)
#' @param lambda.z The elimination rate constant.  \code{NA} will 
#'   prevent extrapolation.
#' @param clast The last observed concentration above the limit of 
#'   quantification.  If not given, \code{clast} is calculated from 
#'   \code{\link{pk.calc.clast.obs}}
#' @param conc.origin The concentration before the first measurement. 
#'   \code{conc.origin} is typically used to set predose values to zero 
#'   (default), set a predose concentration for endogenous compounds, or
#'   set predose concentrations to \code{NA} if otherwise unknown.
#' @param options List of changes to the default 
#'   \code{\link{PKNCA.options}} for calculations.
#' @param interp.method The method for interpolation (either 'lin up/log
#'   down' or 'linear')
#' @param extrap.method The method for extrapolation: "AUCinf", 
#'   "AUClast", or "AUCall".  See details for usage.
#' @param conc.blq How to handle BLQ values. (See 
#'   \code{\link{clean.conc.blq}} for usage instructions.)
#' @param conc.na How to handle NA concentrations.  (See 
#'   \code{\link{clean.conc.na}})
#' @param route.dose What is the route of administration 
#'   ("intravascular" or "extravascular").  See the details below for 
#'   how this parameter is used.
#' @param duration.dose What is the duration of administration? See the 
#'   details below for how this parameter is used.
#' @param out.after Should interpolation occur from the data before 
#'   (\code{FALSE}) or after (\code{TRUE}) the interpolated point?  See 
#'   the details below for how this parameter is used.  It only has a
#'   meaningful effect at the instant of an IV bolus dose.
#' @param check Run \code{\link{check.conc.time}}, 
#'   \code{\link{clean.conc.blq}}, and \code{\link{clean.conc.na}}?
#' @param ... Additional arguments passed to \code{interpolate.conc} or 
#'   \code{extrapolate.conc}.
#' @return The interpolated or extrapolated concentration value as a 
#'   scalar float.
#' @details
#' \describe{
#'   \item{extrap.method}{
#'     \describe{
#'       \item{'AUCinf'}{Use lambda.z to extrapolate beyond the last point with the half-life.}
#'       \item{'AUCall'}{If the last point is above the limit of quantification or missing, this is identical to 'AUCinf'. If the last point is below the limit of quantification, then linear interpolation between the Clast and the next BLQ is used for that  interval and all additional points are extrapolated as 0.}
#'       \item{'AUClast'}{Extrapolates all points after the last above the  limit of quantification as 0.}
#'     }
#'   }
#' }
#' 
#' \code{duration.dose} and \code{direction.out} are ignored if
#' \code{route.dose == "extravascular"}.  \code{direction.out} is ignored
#' if \code{duration.dose > 0}.
#'
#' \code{route.dose} and \code{duration.dose} affect how 
#' interpolation/extrapolation of the concentration occurs at the time 
#' of dosing.  If \code{route.dose == "intravascular"} and 
#' \code{duration.dose == 0} then extrapolation occurs for an IV bolus 
#' using \code{\link{pk.calc.c0}} with the data after dosing.  Otherwise
#' (either \code{route.dose == "extravascular"} or \code{duration.dose >
#' 0}), extrapolation occurs using the concentrations before dosing and 
#' estimating the half-life (or more precisely, estimating 
#' \code{lambda.z}).  Finally, \code{direction.out} can change the
#' direction of interpolation in cases with \code{route.dose ==
#' "intravascular"} and \code{duration.dose == 0}.  When
#' \code{direction.out == "before"} interpolation occurs only with data
#' before the dose (as is the case for \code{route.dose ==
#' "extravascular"}), but if \code{direction.out == "after"}
#' interpolation occurs from the data after dosing.
#' 
#' @seealso \code{\link{pk.calc.clast.obs}},
#'   \code{\link{pk.calc.half.life}}, \code{\link{pk.calc.c0}}
#' @export
interp.extrap.conc <- function(conc, time, time.out,
                               lambda.z=NA,
                               clast=pk.calc.clast.obs(conc, time),
                               options=list(),
                               interp.method=PKNCA.choose.option("auc.method", options),
                               extrap.method="AUCinf",
                               ...,
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

#' @describeIn interp.extrap.conc Interpolate concentrations through Tlast (inclusive)
#' @export
interpolate.conc <- function(conc, time, time.out,
                             options=list(),
                             interp.method=PKNCA.choose.option("auc.method", options),
                             conc.blq=PKNCA.choose.option("conc.blq", options),
                             conc.na=PKNCA.choose.option("conc.na", options),
                             conc.origin=0,
                             ...,
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
  # Ensure that conc.origin is valid
  if (length(conc.origin) != 1) {
    stop("conc.origin must be a scalar")
  }
  if (!(is.na(conc.origin) | (is.numeric(conc.origin) & !is.factor(conc.origin)))) {
    stop("conc.origin must be NA or a number (and not a factor)")
  }
  ## Verify that we are interpolating between the first concentration
  ## and the last above LOQ concentration
  if (length(time.out) != 1) {
    stop("Can only interpolate for one time point per function call")
  }
  tlast <- pk.calc.tlast(data$conc, data$time, check=FALSE)
  if (time.out < min(data$time)) {
    ret <- conc.origin
  } else if (time.out > tlast) {
    stop("interpolate.conc can only works through Tlast, please use interp.extrap.conc to combine both interpolation and extrapolation.")
  } else if (!(tolower(interp.method) %in% c("lin up/log down", "linear"))) {
    stop("interp.method must be one of 'linear' or 'lin up/log down'")
  } else if (time.out %in% data$time) {
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
      stop("You should never see this error.  Please report this as a bug with a reproducible example.") # nocov
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
                             ...,
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
      stop("Invalid extrap.method caught too late (seeing this error indicates a software bug)") # nocov
    }
  }
  ret
}

#' @describeIn interp.extrap.conc Interpolate and extrapolate 
#'   concentrations without interpolating or extrapolating beyond doses.
#' @export
interp.extrap.conc.dose <- function(conc, time,
                                    time.dose, route.dose="extravascular", duration.dose=NA,
                                    time.out, out.after=FALSE,
                                    conc.blq=PKNCA.choose.option("conc.blq", options),
                                    conc.na=PKNCA.choose.option("conc.na", options),
                                    ...,
                                    check=TRUE) {
  if (check) {
    check.conc.time(conc, time)
    data.conc <-
      clean.conc.blq(conc, time,
                     conc.blq=conc.blq, conc.na=conc.na,
                     check=FALSE)
  } else {
    data.conc <- data.frame(conc, time)
  }
  # Check other inputs
  if (!is.character(route.dose)) {
    route.dose <- as.character(route.dose)
  }
  if (!(all(route.dose %in% c("extravascular", "intravascular")))) {
    stop("route.dose must be either 'extravascular' or 'intravascular'")
  }
  if (!(length(route.dose) %in% c(1, length(time.dose)))) {
    stop("route.dose must either be a scalar or the same length as time.dose")
  }
  if (!all(is.na(duration.dose) | (is.numeric(duration.dose) & !is.factor(duration.dose)))) {
    stop("duration.dose must be NA or a number.")
  }
  if (!(length(duration.dose) %in% c(1, length(time.dose)))) {
    stop("duration.dose must either be a scalar or the same length as time.dose")
  }

  ## Separate each data set by time interval
  dosebreaks <- c(-Inf, time.dose, Inf)
  data.conc$cut <- cut(data.conc$time, breaks=dosebreaks, right=FALSE)
  data.conc$after <- FALSE
  data.dose <-
    data.frame(time=time.dose,
               route=route.dose,
               duration=duration.dose,
               cut=cut(time.dose, breaks=dosebreaks, right=FALSE))
  data.dose$iv.bolus <-
    data.dose$route %in% "intravascular" & data.dose$duration %in% 0
  data.out <-
    data.frame(
      time=time.out,
      after=out.after,
      cut=cut(time.out, breaks=dosebreaks, right=FALSE),
      conc=NA,
      # "conc", "dose", "both", or NA for the type of data that is
      # before, at the same time of, or after the output requested.
      event.before=NA,
      event.at=NA,
      event.after=NA,
      # The time that the event before or after occurs on.
      time.before=NA,
      time.after=NA,
      # If the event.before, event.at, or event.after was either "dose"
      # or "both", was the dose an IV bolus?
      iv.bolus.before=NA,
      iv.bolus.at=NA,
      iv.bolus.after=NA,
      Method="",
      stringsAsFactors=FALSE)
  # Find the type of event that occurs at each time point. ####
  all.input.times <-
    rbind(
      data.frame(
        eventtype="conc",
        time=unique(data.conc$time),
        stringsAsFactors=FALSE),
      data.frame(
        eventtype="dose",
        time=unique(data.dose$time),
        stringsAsFactors=FALSE))
  all.input.times$eventtype[
    all.input.times$time %in%
      all.input.times$time[duplicated(all.input.times$time)]] <-
    "both"
  all.input.times$iv.bolus <-
    all.input.times$time %in% data.dose$time[data.dose$iv.bolus]
  all.input.times <- unique(all.input.times)
  # Find the events that occur before, at, and after the currently requested output time. ####
  for (i in seq_len(nrow(data.out))) {
    if (any(mask.before <- all.input.times$time < data.out$time[i])) {
      # Determine what (if anything) happened before the current time.
      data.out$time.before[i] <- max(all.input.times$time[mask.before])
      data.out$event.before[i] <-
        all.input.times$eventtype[all.input.times$time %in% data.out$time.before[i]]
      if (data.out$event.before[i] %in% c("dose", "both")) {
        data.out$iv.bolus.before[i] <-
          data.dose$iv.bolus[data.dose$time %in% data.out$time.before[i]]
      }
    }
    if (any(mask.at <- all.input.times$time %in% data.out$time[i])) {
      # Determine what (if anything) happened at the current time.
      data.out$event.at[i] <- all.input.times$eventtype[mask.at]
      if (data.out$event.at[i] %in% c("dose", "both")) {
        data.out$iv.bolus.at[i] <-
          data.dose$iv.bolus[data.dose$time %in% data.out$time[i]]
      }
    }
    if (any(mask.after <- all.input.times$time > data.out$time[i])) {
      # Determine what (if anything) happened after the current time.
      data.out$time.after[i] <- min(all.input.times$time[mask.after])
      data.out$event.after[i] <-
        all.input.times$eventtype[all.input.times$time %in% data.out$time.after[i]]
      if (data.out$event.after[i] %in% c("dose", "both")) {
        data.out$iv.bolus.after[i] <-
          data.dose$iv.bolus[data.dose$time %in% data.out$time.after[i]]
      }
    }
  }
  # Loop through the methods until all have been tested or no missing
  # values remain.
  for (n in names(interp.extrap.conc.dose.select)) {
    mask <- is.na(data.out$conc) &
      interp.extrap.conc.dose.select[[n]]$select(data.out)
    if (any(mask)) {
      if ("warning" %in% names(interp.extrap.conc.dose.select[[n]])) {
        warning(sprintf("%s: %d data points",
                        interp.extrap.conc.dose.select[[n]]$warning,
                        sum(mask)))
        data.out$Method[mask] <- n
      } else {
        for (current.idx in which(mask)) {
          data.out$conc[current.idx] <-
            interp.extrap.conc.dose.select[[n]]$value(data.out[current.idx,],
                                                      data.conc, data.dose,
                                                      ...)
          data.out$Method[current.idx] <- n
        }
      }
    }
  }
  ret <- data.out$conc
  attr(ret, "Method") <- data.out$Method
  ret
}

interp.extrap.conc.dose.select <-
  list(

    # Copy ####
    "Copy"=
      list(
        select=function(x) {
          (x$event.at %in% "conc" |
             (x$event.at %in% "both" &
                (!x$iv.bolus.at |
                   (x$iv.bolus.at & !x$after))))
        },
        value=function(data.out, data.conc, data.dose, ...) {
          data.conc$conc[data.conc$time %in% data.out$time]
        },
        description="Copy the input concentration at the given time to the output."),

    # Before all events ####
    "Before all events"=
      list(
        select=function(x) {
          is.na(x$event.before)
        },
        value=function(data.out, data.conc, data.dose, conc.origin=0, ...) {
          conc.origin
        },
        description=paste("Interpolation before any events is NA or zero (0)",
                          "depending on the value of conc.origin.  conc.origin",
                          "defaults to zero which is the implicit assumption",
                          "that a complete washout occurred and there is no",
                          "endogenous source of the analyte.")),

    # Interpolation ####
    "Interpolation"=
      list(
        select=function(x) {
          ((x$event.before %in% "conc" |
              (x$event.before %in% "both" &
                 !(x$iv.bolus.before %in% TRUE))) &
             is.na(x$event.at) &
             (x$event.after %in% c("conc", "both")))
        },
        value=function(data.out, data.conc, data.dose, ...) {
          tmp.conc <- data.conc[data.conc$cut %in% data.out$cut,]
          interpolate.conc(conc=tmp.conc$conc, time=tmp.conc$time,
                           time.out=data.out$time,
                           check=FALSE, ...)
        },
        description=paste("With concentrations before and after and not an IV",
                          "bolus before, interpolate between observed concentrations.")),

    # Immediately after an IV bolus with a concentration next ####
    "Immediately after an IV bolus with a concentration next"=
      list(
        select=function(x) {
          (x$iv.bolus.at %in% TRUE) & x$after & (x$event.after %in% "conc")
        },
        value=function(data.out, data.conc, data.dose, ...) {
          tmp.conc <- data.conc[data.conc$cut %in% data.out$cut,]
          tmp.dose <- data.dose[data.dose$cut %in% data.out$cut,]
          pk.calc.c0(conc=tmp.conc$conc, time=tmp.conc$time,
                     time.dose=tmp.dose$time,
                     method=c("logslope", "c1"))
        },
        description=paste("Calculate C0 for the time immediately after an IV",
                          "bolus.  First, attempt using log slope",
                          "back-extrapolation.  If that fails, use the first",
                          "concentration after the dose as C0.")),

    # Immediately after an IV bolus without a concentration next ####
    "Immediately after an IV bolus without a concentration next"=
      list(
        select=function(x) {
          (x$iv.bolus.at %in% TRUE) & x$after & !(x$event.after %in% "conc")
        },
        warning="Cannot interpolate immediately after an IV bolus without a concentration next.",
        description="Cannot calculate C0 without a concentration after an IV bolus; return NA."),

    # After an IV bolus with a concentration next ####
    "After an IV bolus with a concentration next"=
      list(
        select=function(x) {
          (x$iv.bolus.before %in% TRUE) & is.na(x$event.at) & (x$event.after %in% "conc")
        },
        value=function(data.out, data.conc, data.dose, ...) {
          tmp.conc <- data.conc[data.conc$cut %in% data.out$cut, c("conc", "time")]
          tmp.dose <- data.dose[data.dose$cut %in% data.out$cut, "time", drop=FALSE]
          tmp.conc <- rbind(
            data.frame(
              time=tmp.dose$time,
              conc=
                pk.calc.c0(conc=tmp.conc$conc, time=tmp.conc$time,
                           time.dose=tmp.dose$time,
                           method=c("logslope", "c1"))),
            tmp.conc)
          interpolate.conc(conc=tmp.conc$conc, time=tmp.conc$time,
                           time.out=data.out$time,
                           check=FALSE, ...)
        },
        description=paste("First, calculate C0 using log slope back-extrapolation",
                          "(falling back to the first post-dose concentration",
                          "if that fails).  Then, interpolate between C0 and",
                          "the first post-dose concentration.")),

    # After an IV bolus without a concentration next ####
    "After an IV bolus without a concentration next"=
      list(
        select=function(x) {
          (x$iv.bolus.before %in% TRUE) & is.na(x$event.at) & !(x$event.after %in% "conc")
        },
        warning="Cannot interpolate after an IV bolus without a concentration next.",
        description=paste("Between an IV bolus and anything other than a",
                          "concentration, interpolation cannot occur.  Return NA")),

    # After a dose with no event and a concentration after ####
    "After a dose with no event and a concentration after"=
      list(
        select=function(x) {
          (x$event.before %in% "dose") & (x$iv.bolus.before %in% FALSE) &
            is.na(x$event.at) &
            (x$event.after %in% "conc")
        },
        value=function(data.out, data.conc, data.dose, ...) {
          # Extrapolate the concentration to the dose time
          tmp.conc <- data.conc[as.numeric(data.conc$cut) %in% (as.numeric(data.out$cut) - 1),]
          if (tmp.conc$conc[nrow(tmp.conc)] == 0) {
            tmp.conc2 <-
              rbind(data.frame(conc=0,
                               time=data.dose$time),
                    data.conc[data.conc$cut %in% data.out$cut, c("conc", "time")])
          } else {
            tmp.conc2 <-
              rbind(data.frame(conc=extrapolate.conc(conc=tmp.conc$conc, time=tmp.conc$time,
                                                     time.out=data.dose$time,
                                                     ..., check=FALSE),
                               time=data.dose$time),
                    data.conc[data.conc$cut %in% data.out$cut, c("conc", "time")])
          }
          # Interpolate between the dose time and the current time
          interpolate.conc(conc=tmp.conc2$conc, time=tmp.conc2$time, time.out=data.out$time,
                           ..., check=FALSE)
        },
        description=paste("If there is a non-IV bolus dose before, nothing",
                          "at the current time, and a concentration after,",
                          "extrapolate from prior data to the concentration",
                          "at the time of dosing then interpolate between the",
                          "dosing time and the requested output time.")
        ),
    # Doses with no concentrations between ####
    "Doses with no concentrations between"=
      list(
        select=function(x) {
          (x$event.before %in% c("dose", "both")) &
            (x$event.at %in% c("dose", NA))
        },
        warning="Cannot interpolate between two doses or after a dose without a concentration after the first dose.",
        description="Two doses with no concentrations between them, return NA."),

    # Dose as the last event and requesting a concentration after ####
    "Dose as the last event and requesting a concentration after"=
      list(
        select=function(x) {
          (x$event.before %in% c("dose", "both")) &
            is.na(x$event.at)
        },
        warning="Cannot extrapolate from a dose without any concentrations after it.",
        description=paste("Cannot estimate the concentration after a dose",
                          "without concentrations after the dose, return NA.")),

    # After all events, concentration before ####
    "After all events, concentration before"=
      list(
        select=function(x) {
          is.na(x$event.after) & (x$event.before %in% "conc")
        },
        value=function(data.out, data.conc, data.dose, ...) {
          tmp.conc <- data.conc[data.conc$cut %in% data.out$cut,]
          extrapolate.conc(conc=tmp.conc$conc, time=tmp.conc$time,
                           time.out=data.out$time, check=FALSE, ...)
        },
        description="Extrapolate from the concentration data after the last dose."),

    # Concentration, nothing or dose but not(dose, IV bolus, After), nothing or dose ####
    "Concentration, nothing or dose but not(dose, IV bolus, After), nothing or dose"=
      list(
        select=function(x) {
          (x$event.before %in% "conc") &
            (is.na(x$event.at) |
               ((x$event.at %in% "dose") & !((x$iv.bolus.at %in% TRUE) & x$after)))
        },
        value=function(data.out, data.conc, data.dose, ...) {
          if (data.out$event.at %in% "dose") {
            ## Shift backward one interval because the dose is the start
            ## of an interval.
            tmp.conc <- data.conc[as.numeric(data.conc$cut) %in% (as.numeric(data.out$cut) - 1),]
          } else {
            tmp.conc <- data.conc[data.conc$cut %in% data.out$cut,]
          }
          if (tmp.conc$conc[nrow(tmp.conc)] == 0) {
            0
          } else {
            # Solve for the half-life if required
            args <- list(conc=tmp.conc$conc, time=tmp.conc$time,
                         time.out=data.out$time, check=FALSE, ...)
            if (!("lambda.z" %in% names(args))) {
              args$lambda.z <- pk.calc.half.life(conc=tmp.conc$conc, time=tmp.conc$time)$lambda.z
            }
            do.call(extrapolate.conc, args)
          }
        },
        description=paste("Extrapolate when the previous measurement is a",
                          "concentration, and the current event is a dose or",
                          "nothing.  If the previous concentration measurement",
                          "is below the limit of quantification (zero), then",
                          "return zero.")))

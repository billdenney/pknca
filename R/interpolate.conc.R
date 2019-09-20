#' Interpolate concentrations between measurements or extrapolate concentrations
#' after the last measurement.
#' 
#' \code{interpolate.conc()} and \code{extrapolate.conc()} returns an
#' interpolated (or extrapolated) concentration. \code{interp.extrap.conc()}
#' will choose whether interpolation or extrapolation is required and will also
#' operate on many concentrations.  These will typically be used to estimate the
#' concentration between two measured concentrations or after the last measured
#' concentration. Of note, these functions will not extrapolate prior to the
#' first point.
#' 
#' @param conc Measured concentrations
#' @param time Time of the concentration measurement
#' @param time.dose Time of the dose
#' @param time.out Time when interpolation is requested (vector for
#'   \code{interp.extrap.conc()}, scalar otherwise)
#' @param lambda.z The elimination rate constant.  \code{NA} will prevent
#'   extrapolation.
#' @param clast The last observed concentration above the limit of
#'   quantification.  If not given, \code{clast} is calculated from
#'   \code{\link{pk.calc.clast.obs}()}
#' @param conc.origin The concentration before the first measurement.
#'   \code{conc.origin} is typically used to set predose values to zero
#'   (default), set a predose concentration for endogenous compounds, or set
#'   predose concentrations to \code{NA} if otherwise unknown.
#' @param options List of changes to the default \code{\link{PKNCA.options}()}
#'   for calculations.
#' @param interp.method The method for interpolation (either "lin up/log down"
#'   or "linear")
#' @param extrap.method The method for extrapolation: "AUCinf", "AUClast", or
#'   "AUCall".  See details for usage.
#' @param conc.blq How to handle BLQ values. (See \code{\link{clean.conc.blq}()}
#'   for usage instructions.)
#' @param conc.na How to handle NA concentrations.  (See
#'   \code{\link{clean.conc.na}()})
#' @param route.dose What is the route of administration ("intravascular" or
#'   "extravascular").  See the details for how this parameter is used.
#' @param duration.dose What is the duration of administration? See the details
#'   for how this parameter is used.
#' @param out.after Should interpolation occur from the data before
#'   (\code{FALSE}) or after (\code{TRUE}) the interpolated point?  See the
#'   details for how this parameter is used.  It only has a meaningful effect at
#'   the instant of an IV bolus dose.
#' @param check Run \code{\link{check.conc.time}()},
#'   \code{\link{clean.conc.blq}()}, and \code{\link{clean.conc.na}()}?
#' @param ... Additional arguments passed to \code{interpolate.conc()} or
#'   \code{extrapolate.conc()}.
#' @return The interpolated or extrapolated concentration value as a scalar
#'   double (or vector for \code{interp.extrap.conc()}).
#'
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
#' \code{duration.dose} and \code{direction.out} are ignored if \code{route.dose
#' == "extravascular"}.  \code{direction.out} is ignored if \code{duration.dose
#' > 0}.
#'
#' \code{route.dose} and \code{duration.dose} affect how
#' interpolation/extrapolation of the concentration occurs at the time of
#' dosing.  If \code{route.dose == "intravascular"} and \code{duration.dose ==
#' 0} then extrapolation occurs for an IV bolus using \code{\link{pk.calc.c0}()}
#' with the data after dosing.  Otherwise (either \code{route.dose ==
#' "extravascular"} or \code{duration.dose > 0}), extrapolation occurs using the
#' concentrations before dosing and estimating the half-life (or more precisely,
#' estimating \code{lambda.z}).  Finally, \code{direction.out} can change the
#' direction of interpolation in cases with \code{route.dose == "intravascular"}
#' and \code{duration.dose == 0}.  When \code{direction.out == "before"}
#' interpolation occurs only with data before the dose (as is the case for
#' \code{route.dose == "extravascular"}), but if \code{direction.out == "after"}
#' interpolation occurs from the data after dosing.
#' 
#' @seealso \code{\link{pk.calc.clast.obs}()},
#'   \code{\link{pk.calc.half.life}()}, \code{\link{pk.calc.c0}()}
#' @export
interp.extrap.conc <- function(conc, time, time.out,
                               lambda.z=NA,
                               clast=pk.calc.clast.obs(conc, time),
                               options=list(),
                               interp.method=NULL,
                               extrap.method="AUCinf",
                               ...,
                               conc.blq=NULL,
                               conc.na=NULL,
                               check=TRUE) {
  # Check inputs
  interp.method <- PKNCA.choose.option(name="auc.method", value=interp.method, options=options)
  conc.blq <- PKNCA.choose.option(name="conc.blq", value=conc.blq, options=options)
  conc.na <- PKNCA.choose.option(name="conc.na", value=conc.na, options=options)
  if (check) {
    check.conc.time(conc, time)
    data <-
      clean.conc.blq(
        conc, time,
        conc.blq=conc.blq, conc.na=conc.na,
        check=FALSE
      )
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
      ret[i] <-
        interpolate.conc(
          conc=data$conc, time=data$time,
          time.out=time.out[i],
          interp.method=interp.method,
          conc.blq=conc.blq,
          conc.na=conc.na,
          check=FALSE
        )
    } else {
      ret[i] <-
        extrapolate.conc(
          conc=data$conc, time=data$time,
          time.out=time.out[i],
          lambda.z=lambda.z,
          clast=clast,
          extrap.method=extrap.method,
          check=FALSE
        )
    }
  ret
}

#' @describeIn interp.extrap.conc Interpolate concentrations through Tlast (inclusive)
#' @export
interpolate.conc <- function(conc, time, time.out,
                             options=list(),
                             interp.method=NULL,
                             conc.blq=NULL,
                             conc.na=NULL,
                             conc.origin=0,
                             ...,
                             check=TRUE) {
  ## Check the inputs
  interp.method <- PKNCA.choose.option(name="auc.method", value=interp.method, options=options)
  conc.blq <- PKNCA.choose.option(name="conc.blq", value=conc.blq, options=options)
  conc.na <- PKNCA.choose.option(name="conc.na", value=conc.na, options=options)
  if (check) {
    check.conc.time(conc, time)
    data <-
      clean.conc.blq(
        conc=conc, time=time,
        conc.blq=conc.blq, conc.na=conc.na,
        check=FALSE
      )
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
  tlast <- pk.calc.tlast(conc=data$conc, time=data$time, check=FALSE)
  if (time.out < min(data$time)) {
    ret <- conc.origin
  } else if (time.out > tlast) {
    stop("`interpolate.conc()` can only works through Tlast, please use `interp.extrap.conc()` to combine both interpolation and extrapolation.")
  } else if (time.out %in% data$time) {
    ## See if there is an exact time match and return that if it
    ## exists.
    ret <- conc[time.out == data$time]
  } else {
    ## Find the last time before and the first time after the output
    ## time.
    time_1 <- max(data$time[data$time <= time.out])
    time_2 <- min(data$time[time.out <= data$time])
    conc_1 <- data$conc[data$time == time_1]
    conc_2 <- data$conc[data$time == time_2]
    interp.method <- tolower(interp.method)
    if (is.na(conc_1) | is.na(conc_2)) {
      ret <- NA_real_
    } else if ((interp.method == "linear") |
        (interp.method == "lin up/log down" &
         ((conc_1 <= 0 | conc_2 <= 0) |
          (conc_1 <= conc_2)))) {
      ## Do linear interpolation if:
      ##   linear interpolation is selected or
      ##   lin up/log down interpolation is selected and
      ##     one concentration is 0 or
      ##     the concentrations are equal
      ret <- conc_1+(time.out-time_1)/(time_2-time_1)*(conc_2-conc_1)
    } else if (interp.method == "lin up/log down") {
      ret <-
        exp(
          log(conc_1)+
            (time.out-time_1)/(time_2-time_1)*(log(conc_2)-log(conc_1))
        )
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
                             conc.na=NULL,
                             conc.blq=NULL,
                             ...,
                             check=TRUE) {
  conc.na <- PKNCA.choose.option(name="conc.na", value=conc.na, options=options)
  conc.blq <- PKNCA.choose.option(name="conc.blq", value=conc.blq, options=options)
  if (check) {
    check.conc.time(conc, time)
    data <-
      clean.conc.blq(
        conc=conc, time=time,
        conc.na=conc.na,
        check=FALSE
      )
  } else {
    data <- data.frame(conc, time)
  }
  extrap.method <- tolower(extrap.method)
  if (!(extrap.method %in% c("aucinf", "aucall", "auclast")))
    stop("extrap.method must be one of 'AUCinf', 'AUClast', or 'AUCall'")
  if (length(time.out) != 1)
    stop("Only one time.out value may be estimated at once.")
  tlast <- pk.calc.tlast(conc=data$conc, time=data$time, check=FALSE)
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
      time_prev <- max(data$time[data$time <= time.out])
      conc_prev <- data$conc[data$time %in% time_prev]
      if (conc_prev %in% 0) {
        ## If we are already BLQ, then we have confirmed that there
        ## are no more ALQ measurements (because we are beyond
        ## Tlast) and therefore we can extrapolate as 0.
        ret <- 0
      } else {
        if (time_prev != max(data$time)) {
          time_next <- min(data$time[data$time >= time.out])
        }
        ## If we are not already BLQ, then we have confirmed that we
        ## are in the triangle extrapolation region and need to draw
        ## a line.
        ret <- (time.out - time_prev)/(time_next - time_prev)*conc_prev
      }
    } else {
      stop("Invalid extrap.method caught too late (seeing this error indicates a software bug)") # nocov
    }
  }
  ret
}

# Choices for events in interp.extrap.conc.dose.  This is included here to assist with testing later.
event_choices_interp.extrap.conc.dose <-
  list(conc_dose_iv_bolus_after="conc_dose_iv_bolus_after",
       conc_dose="conc_dose",
       dose_iv_bolus_after="dose_iv_bolus_after",
       dose="dose",
       conc="conc",
       output_only="output_only",
       none="none")

#' @importFrom dplyr case_when
#' @describeIn interp.extrap.conc Interpolate and extrapolate 
#'   concentrations without interpolating or extrapolating beyond doses.
#' @export
interp.extrap.conc.dose <- function(conc, time,
                                    time.dose, route.dose="extravascular", duration.dose=NA,
                                    time.out, out.after=FALSE,
                                    options=list(),
                                    conc.blq=NULL,
                                    conc.na=NULL,
                                    ...,
                                    check=TRUE) {
  # Check inputs
  conc.na <- PKNCA.choose.option(name="conc.na", value=conc.na, options=options)
  conc.blq <- PKNCA.choose.option(name="conc.blq", value=conc.blq, options=options)
  if (check) {
    check.conc.time(conc, time)
    data_conc <-
      clean.conc.blq(conc, time,
                     conc.blq=conc.blq, conc.na=conc.na,
                     check=FALSE)
  } else {
    data_conc <- data.frame(conc, time)
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

  # Generate a single timeline
  
  # Concentrations are assumed to occur before dosing
  data_conc$out_after <- FALSE
  data_dose <-
    merge(
      data.frame(dose=TRUE,
                 time=time.dose,
                 route=route.dose,
                 duration=duration.dose,
                 iv_bolus=route.dose %in% "intravascular" & duration.dose %in% 0,
                 stringsAsFactors=FALSE),
      # Expand IV bolus dosing to have a before and after concentration
      data.frame(iv_bolus=c(FALSE, TRUE, TRUE),
                 out_after=c(FALSE, FALSE, TRUE)),
      all.x=TRUE)
  data_out <-
    data.frame(out=TRUE,
               out_after=out.after,
               out_order=1:length(time.out),
               time=time.out)
  data_all <-
    merge(merge(data_conc,
                data_dose,
                all=TRUE),
          data_out,
          all=TRUE)
  data_all$dose_event <- data_all$dose %in% TRUE
  data_all$conc_event <- !is.na(data_all$conc)
  data_all$iv_bolus <- data_all$iv_bolus %in% TRUE
  data_all$out <- data_all$out %in% TRUE
  data_all$dose_count <- cumsum(data_all$dose_event)
  mask_include_before <- data_all$dose_event & data_all$conc_event & !data_all$out_after
  data_all$dose_count_prev <- data_all$dose_count - mask_include_before
  data_all$event <- dplyr::case_when(
    data_all$dose_event & data_all$conc_event & data_all$iv_bolus & data_all$out_after~event_choices_interp.extrap.conc.dose$conc_dose_iv_bolus_after,
    data_all$dose_event & data_all$conc_event~event_choices_interp.extrap.conc.dose$conc_dose,
    data_all$dose_event & data_all$iv_bolus & data_all$out_after~event_choices_interp.extrap.conc.dose$dose_iv_bolus_after,
    data_all$dose_event~event_choices_interp.extrap.conc.dose$dose,
    data_all$conc_event~event_choices_interp.extrap.conc.dose$conc,
    data_all$out~event_choices_interp.extrap.conc.dose$output_only, # interpolation/extrapolation-only row
    TRUE~"unknown") # should never happen
  if (any(mask_unknown <- data_all$event %in% "unknown")) {
    # All events should be accounted for
    stop("Unknown event in interp.extrap.conc.dose at time(s): ",
         paste(unique(data_all$time[mask_unknown]), collapse=", ")) # nocov
  }
  # Remove "output_only" from event_before and event_after
  simple_locf <- function(x, missing_val) {
    mask_found <- !(x %in% missing_val)
    ret <- x[mask_found]
    ret[cumsum(mask_found)]
  }
  data_all$event_before <- simple_locf(c(event_choices_interp.extrap.conc.dose$none, data_all$event[-nrow(data_all)]),
                                       "output_only")
  data_all$event_after <- rev(simple_locf(rev(c(data_all$event[-1], event_choices_interp.extrap.conc.dose$none)),
                                          "output_only"))

  # Loop through the methods until all have been tested or no missing
  # values remain.
  data_all$method <- NA_character_
  for (nm in names(interp.extrap.conc.dose.select)) {
    mask <- is.na(data_all$method) &
      do.call(interp.extrap.conc.dose.select[[nm]]$select, list(x=data_all))
    if (any(mask)) {
      if ("warning" %in% names(interp.extrap.conc.dose.select[[nm]])) {
        warning(sprintf("%s: %d data points",
                        interp.extrap.conc.dose.select[[nm]]$warning,
                        sum(mask)))
        data_all$method[mask] <- nm
      } else {
        for (current_idx in which(mask)) {
          data_all$conc[current_idx] <-
            do.call(interp.extrap.conc.dose.select[[nm]]$value,
                    list(data_all=data_all,
                         current_idx=current_idx,
                         options=options,
                         ...))
          data_all$method[current_idx] <- nm
        }
      }
    }
  }
  if (any(mask_no_method <- is.na(data_all$method))) {
    # This should never happen, all eventualities should be covered
    stop("No method for imputing concentration at time(s): ",
         paste(unique(data_all$time[mask_no_method]), collapse=", ")) # nocov
  }
  # Filter to the requested time points and output
  data_out <- data_all[data_all$out,,drop=FALSE]
  data_out <- data_out[order(data_out$out_order),,drop=FALSE]
  ret <- data_out$conc
  attr(ret, "Method") <- data_out$method
  ret
}

# Dose-aware interpolation/extrapolation functions ####

# Impossible combinations ####
iecd_impossible_select <- function(x) {
  x$event_before %in% "output_only" | # Restricted in code
    x$event %in% "none" | # "none" events do not occur, they are before or after the timeline
    x$event_after %in% "output_only" | # Restricted in code
    (x$event %in% "output_only" &
       x$event_after %in% c("conc_dose_iv_bolus_after", "dose_iv_bolus_after")) |
    (x$event_before %in% c("conc_dose_iv_bolus_after", "dose_iv_bolus_after") &
       x$event %in% c("conc_dose_iv_bolus_after", "dose_iv_bolus_after")) |
    (x$event %in% c("conc_dose_iv_bolus_after", "dose_iv_bolus_after") &
       x$event_after %in% c("conc_dose_iv_bolus_after", "dose_iv_bolus_after")) |
    (x$event_before %in% c("conc_dose_iv_bolus_after", "dose_iv_bolus_after") &
       x$event %in% c("none", "output_only") &
       x$event_after %in% c("conc_dose_iv_bolus_after", "dose_iv_bolus_after"))
}
iecd_impossible_value <- function(data_all, current_idx, ...) {
  stop(
    sprintf(
      "Impossible combination requested for interp.extrap.conc.dose.  event_before: %s, event: %s, event_after: %s",
      data_all$event_before[current_idx],
      data_all$event[current_idx],
      data_all$event_after[current_idx])) # nocov
}


# Observed concentration ####
iecd_observed_select <- function(x) {
  x$event %in% c("conc_dose_iv_bolus_after", "conc_dose", "conc") &
    !do.call(interp.extrap.conc.dose.select[["Impossible combinations"]]$select, list(x=x))
}
iecd_observed_value <- function(data_all, current_idx, ...) {
  data_all$conc[current_idx]
}

# Before all events ####
iecd_before_select <- function(x) {
  x$event_before %in% "none" &
    !(x$event %in% c("conc_dose_iv_bolus_after", "conc_dose", "conc", "dose_iv_bolus_after", "none")) &
    !(x$event_after %in% "output_only" |
        (x$event %in% "output_only" &
           x$event_after %in% c("conc_dose_iv_bolus_after", "dose_iv_bolus_after"))) # because these are impossible
}
iecd_before_value <- function(data_all, current_idx, conc.origin=0, ...) {
  conc.origin
}

# Interpolation ####
iecd_interp_select <- function(x) {
  x$event_before %in% c("conc_dose_iv_bolus_after", "conc_dose", "conc") &
    x$event %in% c("output_only") &
    x$event_after %in% c("conc_dose", "conc") &
    !(x$event_before %in% "conc_dose" &
        x$event_after %in% "conc_dose")
}
iecd_interp_value <- function(data_all, current_idx, ...) {
  tmp_conc <- data_all[!is.na(data_all$conc) &
                         data_all$dose_count %in% data_all$dose_count[current_idx],]
  interpolate.conc(conc=tmp_conc$conc, time=tmp_conc$time,
                   time.out=data_all$time[current_idx],
                   check=FALSE, ...)
}

# Extrapolation ####
iecd_extrap_select <- function(x) {
  extrap_output_only <-
    x$event_before %in% c("conc_dose_iv_bolus_after", "conc") &
    x$event %in% "output_only" &
    x$event_after %in% c("dose", "none")
  extrap_dose <-
    x$event_before %in% c("conc_dose_iv_bolus_after", "conc") &
    x$event %in% "dose" &
    !(x$event_after %in% "output_only")
  extrap_output_only | extrap_dose
}
iecd_extrap_value <- function(data_all, current_idx, lambda.z, ...) {
  last_conc <- data_all[data_all$time < data_all$time[current_idx] &
                          !is.na(data_all$conc),]
  last_conc <- last_conc[nrow(last_conc),]
  if (last_conc$conc %in% 0) {
    # BLQ continues to be BLQ
    0
  } else {
    if (missing(lambda.z)) {
      lambda.z <- NA_real_
    }
    args <- list(conc=last_conc$conc[nrow(last_conc)],
                 time=last_conc$time[nrow(last_conc)],
                 time.out=data_all$time[current_idx], lambda.z=lambda.z,
                 ...)
    if (!("clast" %in% names(args))) {
      args$clast <- last_conc$conc[nrow(last_conc)]
    }
    do.call(extrapolate.conc, args)
  }
}

# Immediately after an IV bolus with a concentration next ####
iecd_iv_conc_select <- function(x) {
  !(x$event_before %in% c("conc_dose_iv_bolus_after", "dose_iv_bolus_after", "output_only")) &
    x$event %in% "dose_iv_bolus_after" &
    x$event_after %in% c("conc_dose", "conc")
}
iecd_iv_conc_value <- function(data_all, current_idx, ...) {
  tmp_conc <- data_all[data_all$conc_event &
                         (data_all$dose_count %in% data_all$dose_count[current_idx] |
                            data_all$dose_count_prev %in%  data_all$dose_count[current_idx]),]
  tmp_dose <- data_all[data_all$dose_event & 
                         data_all$dose_count %in% data_all$dose_count[current_idx],]
  pk.calc.c0(conc=tmp_conc$conc, time=tmp_conc$time,
             time.dose=tmp_dose$time,
             method=c("logslope", "c1"))
}

# Immediately after an IV bolus without a concentration next ####
iecd_iv_noconc_select <- function(x) {
  bolus_is_current <-
    x$event_before %in% c("conc", "conc_dose", "dose", "none") &
    x$event %in% "dose_iv_bolus_after" &
    x$event_after %in% c("dose", "none")
  bolus_is_previous <-
    x$event_before %in% "dose_iv_bolus_after" &
    x$event %in% "dose" &
    x$event_after %in% c("conc_dose_iv_bolus_after", "conc_dose", "dose_iv_bolus_after", "dose", "conc", "none")
  bolus_is_current | bolus_is_previous
}

# After an IV bolus with a concentration next ####
iecd_afteriv_conc_select <- function(x) {
  x$event_before %in% c("dose_iv_bolus_after") &
    x$event %in% c("output_only") &
    x$event_after %in% c("conc_dose", "conc")
}
iecd_afteriv_conc_value <- function(data_all, current_idx, ...) {
  tmp_conc <- data_all[data_all$conc_event &
                         (data_all$dose_count %in% data_all$dose_count[current_idx] |
                            data_all$dose_count_prev %in%  data_all$dose_count[current_idx]),
                       c("conc", "time")]
  tmp_dose <- data_all[data_all$dose_event & 
                         data_all$dose_count %in% data_all$dose_count[current_idx],]
  tmp_conc <- rbind(
    data.frame(
      time=tmp_dose$time,
      conc=
        pk.calc.c0(conc=tmp_conc$conc, time=tmp_conc$time,
                   time.dose=tmp_dose$time,
                   method=c("logslope", "c1"))),
    tmp_conc)
  interpolate.conc(conc=tmp_conc$conc,
                   time=tmp_conc$time,
                   time.out=data_all$time[current_idx],
                   check=FALSE, ...)
}

# After an IV bolus without a concentration next ####
iecd_afteriv_noconc_select <- function(x) {
  x$event_before %in% c("dose_iv_bolus_after") &
    x$event %in% "output_only" &
    x$event_after %in% c("dose", "none")
}

# Doses with no concentrations between ####
iecd_dose_noconc_select <- function(x) {
  dose_current <-
    x$event_before %in% c("conc_dose", "dose") &
    x$event %in% "dose" &
    !(x$event_after %in% "output_only")
  dose_around <-
    x$event_before %in% c("dose", "conc_dose") &
    x$event %in% "output_only" &
    x$event_after %in% c("dose", "conc_dose")
  dose_current | dose_around
}

# Dose as the last event in the timeline and requesting a concentration after ####
iecd_dose_last_select <- function(x) {
  x$event_before %in% c("conc_dose", "dose") &
    x$event %in% "output_only" &
    x$event_after %in% "none"
}

# Dose before, concentration after without a dose ####
iecd_dose_conc_select <- function(x) {
  x$event_before %in% "dose" &
    x$event %in% "output_only" &
    x$event_after %in% "conc"
}
iecd_dose_conc_value <- function(data_all, current_idx, ...) {
  data_tmp <- data_all[data_all$dose_event | data_all$conc_event,]
  interpolate.conc(conc=data_tmp$conc, time=data_tmp$time,
                   time.out=data_all$time[current_idx], ...,
                   check=FALSE)
}
  
interp.extrap.conc.dose.select <-
  list(
    "Impossible combinations"=
      list(
        select="iecd_impossible_select",
        value="iecd_impossible_value",
        description="The event combination cannot exist."),
    "Observed concentration"=
      list(
        select="iecd_observed_select",
        value="iecd_observed_value",
        description="Copy the input concentration at the given time to the output."),
    "Before all events"=
      list(
        select="iecd_before_select",
        value="iecd_before_value",
        description=paste("Interpolation before any events is NA or zero (0)",
                          "depending on the value of conc.origin.  conc.origin",
                          "defaults to zero which is the implicit assumption",
                          "that a complete washout occurred and there is no",
                          "endogenous source of the analyte.")),
    "Interpolation"=
      list(
        select="iecd_interp_select",
        value="iecd_interp_value",
        description=paste("With concentrations before and after and not an IV",
                          "bolus before, interpolate between observed concentrations.")),
    "Extrapolation"=
      list(
        select="iecd_extrap_select",
        value="iecd_extrap_value",
        description="Extrapolate from a concentration to a dose"),
    "Immediately after an IV bolus with a concentration next"=
      list(
        select="iecd_iv_conc_select",
        value="iecd_iv_conc_value",
        description=paste("Calculate C0 for the time immediately after an IV",
                          "bolus.  First, attempt using log slope",
                          "back-extrapolation.  If that fails, use the first",
                          "concentration after the dose as C0.")),
    "Immediately after an IV bolus without a concentration next"=
      list(
        select="iecd_iv_noconc_select",
        warning="Cannot interpolate immediately after an IV bolus without a concentration next.",
        description="Cannot calculate C0 without a concentration after an IV bolus; return NA."),
    "After an IV bolus with a concentration next"=
      list(
        select="iecd_afteriv_conc_select",
        value="iecd_afteriv_conc_value",
        description=paste("First, calculate C0 using log slope back-extrapolation",
                          "(falling back to the first post-dose concentration",
                          "if that fails).  Then, interpolate between C0 and",
                          "the first post-dose concentration.")),
    "After an IV bolus without a concentration next"=
      list(
        select="iecd_afteriv_noconc_select",
        warning="Cannot interpolate after an IV bolus without a concentration next.",
        description=paste("Between an IV bolus and anything other than a",
                          "concentration, interpolation cannot occur.  Return NA")),
    "Doses with no concentrations between"=
      list(
        select="iecd_dose_noconc_select",
        warning="Cannot interpolate between two doses or after a dose without a concentration after the first dose.",
        description="Two doses with no concentrations between them, return NA."),
    "Dose as the last event in the timeline and requesting a concentration after"=
      list(
        select="iecd_dose_last_select",
        warning="Cannot extrapolate from a dose without any concentrations after it.",
        description=paste("Cannot estimate the concentration after a dose",
                          "without concentrations after the dose, return NA.")),
    "Dose before, concentration after without a dose"=list(
      select="iecd_dose_conc_select",
      value="iecd_dose_conc_value",
      description="If the concentration at the dose is estimable, interpolate.  Otherwise, NA."))

# Internal methods for AUC integration

aucintegrate_linear <- function(conc.1, conc.2, time.1, time.2) {
  (time.2-time.1) * (conc.2+conc.1)/2
}

aucintegrate_log <- function(conc.1, conc.2, time.1, time.2) {
  (time.2-time.1) * (conc.2-conc.1)/log(conc.2/conc.1)
}

aucintegrate_inf <- function(clast, tlast, lambda.z) {
  clast/lambda.z
}

aumcintegrate_linear <- function(conc.1, conc.2, time.1, time.2) {
  (time.2-time.1) * (conc.2*time.2+conc.1*time.1)/2
}

aumcintegrate_log <- function(conc.1, conc.2, time.1, time.2) {
  ((time.2-time.1) * (conc.2*time.2-conc.1*time.1) / log(conc.2/conc.1)-
     (time.2-time.1)^2 * (conc.2-conc.1) / (log(conc.2/conc.1)^2))
}

aumcintegrate_inf <- function(conc.last, time.last, lambda.z) {
  (conc.last*time.last/lambda.z) + conc.last/(lambda.z^2)
}

#' Interpolate or extrapolate concentrations using the provided method
#'
#' @param conc_1,conc_2 The concentration at time1 and time2
#' @param time_1,time_2 The time value associated with conc1 and conc2
#' @param time_out Time when interpolation is requested
#' @param tlast The time of the last concentration above the lower limit of
#'   quantification (LOQ)
#' @returns The interpolated or extrapolated value using the correct method
#' @keywords Internal
#' @name interp_extrap_conc_method
NULL

#' @rdname interp_extrap_conc_method
interpolate_conc_linear <- function(conc_1, conc_2, time_1, time_2, time_out) {
  conc_1+(time_out-time_1)/(time_2-time_1)*(conc_2-conc_1)
}

#' @rdname interp_extrap_conc_method
interpolate_conc_log <- function(conc_1, conc_2, time_1, time_2, time_out) {
  exp(
    log(conc_1)+
      (time_out-time_1)/(time_2-time_1)*(log(conc_2)-log(conc_1))
  )
}

#' @inheritParams assert_lambdaz
#' @param clast The concentration at the last time above the lower LOQ
#' @rdname interp_extrap_conc_method
extrapolate_conc_lambdaz <- function(clast, lambda.z, tlast, time_out) {
  clast*exp(-lambda.z*(time_out - tlast))
}

#' Choose how to interpolate, extrapolate, or integrate data in each
#' concentration interval
#'
#' @inheritParams assert_conc_time
#' @inheritParams PKNCA.choose.option
#' @inheritParams assert_aucmethod
#' @param auc.type The type of AUC to compute.  Choices are 'AUCinf', 'AUClast',
#'   and 'AUCall'.
#' @param tlast Time of last concentration above the limit of quantification
#'   (will be calculated, if not provided)
#' @keywords Internal
#' @returns A character vector of methods for interpolation/extrapolation
#'   methods that is the same length as `conc` which indicates how to
#'   interpolate/integrate between each of the concentrations (all but the last
#'   value in the vector) and how to extrapolate after `tlast` (the last item in
#'   the vector).  Possible values in the vector are: 'zero', 'linear', 'log',
#'   and 'extrap_log'
choose_interval_method <- function(conc, time, tlast, method, auc.type, options) {
  # Input checking
  stopifnot(is.numeric(conc))
  stopifnot(is.numeric(time))
  stopifnot(!any(is.na(time)))
  stopifnot(!any(is.na(conc)))
  stopifnot(length(conc) == length(time))
  assert_aucmethod(method)
  stopifnot(length(auc.type) == 1)
  stopifnot(auc.type %in% c("AUCinf", "AUClast", "AUCall"))

  if (missing(tlast)) {
    tlast <- pk.calc.tlast(conc, time, check=FALSE)
  } else {
    stopifnot(is.numeric(tlast))
    stopifnot(length(tlast) == 1)
  }

  # Where is tlast in the data?
  idx_tlast <- which(time == tlast)

  ret <- rep(NA_character_, length(conc))
  # Handle all interpolation
  idx_1 <- seq(1, length(conc) - 1)
  idx_2 <- idx_1 + 1
  mask_zero <- conc[idx_1] == 0 & conc[idx_2] == 0
  if (all(conc %in% 0)) {
    ret[] <- "zero"
    # short circuit other options
    return(ret)
  } else if (method == "linear") {
    ret[seq_len(idx_tlast - 1)] <- "linear"
  } else if (method == "lin up/log down") {
    mask_down <- conc[idx_2] < conc[idx_1] & conc[idx_2] != 0
    mask_up <- !(mask_down | mask_zero)
    ret[c(mask_down, FALSE)] <- "log"
    ret[c(mask_up, FALSE)] <- "linear"
  } else if (method == "lin-log") {
    # We only need tmax for lin-log, and it is always recalculated to prevent
    # potential issues where tmax may be relative to the start of an interval vs
    # relative to absolute time.
    tmax <- pk.calc.tmax(conc = conc, time = time, options = options)
    mask_pre_tmax <- time[idx_2] <= tmax
    mask_post_tmax <- !mask_pre_tmax
    mask_zero_start_end <- conc[idx_1] %in% 0 | conc[idx_2] %in% 0
    mask_linear <- mask_pre_tmax | mask_zero_start_end
    mask_log <- mask_post_tmax & !mask_zero_start_end
    ret[c(mask_linear, FALSE)] <- "linear"
    ret[c(mask_log, FALSE)] <- "log"
  } else {
    stop("Unknown integration method, please report a bug: ", method) # nocov
  }
  ret[c(mask_zero, FALSE)] <- "zero"
  # What happens after tlast?
  if (idx_tlast < (length(ret) - 1)) {
    ret[seq(idx_tlast + 1, length(ret) - 1)] <- "zero"
  }
  # what happens at tlast?
  if (idx_tlast < length(conc)) {
    if (auc.type == "AUCall") {
      ret[idx_tlast] <- "linear"
    } else {
      ret[idx_tlast] <- "zero"
    }
  }

  # Handle extrapolation
  if (auc.type == "AUCinf") {
    ret[length(ret)] <- "extrap_log"
  } else {
    # Extrapolation is zero, except for AUCinf
    ret[length(ret)] <- "zero"
  }
  ret
}

#' Support function for AUC integration
#'
#' @inheritParams choose_interval_method
#' @inheritParams pk.calc.auxc
#' @param clast The last concentration above the limit of quantification
#' @param interval_method The method for integrating each interval of `conc`
#' @keywords Internal
auc_integrate <- function(conc, time, clast, tlast, lambda.z, interval_method, fun_linear, fun_log, fun_inf) {
  assert_lambdaz(lambda.z = lambda.z)
  interval_method_within <- interval_method[-length(interval_method)]
  interval_method_extrap <- interval_method[length(interval_method)]
  idx_1 <- seq_len(length(conc) - 1)
  idx_1_linear <- idx_1[interval_method_within == "linear"]
  idx_1_log <- idx_1[interval_method_within == "log"]

  ret <-
    c(
      fun_linear(conc[idx_1_linear], conc[idx_1_linear + 1],
                 time[idx_1_linear], time[idx_1_linear + 1]),
      fun_log(conc[idx_1_log], conc[idx_1_log + 1],
              time[idx_1_log], time[idx_1_log + 1])
    )

  if (interval_method_extrap %in% "extrap_log") {
    # Whether AUCinf,obs or AUCinf,pred is calculated depends on if clast,obs
    # or clast,pred is passed in.
    ret[length(ret)+1] <- fun_inf(clast, tlast, lambda.z)
  } else if (interval_method_extrap != "zero") {
    stop("Invalid interval_method_extrap, please report a bug: ", interval_method_extrap) # nocov
  }
  ret <- sum(ret)
  ret
}

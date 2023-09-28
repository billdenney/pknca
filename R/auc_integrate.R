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

#' Choose how to interpolate, extrapolate, or integrate data in each
#' concentration interval
#'
#' @param conc Concentration measured
#' @param time Time of concentration measurement (must be monotonically
#'   increasing and the same length as the concentration data)
#' @param method The method for integration (either 'lin up/log down' or
#'   'linear')
#' @param auc.type The type of AUC to compute.  Choices are 'AUCinf', 'AUClast',
#'   and 'AUCall'.
#' @param tlast Time of last concentration above the limit of quantification
#'   (will be calculated, if not provided)
#' @keywords Internal
#' @return A character vector of methods for interpolation/extrapolation methods
#'   that is the same length as `conc` which indicates how to
#'   interpolate/integrate between each of the concentrations (all but the last
#'   value in the vector) and how to extrapolate after `tlast` (the last item in
#'   the vector).  Possible values in the vector are: 'zero', 'linear', 'log',
#'   and 'extrap_log'
choose_interval_method <- function(conc, time, tlast, method, auc.type) {
  # Input checking
  stopifnot(is.numeric(conc))
  stopifnot(is.numeric(time))
  stopifnot(!any(is.na(time)))
  stopifnot(!any(is.na(conc)))
  stopifnot(length(conc) == length(time))
  stopifnot(length(method) == 1)
  stopifnot(method %in% c("lin up/log down", "linear"))
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

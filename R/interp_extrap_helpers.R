#' Choose a method for calculation in the interval between concentrations
#'
#' This function should be used for any interpolation/extrapolation function. It
#' will standardize the method of choosing which method to use for interpolation
#' and extrapolation.
#'
#' @param conc A vector of concentrations (\code{NA} values are not allowed)
#' @param time A vector of times (\code{NA} values are not allowed)
#' @param interp_method Method to use for interpolation between time points
#' @param extrap_method Method to use for extrapolation after the last time
#'   point above (an AUC calculation method)
#' @param tmax Time of maximum concentration
#' @return A character vector of extrapolation methods to use between each
#'   \code{conc} and after the last \code{conc}.  Values will be one or more of
#'   "linear" (use linear interpolation), "log" (use log interpolation), "zero"
#'   (the value is zero), and the last value may be "clastpred", "clastobs", or
#'   "zero" indicating extrapolation from tlast using lambda.z and clast,pred or
#'   clast,obs, or zero.
#' @keywords Internal
#' @examples
#' PKNCA:::choose_interp_extrap_method(
#'   conc=c(1, 2, 4, 2, 1, 0, 0),
#'   time=0:6,
#'   interp_method="lin up/log down",
#'   extrap_method="aucinf.obs"
#' )
choose_interp_extrap_method <- function(conc, time, interp_method, extrap_method, tmax) {
  checkmate::assert_numeric(x=conc, any.missing=FALSE, min.len=1)
  checkmate::assert_numeric(x=time, any.missing=FALSE, len=length(conc), unique=TRUE, sorted=TRUE)
  checkmate::assert_choice(x=interp_method, choices=c("linear", "lin up/log down", "lin/log", "log"))
  checkmate::assert_choice(x=extrap_method, choices=c("aucinf.obs", "aucinf.pred", "auclast", "aucall"))
  
  if (length(conc) == 1) {
    # only extrapolate (consider adding a warning, but that may be too noisy)
    has_conc_zero_after <- conc == 0
    ret_interp <- character(0)
    after_tlast <- FALSE
  } else if (length(conc) > 1) {
    # Evaluate interpolation options
    conc_1 <- conc[-length(conc)]
    conc_2 <- conc[-1]

    has_conc_zero_after <- conc_2 <= 0

    tlast <- pk.calc.tlast(conc, time)
    before_tlast <- time[-1] <= tlast
    after_tlast <- !before_tlast
  
    ret_interp <- rep(NA_character_, length(conc) - 1)
    ret_interp[after_tlast] <- "zero"
    if (all(has_conc_zero_after)) {
      ret_interp[] <- "zero"
    } else if (interp_method %in% c("linear", "log")) {
      ret_interp[c(!has_conc_zero_after & before_tlast)] <- interp_method
      ret_interp[c(has_conc_zero_after & before_tlast)] <- "linear"
    } else if (interp_method == "lin up/log down") {
      has_conc_increasing <- !has_conc_zero_after & (conc_1 <= conc_2)
      use_linear <- has_conc_zero_after | has_conc_increasing
      ret_interp[c(use_linear & before_tlast)] <- "linear"
      ret_interp[c(!use_linear & before_tlast)] <- "log"
    } else if (interp_method == "lin/log") {
      has_conc_zero_before <- conc_1 <= 0
      before_tmax <- time[-1] <= tmax
      use_linear <- has_conc_zero_after | has_conc_zero_before | before_tmax
      ret_interp[c(use_linear & before_tlast)] <- "linear"
      ret_interp[c(!use_linear & before_tlast)] <- "log"
    } else {
      # I don't think you can get here with the assert_choice above
      stop("Please report a bug: Invalid interp_method: ", interp_method) # nocov
    }
  }
  # extrapolation
  if (all(has_conc_zero_after)) {
    # everything is zero
    ret_extrap <- "zero"
  } else {
    # All extrapolation methods happen from tlast to infinity (or tlast to the
    # next point), so everything else is zero
    if (extrap_method == "aucinf.obs") {
      ret_extrap <- "clastobs"
    } else if (extrap_method == "aucinf.pred") {
      ret_extrap <- "clastpred"
    } else if (extrap_method == "auclast") {
      ret_extrap <- "zero"
    } else if (extrap_method == "aucall") {
      at_tlast <- which(after_tlast)[1]
      if (!is.na(at_tlast)) {
        # If at_tlast is NA, then the last measurement is tlast.
        ret_interp[at_tlast] <- "linear"
      }
      ret_extrap <- "zero"
    } else {
      # I don't think you can get here with the assert_choice above
      stop("Please report a bug: Invalid extrap_method: ", extrap_method) # nocov
    }
  }
  c(ret_interp, ret_extrap)
}

#' Interpolate or extrapolate concentrations using the provided method
#'
#' @param conc_1,conc_2 The concentration at time1 and time2
#' @param time_1,time_2 The time value associated with conc1 and conc2
#' @param time_out Time when interpolation is requested
#' @param tlast The time of the last concentration above the lower limit of
#'   quantification (LOQ)
#' @return The interpolated or extrapolated value using the correct method
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

#' @param clast The concentration at the last time above the lower LOQ
#' @param lambda.z The elimination rate
#' @rdname interp_extrap_conc_method
extrapolate_conc_lambdaz <- function(clast, lambda.z, tlast, time_out) {
  clast*exp(-lambda.z*(time_out - tlast))
}

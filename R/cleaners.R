#' Handle NA values in the concentration measurements as requested by the user.
#'
#' NA concentrations (and their associated times) will be removed then the BLQ
#' values in the middle
#'
#' @inheritParams assert_conc_time
#' @inheritParams PKNCA.choose.option
#' @param \dots Additional items to add to the data frame
#' @param conc.na How to handle NA concentrations?  Either 'drop' or a number to
#'   impute.
#' @param check Run [assert_conc_time()]?
#' @returns The concentration and time measurements (data frame) filtered
#'   and cleaned as requested relative to NA in the concentration.
#' @family Data cleaners
#' @export
clean.conc.na <- function(conc, time, ...,
                          options=list(),
                          conc.na=NULL,
                          check=TRUE) {
  conc.na <- PKNCA.choose.option(name="conc.na", value=conc.na, options=options)
  if (check)
    assert_conc_time(conc, time)
  # Prep it as a data frame
  ret <- data.frame(conc, time, ..., stringsAsFactors=FALSE)
  if (conc.na %in% "drop") {
    # If it is set to "drop" then omit the NA concentrations
    ret <- ret[!is.na(conc),]
  } else if (is.numeric(conc.na)) {
    ret$conc[is.na(conc)] <- conc.na
  } else {
    # This case should already have been captured by the PKNCA.options
    # call above.
    stop("Unknown how to handle conc.na") # nocov
  }
  ret
}

#' Handle BLQ values in the concentration measurements as requested by the user.
#'
#' @inheritParams assert_conc_time
#' @inheritParams PKNCA.choose.option
#' @param \dots Additional arguments passed to clean.conc.na
#' @param conc.blq How to handle a BLQ value that is between above LOQ values?
#'   See details for description.
#' @param conc.na How to handle NA concentrations.  (See [clean.conc.na()])
#' @param check Run [assert_conc_time()]?
#' @returns The concentration and time measurements (data frame) filtered and
#'   cleaned as requested relative to BLQ in the middle.
#'
#' @details `NA` concentrations (and their associated times) will be handled as
#'   described in [clean.conc.na()] before working with the BLQ values.  The
#'   method for handling NA concentrations can affect the output of which points
#'   are considered BLQ and which are considered "middle".  Values are
#'   considered BLQ if they are 0.
#'
#'   `conc.blq` can be set either a scalar indicating what should be done for
#'   all BLQ values or a list with elements either named "first", "middle" and "last" or "before.tmax" and "after.tmax"
#'   each set to a scalar.
#'
#' The meaning of each of the list elements is:
#' \describe{
#'   \item{first}{Values up to the first non-BLQ value.  Note
#'     that if all values are BLQ, this includes all values.}
#'   \item{middle}{Values that are BLQ between the first and last
#'     non-BLQ values.}
#'   \item{last}{Values that are BLQ after the last non-BLQ value}
#'   \item{before.tmax}{Values that are BLQ before the time at first maximum concentration}
#'   \item{after.tmax}{Values that are BLQ after the time at first maximum concentration}
#' }
#'
#' The valid settings for each are:
#' \describe{
#'   \item{"drop"}{Drop the BLQ values}
#'   \item{"keep"}{Keep the BLQ values}
#'   \item{a number}{Set the BLQ values to that number}
#' }
#'
#' @family Data cleaners
#' @export
clean.conc.blq <- function(conc, time,
                           ...,
                           options=list(),
                           conc.blq=NULL,
                           conc.na=NULL,
                           check=TRUE) {
  conc.blq <- PKNCA.choose.option(name="conc.blq", value=conc.blq, options=options)
  conc.na <- PKNCA.choose.option(name="conc.na", value=conc.na, options=options)
  if (check) {
    assert_conc_time(conc, time)
  }
  # Handle NA concentrations and make the data frame
  ret <- clean.conc.na(conc, time, ..., conc.na=conc.na, check=FALSE)
  # If all data has been excluded, then don't do anything
  if (nrow(ret) > 0) {
    tfirst <- pk.calc.tfirst(ret$conc, ret$time, check=FALSE)
    tlast <- pk.calc.tlast(ret$conc, ret$time, check=FALSE)
    tmax <- pk.calc.tmax(ret$conc, ret$time, check=FALSE)

    # If all measurements are BLQ
    if (all(ret$conc == 0)){
      # Apply "first" BLQ rule to everything for tfirst/tlast
      tfirst <- max(ret$time)
      tlast <- tfirst + 1

      # Apply "before.tmax" BLQ rule to everything for tmax
      tmax <- max(ret$time)
    }

    # Depending on the specified argument perform the corresponding action
    for (i in seq_len(length(conc.blq))) {
      # Set the mask to apply the rule to
      time_type <- names(conc.blq)[i]
      if (is.null(time_type) & length(conc.blq) == 1) {
        mask <- ret$conc %in% 0
      } else if (time_type == "first") {
        mask <- ret$time <= tfirst & ret$conc %in% 0
      } else if (time_type == "middle") {
        mask <- tfirst < ret$time & ret$time < tlast & ret$conc %in% 0
      } else if (time_type == "last") {
        mask <- tlast <= ret$time & ret$conc %in% 0
      } else if (time_type == "before.tmax") {
        mask <- ret$time < tmax & ret$conc %in% 0
      } else if (time_type == "after.tmax") {
        mask <- tmax <= ret$time & ret$conc %in% 0
      } else {
        stop("There is a bug in cleaning the conc.blq with position names") # nocov
      }
      # Choose the rule to apply
      this_rule <- unname(conc.blq)[[i]]

      if (this_rule %in% "keep") {
        # Do nothing
      } else if (this_rule %in% "drop") {
        ret <- ret[!mask,]
      } else if (is.numeric(this_rule)) {
        ret$conc[mask] <- this_rule
      } else {
        # This case should already have been captured by the PKNCA.options
        # call above.
        stop(sprintf("Unknown how to handle conc.blq rule %s", # nocov
                     as.character(this_rule)))                 # nocov
      }
    }
  }
  ret
}

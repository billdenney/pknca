#' Handle NA values in the concentration measurements as requested by 
#' the user.
#' 
#' NA concentrations (and their associated times) will be removed then
#' the BLQ values in the middle
#' 
#' @param conc Measured concentrations
#' @param time Time of the concentration measurement
#' @param \dots Additional items to add to the data frame
#' @param options List of changes to the default 
#'   \code{\link{PKNCA.options}} for calculations.
#' @param conc.na How to handle NA concentrations?  Either 'drop' or a 
#'   number to impute.
#' @param check Run \code{\link{check.conc.time}}?
#' @return The concentration and time measurements (data frame) filtered
#'   and cleaned as requested relative to NA in the concentration.
#' @family Data cleaners
#' @export
clean.conc.na <- function(conc, time, ...,
                          options=list(),
                          conc.na=NULL,
                          check=TRUE) {
  conc.na <- PKNCA.choose.option(name="conc.na", value=conc.na, options=options)
  if (check)
    check.conc.time(conc, time)
  ## Prep it as a data frame
  ret <- data.frame(conc, time, ..., stringsAsFactors=FALSE)
  if (conc.na %in% "drop") {
    ## If it is set to "drop" then omit the NA concentrations
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

#' Handle BLQ values in the concentration measurements as requested by 
#' the user.
#' 
#' @param conc Measured concentrations
#' @param time Time of the concentration measurement
#' @param \dots Additional arguments passed to clean.conc.na
#' @param options List of changes to the default 
#'   \code{\link{PKNCA.options}} for calculations.
#' @param conc.blq How to handle a BLQ value that is between above LOQ 
#'   values?  See details for description.
#' @param conc.na How to handle NA concentrations.  (See 
#'   \code{\link{clean.conc.na}})
#' @param check Run \code{\link{check.conc.time}}?
#' @return The concentration and time measurements (data frame) filtered
#'   and cleaned as requested relative to BLQ in the middle.
#'
#' @details NA concentrations (and their associated times) will be 
#'   handled as described in \code{\link{clean.conc.na}} before working 
#'   with the BLQ values.  The method for handling NA concentrations can
#'   affect the output of which points are considered BLQ and which are 
#'   considered "middle".  Values are considered BLQ if they are 0.
#'
#' \code{conc.blq} can be set either a scalar indicating what
#' should be done for all BLQ values or a list with elements named
#' "first", "middle", and "last" each set to a scalar.
#'
#' The meaning of each of the list elements is:
#' \describe{
#'   \item{first}{Values up to the first non-BLQ value.  Note
#'     that if all values are BLQ, this includes all values.}
#'   \item{middle}{Values that are BLQ between the first and last
#'     non-BLQ values.}
#'   \item{last}{Values that are BLQ after the last non-BLQ value}
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
  if (check)
    check.conc.time(conc, time)
  ## Handle NA concentrations and make the data frame
  ret <- clean.conc.na(conc, time, ..., conc.na=conc.na, check=FALSE)
  ## If all data has been excluded, then don't do anything
  if (nrow(ret) > 0) {
    tfirst <- pk.calc.tfirst(ret$conc, ret$time, check=FALSE)
    if (is.na(tfirst)) {
      ## All measurements are BLQ; so apply the "first" BLQ rule to
      ## everyting.
      tfirst <- max(ret$time)
      tlast <- tfirst + 1
    } else {
      ## There is at least one above LOQ concentration
      tlast <- pk.calc.tlast(ret$conc, ret$time, check=FALSE)
    }
    ## For each of the first, middle, and last, do the right thing to
    ## the values in that set.
    for (n in c("first", "middle", "last")) {
      ## Set the mask to apply the rule to
      if (n == "first") {
        mask <- (ret$time <= tfirst &
                   ret$conc %in% 0)
      } else if (n == "middle") {
        mask <- (tfirst < ret$time &
                   ret$time < tlast &
                     ret$conc %in% 0)
      } else if (n == "last") {
        mask <- (tlast <= ret$time &
                   ret$conc %in% 0)
      } else {
        stop("There is a bug in cleaning the conc.blq with position names") # nocov
      }
      ## Choose the rule to apply
      this_rule <-
        if (is.list(conc.blq)) {
          conc.blq[[n]]
        } else {
          conc.blq
        }
      if (this_rule %in% "keep") {
        ## Do nothing
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

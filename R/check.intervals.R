#' Check the formatting of a calculation interval specification data
#' frame.
#'
#' Calculation interval specifications are data frames defining what
#' calculations will be required and summarized from all time
#' intervals.  Note: parameters which are not requested may be
#' calculated if it is required for (or computed at the same time as)
#' a requested parameter.
#' 
#' The data frame columns (with variable type in parentheses) are:
#' \describe{
#'   \item{\code{start}}{The starting time (numeric)}
#'   \item{\code{end}}{The ending time (numeric including Inf)}
#'   \item{\code{aucinf}}{Compute AUCinf (logical)}
#'   \item{\code{auclast}}{Compute AUClast (logical)}
#'   \item{\code{aucall}}{Compute AUCall (logical)}
#'   \item{\code{aumcinf}}{Compute AUMCinf (logical)}
#'   \item{\code{aumclast}}{Compute AUMClast (logical)}
#'   \item{\code{aumcall}}{Compute AUMCall (logical)}
#'   \item{\code{tfirst}}{The time of the first concentration above the limit
#'     of quantification (logical)}
#'   \item{\code{tmax}}{The time of observed maximum concentration (logical)}
#'   \item{\code{tlast}}{The time of the last concentration above the limit
#'     of quantification (logical)}
#'   \item{\code{cmin}}{The observed minimum concentration during the interval
#'     (logical)}
#'   \item{\code{cmax}}{The observed maximum concentration (logical)}
#'   \item{\code{clast.obs}}{The observed last concentration (logical)}
#'   \item{\code{clast.pred}}{The concentration at \code{tlast} predicted by
#'     the half life (logical)}
#'   \item{\code{half.life}}{The half-life (logical)}
#'   \item{\code{thalf.eff}}{The effective half-life (logical)}
#'   \item{\code{aucpext}}{The percent of the AUCinf that is extrapolated after
#'     the AUClast (logical)}
#'   \item{\code{cl}}{The clearance (force: 'force' indicates that
#'     clearance should be calculated even if it is a multiple-dose study
#'     and the drug has not reached steady-state.)}
#'   \item{\code{mrt}}{The mean residence time (logical)}
#'   \item{\code{vz}}{Terminal volume of distribution (logical)}
#'   \item{\code{vss}}{Steady-state volume of distribution (logical)}
#' }
#'
#' The variable types for each column are:
#' \describe{
#'   \item{logical}{A logical variable.}
#'   \item{numeric}{A numeric (non-factor) column}
#'   \item{force}{logical or the text \code{'force'}.  \code{'force'}
#'     indicates that checking if the calculation is appropriate should be
#'     skipped.}
#'   \item{character or factor}{The text suggested as either a character or
#'     a factor}
#' }
#' 
#' \code{start} and \code{end} time must always be given, and the
#' \code{start} must be before the \code{end}.
#'
#' @param x The data frame specifying what to calculate during each
#' time interval
#' @return x The potentially updated data frame with the interval
#' calculation specification.
#'
#' @seealso \code{\link{check.interval.deps}}
#' @export
check.interval.specification <- function(x) {
  if (!is.data.frame(x)) {
    ## Just a warning and let as.data.frame make it an error if
    ## it can't be coerced.
    warning("AUC specification must be a data.frame")
    x <- as.data.frame(x, stringsAsFactors=FALSE)
  }
  if (nrow(x) == 0)
    stop("interval specification has no rows")
  ## Confirm that the minimal columns (start and end) exist
  if (length(missing.required.cols <- setdiff(c("start", "end"), names(x))) > 0)
    stop(sprintf("Column(s) %s missing from interval specification",
                 paste0("'", missing.required.cols, "'",
                        collapse=", ")))
  interval.cols <- get.interval.cols()
  ## Check the edit of each column
  for (n in names(interval.cols))
    if (!(n %in% names(x))) {
      if (is.vector(interval.cols[[n]]$values)) {
        ## Set missing columns to the default value
        x[,n] <- interval.cols[[n]]$values[1]
      } else {
        stop("Cannot assign default value for interval column", n)
      }
    } else {
      ## Confirm the edits of the given columns
      if (is.vector(interval.cols[[n]]$values)) {
        if (!all(x[,n] %in% interval.cols[[n]]$values))
          stop(sprintf("Invalid value(s) in column %s:", n),
               paste(unique(setdiff(x[,n], interval.cols[[n]]$values)),
                     collapse=", "))
      } else if (is.function(interval.cols[[n]]$values)) {
        if (is.factor(x[,n]))
          stop(sprintf("Interval column '%s' should not be a factor", n))
        interval.cols[[n]]$values(x[,n])
      } else {
        stop("Invalid 'values' for column specification", n)
      }
    }
  ## Now check specific columns
  ## ##############################
  ## start and end
  if (any(x$start %in% NA))
    stop("AUC specification may not have NA for the starting time")
  if (any(x$end %in% NA))
    stop("AUC specification may not have NA for the end time")
  if (any(is.infinite(x$start)))
    stop("start may not be infinite")
  if (any(x$start >= x$end))
    stop("start must be < end")
  ## Confirm that something is being calculated for each interval (and
  ## warn if not)
  mask.calculated <- rep(FALSE, nrow(x))
  for (n in setdiff(names(interval.cols), c("start", "end")))
    mask.calculated <-
      (mask.calculated |
       !(x[,n] %in% c(NA, FALSE)))
  if (any(!mask.calculated))
    warning("Nothing to be calculated in interval specification number(s): ",
            paste((1:nrow(x))[!mask.calculated], collapse=", "))
  ## Put the columns in the right order and return the checked data
  ## frame
  x[,c(names(interval.cols),
       setdiff(names(x), names(interval.cols)))]
}

#' Take in a single row of an interval specification and return that
#' row updated with any additional calculations that must be done to
#' fulfil all dependencies.
#'
#' @param x A data frame with one or morw rows of the PKNCA interval
#' @return The interval specification with additional calculations
#' added where requested outputs require them.
#' @seealso \code{\link{check.interval.specification}}
check.interval.deps <- function(x) {
  ## Ensure that the input is a valid interval specification
  ret <- check.interval.specification(x)
  colspec <- get.interval.cols()
  for (n in names(colspec)) {
    if (is.logical(ret[,n])) {
      ## This is a calculation to complete, otherwise it's something
      ## informative but not caluclated.
      mask.calculated <- ret[,n]
      for (deps in colspec[[n]]$depends)
        ret[,deps] <- mask.calculated | ret[,deps]
    }
  }
  ret
}

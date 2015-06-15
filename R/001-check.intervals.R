#' Check the column edits for the given column.  This is not exported.
interval_checkers <- list(
  logical=make.logical,
  numeric=function(x) {
    if ((!is.factor(x)) & is.numeric(x))
      return(x)
    stop("Must be numeric and not a factor")
  },
  character=function(x) {
    if (is.factor(x))
      x <- as.character(x)
    if (!is.character(x))
      stop("Must be a character or a factor")
    x
  },
  force=function(x) {
    if (is.factor(x))
      x <- as.character(x)
    ret <- rep(NA, length(x))
    if (any(mask.force <- tolower(x) %in% 'force'))
      ret[mask.force] <- 'force'
    if (any(!mask.force))
      ret[!mask.force] <- make.logical(x[!mask.force])
    ret
  })
interval_defaults <- list(
  logical=FALSE,
  numeric=NA,
  character=as.character(NA),
  force=FALSE)

interval_edits <- list(
  start="numeric",
  end="numeric",
  auc.type="character",
  half.life="logical",
  tfirst="logical",
  tmax="logical",
  tlast="logical",
  cmin="logical",
  cmax="logical",
  clast.obs="logical",
  clast.pred="logical",
  thalf.eff="logical",
  aucpext="logical",
  cl="force",
  mrt="logical",
  vz="logical",
  vss="logical")

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
#'   \item{\code{auc.type}}{What type of AUC should be computed: 'AUCinf',
#'     'AUClast', 'AUCall', or NA (no AUC computed during the current
#'     interval)}
#'   \item{\code{half.life}}{The half-life (logical)}
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
#'   \item{logical}{A logical variable, as interpreted with the
#'     \code{\link{make.logical}} function.}
#'   \item{numeric}{A numeric (non-factor) column}
#'   \item{force}{logical or the text \code{'force'}.  \code{'force'}
#'     indicates that checking if the calculation is appropriate should be
#'     skipped.}
#'   \item{character or factor}{The text suggested as either a character or
#'     a factor}
#' }
#' 
#' They are interpreted with the following rules:
#'
#' \itemize{
#'   \item The \code{start} time must always be given.
#'   \item The \code{start} must be before the \code{end} (if \code{end} is
#'         given)
#'   \item If the \code{end} time is given, then the other
#'         specifications (i.e. \code{auc.type} and
#'         \code{half.life} will only be done in that interval.
#'   \item If the \code{end} time is NA, then
#'   \itemize{
#'     \item No AUCs will be calculated.
#'     \item It is an error to set \code{auc.type} to anything other than NA.
#'     \item It is an error to set \code{half.life} to FALSE.
#'     \item \code{half.life} is computed from the \code{start} to the last
#'           available measurement.
#'   }
#' }
#'
#' @param x The data frame specifying what to calculate during each
#' time interval
#' @return x The potentially updated data frame with the interval
#' calculation specification.
#' @export
check.interval.specification <- function(x) {
  if (!is.data.frame(x)) {
    ## Just a warning and let as.data.frame make it an error if
    ## it can't be coerced.
    warning("AUC specification must be a data.frame")
    x <- as.data.frame(x, stringsAsFactors=FALSE)
  }
  if (nrow(x) == 0) {
    ## Return it as is-- nothing is requested
    return(x)
  }
  ## Check for extra columns
  if (length(new.names <- setdiff(names(x), names(interval_edits))) > 0) {
    warning("Some names in the interval specification are given though not required: ",
            paste(new.names, collapse=", "))
  }
  ## Check the edit of each column
  for (n in names(interval_edits))
    if (!(n %in% names(x))) {
      ## Set missing columns to the default value
      x[,n] <- interval_defaults[[interval_edit[[n]]]]
    } else {
      ## Confirm the edits of the given columns
      x[,n] <- interval_checkers[[interval_edits[[n]]]](x[,n])
    }
  ## Now check specific columns
  ## ##############################
  ## start and end
  if (any(x$start %in% NA))
    stop("start may not be NA")
  if (any(x$end %in% NA))
    stop("end may not be NA")
  if (any(is.infinite(x$start)))
    stop("start may not be infinite")
  if (any(x$start >= x$end))
    stop("start must be < end")
  ## auc.type
  if (!all(tolower(x$auc.type) %in% c("aucinf", "auclast", "aucall", NA)))
    stop("auc.type must be one of 'aucinf', 'auclast', 'aucall', or NA")
  ## 
  ## If end, last, all, or half.life is missing, add it as all NA.
  ## Inappropriate combinations will be checked later.
  added.cols <- list(end=NA, auc.type="AUClast", half.life=FALSE)
  for (new.col in names(added.cols))
    if (!(new.col %in% names(x))) {
      warning(sprintf("AUC specification column '%s' is missing.  Adding with all %s",
                      new.col, as.character(added.cols[[new.col]])))
      x[,new.col] <- added.cols[[new.col]]
    }
  required.cols <- c("start", "end", "auc.type", "half.life")
  if (!all(required.cols %in% names(x)))
    stop(paste("AUC specification must have columns for",
               paste(required.cols, collapse=", "),
               "Column(s) missing:",
               paste(setdiff(required.cols, names(x)))))
  ## If there are additional columns, remove them
  addl.cols <- setdiff(names(x), required.cols)
  if (length(addl.cols) > 0) {
    message(paste("AUC specification columns to specify the group(s):",
                  paste(addl.cols, collapse=", ")))
  }
  ## Ensure that all columns have the right edit(s)
  if (!is.logical(x$half.life)) {
    warning("AUC specification for 'half.life' must be a logical vector, attempting conversion")
    x$half.life <- check.conversion(x$half.life, as.logical)
  }
  if (!is.numeric(x$start)) {
    warning("AUC specification for 'start' time must be numeric, attempting conversion")
    x$start <- check.conversion(x$start, as.numeric)
  }
  if (!all(is.na(x$end) | is.numeric(x$end))) {
    warning("AUC specification for 'end' time must be numeric or NA, attempting conversion")
    x$end <- check.conversion(x$end, as.numeric)
  }
  ## Ensure that only valid combinations are provided
  if (any(is.na(x$start)))
    stop("AUC specification may not have NA for the starting time")
  if (!(is.character(x$auc.type) | is.factor(x$auc.type) | all(is.na(x$auc.type)))) {
    stop("auc.type must be either a character or factor (or NA)")
  }
  if (any(is.na(x$end) & !is.na(x$auc.type)))
    stop("AUC specification may not have NA for the end time and request an auc.type")
  if (any(is.na(x$end) & !x$half.life))
    stop("AUC specification may not have NA for the end time and not request half.life")
  if (any(!is.na(x$end) & x$end <= x$start))
    stop("AUC specification end must be after the start when end is given")
  x
}

#' Make a into a logical variable interpreting more than just
#' TRUE/FALSE.
#'
#' This function will take a vector (logical, numeric, character, or
#' factor) and convert it into a logical vector.
#'
#' The values in the vector will be interpreted as follows:
#'
#' \itemize{
#'   \item TRUE
#'   \itemize{
#'     \item TRUE (as a logical value)
#'     \item "TRUE" as a character string or factor level
#'     \item "T" as a character string or factor level
#'     \item "YES" as a character string or factor level
#'     \item "Y" as a character string or factor level
#'     \item nonzero as a numeric (non-factor) value
#'   }
#'   \item FALSE
#'   \itemize{
#'     \item FALSE (as a logical value)
#'     \item "FALSE" as a character string or factor level
#'     \item "F" as a character string or factor level
#'     \item "NO" as a character string or factor level
#'     \item "N" as a character string or factor level
#'     \item zero as a numeric (non-factor) value
#'   }
#' }
#'
#' Factors are converted to characters.  Characters and factor levels
#' are compared case-insensitively.  Character strings that are not in
#' the specified list are considered to be NAs.  Logical vectors are
#' returned with \code{NA} values converted to the \code{na.value}.
#'
#' @param x The vector to convert to a logical value
#' @param na.value What should NA be converted into?
#' @return A logical vector
#' @export
make.logical <- function(x, na.value=FALSE) {
  ## Handle factors
  if (is.factor(x))
    x <- as.character(x)
  if (is.logical(x)) {
    ret <- x
  } else if (is.numeric(x)) {
    ret <- !(x %in% 0)
  } else if (is.character(x)) {
    ret <- rep(NA, length(x))
    ret[toupper(x) %in% c("TRUE", "T", "YES", "Y")] <- TRUE
    ret[toupper(x) %in% c("FALSE", "F", "NO", "N")] <- FALSE
  } else {
    stop("Cannot handle class:", class(x))
  }
  ret[is.na(x)] <- na.value
  ret
}

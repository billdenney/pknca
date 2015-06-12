#' Check the formatting of a calculation interval specification data
#' frame.
#'
#' Calculation interval specifications are data frames with four
#' columns: \code{start}, \code{end}, \code{end.text},
#' \code{half.life}.
#'
#' The columns are defined as:
#' \describe{
#'   \item{\code{start}}{The starting time as a number}
#'   \item{\code{end}}{The ending time as either NA or a number (including Inf)}
#'   \item{\code{auc.type}}{What type of AUC should be computed: 'AUCinf',
#'     'AUClast', 'AUCall', or NA (for half-life only)}
#'   \item{\code{half.life}}{Compute the half-life for the interval (logical variable)}
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

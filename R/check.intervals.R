#' Check the formatting of a calculation interval specification data 
#' frame.
#' 
#' Calculation interval specifications are data frames defining what 
#' calculations will be required and summarized from all time intervals.
#' Note: parameters which are not requested may be calculated if it is 
#' required for (or computed at the same time as) a requested parameter.
#' 
#' \code{start} and \code{end} time must always be given as columns, and
#' the \code{start} must be before the \code{end}.  Other columns define
#' the parameters to be calculated and the groupings to apply the
#' intervals to.
#' 
#' @param x The data frame specifying what to calculate during each time
#'   interval
#' @return x The potentially updated data frame with the interval 
#'   calculation specification.
#'   
#' @seealso \code{\link{check.interval.deps}}, 
#'   \code{\link{get.parameter.deps}}, \code{\link{get.interval.cols}}
#' @export
check.interval.specification <- function(x) {
  if (!is.data.frame(x)) {
    ## Just a warning and let as.data.frame make it an error if
    ## it can't be coerced.
    warning("Interval specification must be a data.frame")
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
        # It would probably take malicious code to get here (altering
        # the intervals without using add.interval.col
        stop("Cannot assign default value for interval column", n) # nocov
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
    stop("Interval specification may not have NA for the starting time")
  if (any(x$end %in% NA))
    stop("Interval specification may not have NA for the end time")
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

#' Get all columns that depend on a parameter
#' 
#' @param x The parameter name (as a character string)
#' @return A character vector of parameter names that depend on the 
#'   parameter \code{x}.  If none depend on \code{x}, then the result 
#'   will be an empty vector.
#' @export
get.parameter.deps <- function(x) {
  all.intervals <- get.interval.cols()
  if (!(x %in% names(all.intervals))) {
    stop("x must be the name of an NCA parameter listed by the function 'get.interval.cols'")
  }
  funmap <-
    lapply(all.intervals,
           function(x) {
             if (is.na(x$FUN) &
                 is.null(x$depends)) {
               # For columnns like "start" and "end"
               retfun <- NA
             } else if (is.na(x$FUN)) {
               if (length(x$depends) == 1) {
                 # When the value is calculated by the same function as
                 # another parameter.
                 retfun <- all.intervals[[x$depends]]$FUN
               } else {
                 # It would probably take malicious code to get here (an
                 # example of malicious code could be altering the
                 # intervals without using add.interval.col)
                 stop("Invalid interval definition with no function and multiple dependencies.") # nocov
               }
             } else {
               retfun <- x$FUN
             }
             # Define a function call by its function name and the
             # changes to the formal arguments made.
             append(list(retfun), x$formalsmap)
           })
  # Find all parameters that are defined by the same function with the same 
  samefun <- function(n, funmap) {
    mask.funmap <- rep(FALSE, length(funmap))
    for (current in n) {
      for (i in seq_len(length(funmap))) {
        mask.funmap[i] <- mask.funmap[i] |
          !any(is.na(funmap[[current]][[1]]), is.na(funmap[[i]][[1]])) &
          identical(funmap[[current]], funmap[[i]])
      }
    }
    names(funmap)[mask.funmap]
  }
  searchdeps <- function(current, funmap) {
    # Find any parameters using the same function
    start <- samefun(current, funmap)
    # Find any parameters that depend on the current parameter
    ret <-
      sapply(all.intervals,
             function(x) {
               any(x$depends %in% start)
             })
    # Extract their names
    added <- setdiff(names(ret)[ret], start)
    if (length(added) > 0) {
      # Find any parameters that depend on any of those parameters
      unique(c(start, added,
               searchdeps(added, funmap)))
    } else {
      c(start, added)
    }
  }
  sort(searchdeps(x, funmap))
}

#' Take in a single row of an interval specification and return that
#' row updated with any additional calculations that must be done to
#' fulfill all dependencies.
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

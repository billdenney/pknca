#' Generate a PKNCAresults object
#'
#' This function should not be run directly.  The object is created
#' for summarization.
#'
#' @param result a data frame with NCA calculation results and groups.
#' Each row is one interval and each column is a group name or the
#' name of an NCA parameter.
#' @param data The PKNCAdata used to generate the result
#' @param exclude (optional) The name of a column with concentrations to
#'   exclude from calculations and summarization.  If given, the column
#'   should have values of \code{NA} or \code{""} for concentrations to
#'   include and non-empty text for concentrations to exclude.
#' @return A PKNCAresults object with each of the above within.
#' @family PKNCA objects
#' @export
PKNCAresults <- function(result, data, exclude) {
  result <- pknca_unit_conversion(result=result, units=data$units)
  # Add all the parts into the object
  ret <- list(result=result,
              data=data)
  if (missing(exclude)) {
    ret <- setExcludeColumn(ret, dataname="result")
  } else {
    ret <- setExcludeColumn(ret, exclude=exclude, dataname="result")
  }
  class(ret) <- c("PKNCAresults", class(ret))
  addProvenance(ret)
}

#' Extract the parameter results from a PKNCAresults and return them
#' as a data frame.
#'
#' @param x The object to extract results from
#' @param ... Ignored (for compatibility with generic
#' \code{\link{as.data.frame}})
#' @param out.format Should the output be 'long' (default) or 'wide'?
#' @return A data frame of results
#' @export
as.data.frame.PKNCAresults <- function(x, ..., out.format=c('long', 'wide')) {
  ret <- x$result
  out.format <- match.arg(out.format)
  if (out.format %in% 'wide') {
    if ("PPSTRESU" %in% names(ret)) {
      # Use standardized results
      ret$PPTESTCD <- sprintf("%s (%s)", ret$PPTESTCD, ret$PPSTRESU)
      ret$PPORRES <- ret$PPSTRES
    } else if ("PPORRESU" %in% names(ret)) {
      # Use original results
      ret$PPTESTCD <- sprintf("%s (%s)", ret$PPTESTCD, ret$PPORRESU)
    }
    # Since we moved the results into PPTESTCD and PPORRES regardless of what
    # they really are in the source data, remove the extra units and unit
    # conversion columns to allow spread to work.
    ret <- ret[, setdiff(names(ret), c("PPSTRES", "PPSTRESU", "PPORRESU"))]
    ret <- tidyr::spread(ret, key="PPTESTCD", value="PPORRES")
  }
  ret
}

#' @rdname getDataName
#' @export
getDataName.PKNCAresults <- function(object) {
  "result"
}

#' @rdname is_sparse_pk
#' @export
is_sparse_pk.PKNCAresults <- function(object) {
  is_sparse_pk(object$data)
}

#' @rdname getGroups.PKNCAconc
#' @export
getGroups.PKNCAresults <- function(object,
                                   form=formula(object$data$conc), level,
                                   data=object$result, sep) {
  # Include the start time as a group; this may be dropped later
  grpnames <- c(unlist(object$data$conc$columns$groups), "start")
  if (is_sparse_pk(object)) {
    grpnames <- setdiff(grpnames, object$data$conc$columns$subject)
  }
  if (!missing(level))
    if (is.factor(level) | is.character(level)) {
      level <- as.character(level)
      if (any(!(level %in% grpnames)))
        stop("Not all levels are listed in the group names.  Missing levels are: ",
             paste(setdiff(level, grpnames), collapse=", "))
      grpnames <- level
    } else if (is.numeric(level)) {
      if (length(level) == 1) {
        grpnames <- grpnames[1:level]
      } else {
        grpnames <- grpnames[level]
      }
    }
  data[, grpnames, drop=FALSE]
}

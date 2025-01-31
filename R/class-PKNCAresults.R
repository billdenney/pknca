#' Generate a PKNCAresults object
#'
#' This function should not be run directly.  The object is created for
#' summarization.
#'
#' @param result a data frame with NCA calculation results and groups. Each row
#'   is one interval and each column is a group name or the name of an NCA
#'   parameter.
#' @param data The PKNCAdata used to generate the result
#' @param exclude (optional) The name of a column with concentrations to exclude
#'   from calculations and summarization.  If given, the column should have
#'   values of `NA` or `""` for concentrations to include and non-empty text for
#'   concentrations to exclude.
#' @returns A PKNCAresults object with each of the above within.
#' @family PKNCA objects
#' @export
PKNCAresults <- function(result, data, exclude = NULL) {
  result <- pknca_unit_conversion(result = result, units = data$units)
  # Add all the parts into the object
  ret <- list(result=result,
              data=data)
  ret <- setExcludeColumn(ret, exclude = exclude, dataname = "result")
  class(ret) <- c("PKNCAresults", class(ret))
  addProvenance(ret)
}

#' Extract the parameter results from a PKNCAresults and return them as a
#' data.frame.
#'
#' @param x The object to extract results from
#' @param ... Ignored (for compatibility with generic [as.data.frame()])
#' @param out_format Should the output be 'long' (default) or 'wide'?
#' @param filter_requested Only return rows with parameters that were
#'   specifically requested?
#' @param filter_excluded Should excluded values be removed?
#' @param out.format Deprecated in favor of `out_format`
#' @returns A data.frame (or usually a tibble) of results
#' @export
as.data.frame.PKNCAresults <- function(x, ..., out_format = c('long', 'wide'), filter_requested = FALSE, filter_excluded = FALSE, out.format = deprecated()) {
  if (!filter_excluded) {
    ret <- x$result
  } else {
    ret <- summarize_PKNCAresults_clean_exclude(x)
    ret <- ret[is.na(ret[[x$columns$exclude]]), ]
  }
  # nocov start
  if (lifecycle::is_present(out.format)) {
    lifecycle::deprecate_warn(
      when = "0.11.0",
      what = "PKNCA::as.data.frame.PKNCAresults(out.format = )",
      with = "PKNCA::as.data.frame.PKNCAresults(out_format = )"
    )
    out_format <- out.format
  }
  # nocov end
  out_format <- match.arg(out_format)

  if (filter_requested) {
    intervals_long <-
      tidyr::pivot_longer(
        x$data$intervals,
        cols = setdiff(names(get.interval.cols()), c("start", "end")),
        names_to = "PPTESTCD",
        values_to = "keep_interval"
      )
    intervals_long_filtered <- intervals_long[intervals_long$keep_interval, , drop = FALSE]
    intervals_long_filtered$keep_interval <- NULL
    ret <-
      dplyr::inner_join(
        ret, intervals_long_filtered,
        by = intersect(names(ret), names(intervals_long_filtered))
      )
  }

  if (out_format %in% 'wide') {
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

#' @describeIn group_vars.PKNCAconc Get group_vars for a PKNCAresults object
#'   from the PKNCAconc object within
#' @exportS3Method dplyr::group_vars
group_vars.PKNCAresults <- function(x) {
  group_vars.PKNCAconc(as_PKNCAconc(x))
}

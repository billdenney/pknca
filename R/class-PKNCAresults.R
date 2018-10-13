#' Generate a PKNCAresults object
#'
#' This function should not be run directly.  The object is created
#' for summarization and plotting.
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
  ## Add all the parts into the object
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
    ret <- tidyr::spread_(ret, "PPTESTCD", "PPORRES")
  }
  ret
}

#' Extract all the original data from a PKNCAconc or PKNCAdose object
#' @param object R object to extract the data from.
#' @export
getData.PKNCAresults <- function(object)
  object$result

#' @rdname getDataName
getDataName.PKNCAresults <- function(object)
  "result"

#' @rdname getGroups.PKNCAconc
#' @export
getGroups.PKNCAresults <- function(object,
                                   form=formula(object$data$conc), level,
                                   data=object$result, sep) {
  ## Include the start time as a group; this may be dropped later
  grpnames <- c(all.vars(parseFormula(form)$groups), "start")
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

#' During the summarization of PKNCAresults, do the rounding of values
#' based on the instructions given.
#'
#' @param x The values to summarize
#' @param name The NCA parameter name (matching a parameter name in
#' \code{\link{PKNCA.set.summary}})
#' @return A string of the rounded value
#' @export
roundingSummarize <- function(x, name) {
  summaryInstructions <- PKNCA.set.summary()
  if (!(name %in% names(summaryInstructions)))
    stop(name, "is not in the summarization instructions from PKNCA.set.summary")
  roundingInstructions <- summaryInstructions[[name]]$rounding
  if (is.function(roundingInstructions)) {
    ret <- roundingInstructions(x)
  } else if (is.list(roundingInstructions)) {
    if (length(roundingInstructions) != 1)
      stop("Cannot interpret rounding instructions for ", name)
    if ("signif" == names(roundingInstructions)) {
      ret <- signifString(x, roundingInstructions$signif)
    } else if ("round" == names(roundingInstructions)) {
      ret <- roundString(x, roundingInstructions$round)
    } else {
      stop("Invalid rounding instruction list name for ", name)
    }
  }
  if (!is.character(ret))
    ret <- as.character(ret)
  ret
}

#' Summarize PKNCA results
#' 
#' @details Excluded results will not be included in the summary.
#' 
#' @param object The results to summarize
#' @param drop.group Which group(s) should be dropped from the formula?
#' @param not.requested.string A character string to use when a parameter 
#'   summary was not requested for a parameter within an interval.
#' @param not.calculated.string A character string to use when a parameter 
#'   summary was requested, but the point estimate AND spread calculations (if 
#'   applicable) returned \code{NA}.
#' @param summarize.n.per.group Should a column for \code{N} be added 
#'   (\code{TRUE} or \code{FALSE})?  Note that \code{N} is maximum number of 
#'   parameter results for any parameter; if no parameters are requested for a
#'   group, then \code{N} will be \code{NA}.
#' @param ... Ignored.
#' @return A data frame of NCA parameter results summarized according to the 
#'   summarization settings.
#' @seealso \code{\link{PKNCA.set.summary}}, \code{\link{print.summary_PKNCAresults}}
#' @export
summary.PKNCAresults <- function(object, ...,
                                 drop.group=object$data$conc$subject,
                                 summarize.n.per.group=TRUE,
                                 not.requested.string=".",
                                 not.calculated.string="NC") {
  allGroups <- getGroups(object)
  if (any(c("start", "end") %in% drop.group)) {
    warning("drop.group including start or end may result in incorrect groupings (such as inaccurate comparison of intervals).  Drop these with care.")
  }
  groups <- unique(setdiff(c("start", "end", names(allGroups)), drop.group))
  exclude_col <- object$exclude
  # Ensure that the exclude_col is NA instead of "" for subsequent processing.
  object$result[[exclude_col]] <- normalize_exclude(object$result[[exclude_col]])
  summaryFormula <- stats::as.formula(paste0("~", paste(groups, collapse="+")))
  summaryInstructions <- PKNCA.set.summary()
  ## Find any parameters that request any summaries
  resultDataCols <- 
    lapply(object$data$intervals[,
                                 setdiff(
                                   intersect(names(object$data$intervals),
                                             names(get.interval.cols())),
                                   c("start", "end")),
                                 drop=FALSE],
           FUN=any)
  resultDataCols <- as.data.frame(resultDataCols[unlist(resultDataCols)])
  ret <- unique(object$result[, groups, drop=FALSE])
  if (summarize.n.per.group) {
    ret$N <- NA
  }
  ret <- cbind(ret,
               resultDataCols)
  ret[,names(resultDataCols)] <- not.requested.string
  # Loop over every group that needs summarization
  for (i in seq_len(nrow(ret)))
    ## Loop over every column that needs summarziation
    for (n in names(resultDataCols)) {
      ## Select the rows of the intervals that match the current row
      ## from the return value.
      current.interval <-
        merge(ret[i, groups, drop=FALSE],
              object$data$intervals[,intersect(names(object$data$intervals),
                                               c(groups, n))])
      if (any(current.interval[,n])) {
        currentData <- merge(
          ret[i, groups, drop=FALSE],
          object$result[object$result$PPTESTCD %in% n,,drop=FALSE])
        currentData$PPORRES[!is.na(currentData[[exclude_col]])] <- NA
        if (nrow(currentData) == 0) {
          warning("No results to summarize for ", n, " in result row ", i)
        } else {
          if (summarize.n.per.group) {
            ret$N[i] <- max(ret$N[i], nrow(currentData), na.rm=TRUE)
          }
          ## Calculation is required
          if (is.null(summaryInstructions[[n]])) {
            stop("No summary function is set for parameter ", n, ".  Please set it with PKNCA.set.summary and report this as a bug in PKNCA.") # nocov
          }
          point <- summaryInstructions[[n]]$point(
            currentData$PPORRES)
          na.point <- is.na(point)
          na.spread <- NA
          ## Round the point estimate
          point <- roundingSummarize(point, n)
          current <- point
          if ("spread" %in% names(summaryInstructions[[n]])) {
            spread <- summaryInstructions[[n]]$spread(
              currentData$PPORRES)
            na.spread <- all(is.na(spread))
            if (na.spread) {
              ## The spread couldn't be calculated, so show that
              spread <- not.calculated.string
            } else {
              ## Round the spread
              spread <- roundingSummarize(spread, n)
            }
            ## Collapse the spread into a usable form if it is
            ## longer than one (e.g. a range or a confidence
            ## interval) and put brackets around it.
            spread <- paste0(" [", paste(spread, collapse=", "), "]")
            current <- paste0(current, spread)
          }
          ## Determine if the results were all missing, and if so, give
          ## the not.calculated.string
          if (na.point & (na.spread %in% c(NA, TRUE))) {
            ret[i,n] <- not.calculated.string
          } else {
            ret[i,n] <- current
          }
        }
      }
    }
  ## If N is requested, but it is not provided, then it should be set to not
  ## calculated.
  if (summarize.n.per.group) {
    if (any(mask.na.N <- is.na(ret$N))) {
      ret$N[mask.na.N] <- not.calculated.string
    }
    ret$N <- as.character(ret$N)
  }
  # Extract the summarization descriptions for the caption
  summary_descriptions <-
    unlist(
      lapply(
        X=summaryInstructions[names(resultDataCols)],
        FUN=`[[`,
        i="description"
      )
    )
  simplified_summary_descriptions <- summary_descriptions[!duplicated(summary_descriptions)]
  for (idx in seq_along(simplified_summary_descriptions)) {
    names(simplified_summary_descriptions)[idx] <-
      paste(names(summary_descriptions)[summary_descriptions %in% simplified_summary_descriptions[idx]],
            collapse=", ")
  }
  as_summary_PKNCAresults(
    ret,
    caption=paste(
      names(simplified_summary_descriptions),
      simplified_summary_descriptions,
      sep=": ",
      collapse="; "
    )
  )
}

as_summary_PKNCAresults <- function(data, caption) {
  structure(
    data,
    caption=caption,
    class=c("summary_PKNCAresults", "data.frame")
  )
}

#' Print the results summary
#' @param x A summary_PKNCAresults object
#' @param ... passed to print.data.frame (\code{row.names} is always set to
#'   \code{FALSE})
#' @return \code{x} invisibly
#' @seealso \code{\link{summary.PKNCAresults}}
#' @export
print.summary_PKNCAresults <- function(x, ...) {
  print.data.frame(x, row.names=FALSE, ...)
  cat("\nCaption: ", attr(x, "caption"))
  invisible(x)
}

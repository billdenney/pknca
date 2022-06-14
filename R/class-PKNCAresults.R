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
#' @importFrom tidyr spread
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
getDataName.PKNCAresults <- function(object) {
  "result"
}

#' @rdname is_sparse_pk
is_sparse_pk.PKNCAresults <- function(object) {
  is_sparse_pk(object$data)
}

#' @rdname getGroups.PKNCAconc
#' @export
getGroups.PKNCAresults <- function(object,
                                   form=formula(object$data$conc), level,
                                   data=object$result, sep) {
  # Include the start time as a group; this may be dropped later
  grpnames <- c(all.vars(parseFormula(form)$groups), "start")
  if (is_sparse_pk(object)) {
    grpnames <- setdiff(grpnames, object$data$conc$subject)
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

#' During the summarization of PKNCAresults, do the rounding of values
#' based on the instructions given.
#'
#' @param x The values to summarize
#' @param name The NCA parameter name (matching a parameter name in
#' \code{\link{PKNCA.set.summary}})
#' @return A string of the rounded value
#' @export
roundingSummarize <- function(x, name) {
  summary_instructions <- PKNCA.set.summary()
  if (!(name %in% names(summary_instructions)))
    stop(name, " is not in the summarization instructions from PKNCA.set.summary")
  roundingInstructions <- summary_instructions[[name]]$rounding
  if (is.function(roundingInstructions)) {
    ret <- roundingInstructions(x)
  } else if (is.list(roundingInstructions)) {
    if (length(roundingInstructions) != 1)
      stop("Cannot interpret rounding instructions for ", name, " (please report this as a bug)") # nocov
    if ("signif" == names(roundingInstructions)) {
      ret <- signifString(x, roundingInstructions$signif)
    } else if ("round" == names(roundingInstructions)) {
      ret <- roundString(x, roundingInstructions$round)
    } else {
      stop("Invalid rounding instruction list name for ", name, " (please report this as a bug)") # nocov
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
#' @param pretty_names Should pretty names (easier to understand in a report) be
#'   used?  \code{TRUE} is yes, \code{FALSE} is no, and \code{NULL} is yes if
#'   units are used an no if units are not used.
#' @param ... Ignored.
#' @return A data frame of NCA parameter results summarized according to the 
#'   summarization settings.
#' @seealso \code{\link{PKNCA.set.summary}}, \code{\link{print.summary_PKNCAresults}}
#' @examples
#' conc_obj <- PKNCAconc(as.data.frame(datasets::Theoph), conc~Time|Subject)
#' d_dose <- unique(datasets::Theoph[datasets::Theoph$Time == 0,
#'                                   c("Dose", "Time", "Subject")])
#' dose_obj <- PKNCAdose(d_dose, Dose~Time|Subject)
#' data_obj_automatic <- PKNCAdata(conc_obj, dose_obj)
#' results_obj_automatic <- pk.nca(data_obj_automatic)
#' # To get standard results run summary
#' summary(results_obj_automatic)
#' # To enable numeric conversion and extraction, do not give a spread function
#' # and subsequently run as.numeric on the result columns.
#' PKNCA.set.summary(
#'   name=c("auclast", "cmax", "half.life", "aucinf.obs"),
#'   point=business.geomean,
#'   description="geometric mean"
#' )
#' PKNCA.set.summary(
#'   name=c("tmax"),
#'   point=business.median,
#'   description="median"
#' )
#' summary(results_obj_automatic, not.requested.string="NA")
#' @export
summary.PKNCAresults <- function(object, ...,
                                 drop.group=object$data$conc$subject,
                                 summarize.n.per.group=TRUE,
                                 not.requested.string=".",
                                 not.calculated.string="NC",
                                 pretty_names=NULL) {
  all_group_cols <- getGroups(object)
  if (any(c("start", "end") %in% drop.group)) {
    warning("drop.group including start or end may result in incorrect groupings (such as inaccurate comparison of intervals).  Drop these with care.")
  }
  group_cols <- unique(setdiff(c("start", "end", names(all_group_cols)), drop.group))
  exclude_col <- object$exclude
  # Ensure that the exclude_col is NA instead of "" for subsequent processing.
  raw_results <- object$result
  raw_results[[exclude_col]] <- normalize_exclude(raw_results[[exclude_col]])
  summary_instructions <- PKNCA.set.summary()
  # Find any parameters that request any summaries
  parameter_cols <-
    setdiff(
      intersect(
        names(object$data$intervals),
        names(get.interval.cols())),
      c("start", "end")
    )
  # Columns that will have reported results
  result_data_cols_list <-
    lapply(
      X=object$data$intervals[, parameter_cols, drop=FALSE],
      FUN=any
    )
  result_data_cols_list <- result_data_cols_list[unlist(result_data_cols_list)]

  # Prepare for unit management
  use_units <- "PPORRESU" %in% names(object$result)
  unit_list <- NULL
  result_number_col <- intersect(c("PPSTRES", "PPORRES"), names(object$result))[1]
  if (use_units) {
    # Choose the preferred units column
    unit_col <- intersect(c("PPSTRESU", "PPORRESU"), names(object$result))[1]
    unit_list <- result_data_cols_list
    for (nm in names(unit_list)) {
      # Get all of the units for a given parameter that will be summarized
      unit_list[[nm]] <- unique(object$result[[unit_col]][object$result$PPTESTCD %in% nm])
    }
  }
  if (is.null(pretty_names)) {
    pretty_names <- !is.null(unit_list)
  }
  
  result_data_cols <- as.data.frame(result_data_cols_list)
  # If no other value is filled in, then the default is that it was not
  # requested.
  result_data_cols[, names(result_data_cols)] <- not.requested.string
  # Rows that will have results
  ret_group_cols <- unique(raw_results[, group_cols, drop=FALSE])
  simplified_results <-
    raw_results[raw_results$PPTESTCD %in% names(result_data_cols), , drop=FALSE]
  ret <- unique(raw_results[, group_cols, drop=FALSE])
  if (summarize.n.per.group) {
    ret$N <- NA_integer_
  }
  ret <- cbind(ret, result_data_cols)
  # Loop over every group that needs summarization
  for (row_idx in seq_len(nrow(ret)))
    # Loop over every column that needs summarization
    for (current_parameter in names(result_data_cols)) {
      # Select the rows of the intervals that match the current row
      # from the return value.
      current_interval <-
        merge(
          ret[row_idx, group_cols, drop=FALSE],
          object$data$intervals[,
                                intersect(names(object$data$intervals),
                                          c(group_cols, current_parameter)),
                                drop=FALSE]
        )
      if (any(current_interval[,current_parameter])) {
        current_data <- merge(
          ret[row_idx, group_cols, drop=FALSE],
          simplified_results[simplified_results$PPTESTCD %in% current_parameter,,drop=FALSE])
        # Exclude value, when required
        current_data[[result_number_col]][!is.na(current_data[[exclude_col]])] <- NA
        if (nrow(current_data) == 0) {
          # I don't think that a user can get here
          warning("No results to summarize for ", current_parameter, " in result row ", row_idx) # nocov
        } else {
          if (summarize.n.per.group) {
            ret$N[row_idx] <- max(ret$N[row_idx], nrow(current_data), na.rm=TRUE)
          }
          # Calculation is required
          if (is.null(summary_instructions[[current_parameter]])) {
            stop("No summary function is set for parameter ", current_parameter, ".  Please set it with PKNCA.set.summary and report this as a bug in PKNCA.") # nocov
          }
          point <- summary_instructions[[current_parameter]]$point(current_data[[result_number_col]])
          na_point <- is.na(point)
          na_spread <- NA
          # Round the point estimate
          point <- roundingSummarize(point, current_parameter)
          current <- point
          if ("spread" %in% names(summary_instructions[[current_parameter]])) {
            spread <- summary_instructions[[current_parameter]]$spread(
              current_data[[result_number_col]])
            na_spread <- all(is.na(spread))
            if (na_spread) {
              # The spread couldn't be calculated, so show that
              spread <- not.calculated.string
            } else {
              # Round the spread
              spread <- roundingSummarize(spread, current_parameter)
            }
            # Collapse the spread into a usable form if it is longer than one
            # (e.g. a range or a confidence interval) and put brackets around
            # it.
            spread <- paste0(" [", paste(spread, collapse=", "), "]")
            current <- paste0(current, spread)
          }
          # Determine if the results were all missing, and if so, give
          # the not.calculated.string
          if (na_point & (na_spread %in% c(NA, TRUE))) {
            ret[row_idx, current_parameter] <- not.calculated.string
          } else {
            if (use_units) {
              if (length(unit_list[[current_parameter]]) > 1) {
                # Need to choose the correct, current unit, and if more than one
                # is present, do not summarize.
                units_to_add <- unique(current_data[[unit_col]])
                if (length(units_to_add) > 1) {
                  stop(
                    "Multiple units cannot be summarized together.  For ",
                    current_parameter, ", trying to combine: ",
                    paste(units_to_add, collapse=", ")
                  )
                }
                current <- paste(current, units_to_add)
              }
            }
            ret[row_idx, current_parameter] <- current
          }
        }
      }
    }
  # If N is requested, but it is not provided, then it should be set to not
  # calculated.
  if (summarize.n.per.group) {
    if (any(mask.na.N <- is.na(ret$N))) {
      #ret$N[mask.na.N] <- not.calculated.string
      stop("Invalid subject count (please report this as a bug)") # nocov
    }
    ret$N <- as.character(ret$N)
  }
  # Extract the summarization descriptions for the caption
  summary_descriptions <-
    unlist(
      lapply(
        X=summary_instructions[names(result_data_cols)],
        FUN=`[[`,
        i="description"
      )
    )
  if (pretty_names) {
    # Make the caption use pretty names if they're used in the header
    all_intervals <- get.interval.cols()
    for (idx in seq_along(summary_descriptions)) {
      names(summary_descriptions)[idx] <- all_intervals[[names(summary_descriptions)[idx]]]$pretty_name
    }
  }
  simplified_summary_descriptions <- summary_descriptions[!duplicated(summary_descriptions)]
  for (idx in seq_along(simplified_summary_descriptions)) {
    names(simplified_summary_descriptions)[idx] <-
      paste(names(summary_descriptions)[summary_descriptions %in% simplified_summary_descriptions[idx]],
            collapse=", ")
  }
  ret_pretty <- rename_summary_PKNCAresults(data=ret, unit_list=unit_list, pretty_names=pretty_names)
  as_summary_PKNCAresults(
    ret_pretty,
    caption=paste(
      names(simplified_summary_descriptions),
      simplified_summary_descriptions,
      sep=": ",
      collapse="; "
    )
  )
}

rename_summary_PKNCAresults <- function(data, unit_list, pretty_names) {
  units_to_use <-
    setNames(rep(NA_character_, ncol(data)), names(data))
  if (!is.null(unit_list)) {
    # add the units to the column header, if applicable
    for (nm in names(unit_list)) {
      if (length(unit_list[[nm]]) == 1) {
        units_to_use[nm] <- unit_list[[nm]]
      }
    }
  }
  pretty_names_to_use <-
    setNames(rep(NA_character_, ncol(data)), names(data))
  if (pretty_names) {
    all_intervals <- get.interval.cols()
    for (nm in names(pretty_names_to_use)) {
      if (!is.null(all_intervals[[nm]]$pretty_name)) {
        pretty_names_to_use[nm] <- all_intervals[[nm]]$pretty_name
      }
    }
  }
  for (idx in seq_len(ncol(data))) {
    current_col <- names(data)[idx]
    first_part <- na.omit(c(pretty_names_to_use[current_col], current_col))[1]
    unit_part <- units_to_use[current_col]
    names(data)[idx] <-
      if (is.na(unit_part)) {
        first_part
      } else {
        sprintf("%s (%s)", first_part, unit_part)
      }
  }
  data
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
  cat(paste0("\nCaption: ", attr(x, "caption"), "\n"), fill=TRUE)
  invisible(x)
}

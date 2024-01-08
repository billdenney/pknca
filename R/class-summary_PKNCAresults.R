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
#'   (\code{TRUE} or \code{FALSE})?  \code{NA} means to automatically detect
#'   adding \code{N} if the data has a subject column indicated.  Note that
#'   \code{N} is maximum number of parameter results for any parameter; if no
#'   parameters are requested for a group, then \code{N} will be \code{NA}.
#' @param pretty_names Should pretty names (easier to understand in a report) be
#'   used?  \code{TRUE} is yes, \code{FALSE} is no, and \code{NULL} is yes if
#'   units are used an no if units are not used.
#' @param ... Ignored.
#' @return A data frame of NCA parameter results summarized according to the
#'   summarization settings.
#' @seealso \code{\link{PKNCA.set.summary}}, \code{\link{print.summary_PKNCAresults}}
#' @examples
#' conc_obj <- PKNCAconc(as.data.frame(datasets::Theoph), conc ~ Time | Subject)
#' d_dose <-
#'   unique(datasets::Theoph[
#'     datasets::Theoph$Time == 0,
#'     c("Dose", "Time", "Subject")
#'   ])
#' dose_obj <- PKNCAdose(d_dose, Dose ~ Time | Subject)
#' data_obj_automatic <- PKNCAdata(conc_obj, dose_obj)
#' results_obj_automatic <- pk.nca(data_obj_automatic)
#' # To get standard results run summary
#' summary(results_obj_automatic)
#' # To enable numeric conversion and extraction, do not give a spread function
#' # and subsequently run as.numeric on the result columns.
#' PKNCA.set.summary(
#'   name = c("auclast", "cmax", "half.life", "aucinf.obs"),
#'   point = business.geomean,
#'   description = "geometric mean"
#' )
#' PKNCA.set.summary(
#'   name = c("tmax"),
#'   point = business.median,
#'   description = "median"
#' )
#' summary(results_obj_automatic, not.requested.string = "NA")
#' @export
summary.PKNCAresults <- function(object, ...,
                                 drop.group = object$data$conc$columns$subject,
                                 summarize.n.per.group = NA,
                                 not.requested.string = ".",
                                 not.calculated.string = "NC",
                                 pretty_names = NULL) {
  all_group_cols <- getGroups(object)
  if (any(c("start", "end") %in% drop.group)) {
    warning("drop.group including start or end may result in incorrect groupings (such as inaccurate comparison of intervals).  Drop these with care.")
  }
  group_cols <- unique(setdiff(c("start", "end", names(all_group_cols)), drop.group))
  exclude_col <- object$columns$exclude
  # Ensure that the exclude_col is NA instead of "" for subsequent processing.
  raw_results <- object$result
  raw_results[[exclude_col]] <- normalize_exclude(raw_results[[exclude_col]])
  summary_instructions <- PKNCA.set.summary()
  # Find any parameters that request any summaries
  parameter_cols <-
    setdiff(
      intersect(
        names(object$data$intervals),
        names(get.interval.cols())
      ),
      c("start", "end")
    )
  # Columns that will have reported results
  result_data_cols_list <-
    lapply(
      X = object$data$intervals[, parameter_cols, drop = FALSE],
      FUN = any
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
  ret_group_cols <- unique(raw_results[, group_cols, drop = FALSE])
  simplified_results <-
    raw_results[raw_results$PPTESTCD %in% names(result_data_cols), , drop = FALSE]
  ret <- unique(raw_results[, group_cols, drop = FALSE])

  subject_col <- object$data$conc$columns$subject
  has_subject_col <- length(subject_col) > 0
  if (is.na(summarize.n.per.group)) {
    summarize.n.per.group <- has_subject_col
  } else if (summarize.n.per.group & !has_subject_col) {
    warning("summarize.n.per.group was requested, but no subject column exists")
    summarize.n.per.group <- FALSE
  }
  if (summarize.n.per.group) {
    ret$N <- NA_integer_
  }

  ret <- cbind(ret, result_data_cols)
  # Loop over every group that needs summarization
  for (row_idx in seq_len(nrow(ret))) {
    # Loop over every column that needs summarization
    for (current_parameter in names(result_data_cols)) {
      # Select the rows of the intervals that match the current row
      # from the return value.
      current_interval <-
        merge(
          ret[row_idx, group_cols, drop = FALSE],
          object$data$intervals[,
            intersect(
              names(object$data$intervals),
              c(group_cols, current_parameter)
            ),
            drop = FALSE
          ]
        )
      if (any(current_interval[, current_parameter])) {
        current_data <- merge(
          ret[row_idx, group_cols, drop = FALSE],
          simplified_results[simplified_results$PPTESTCD %in% current_parameter, , drop = FALSE]
        )
        # Exclude value, when required
        current_data[[result_number_col]][!is.na(current_data[[exclude_col]])] <- NA
        if (nrow(current_data) == 0) {
          # I don't think that a user can get here
          warning("No results to summarize for ", current_parameter, " in result row ", row_idx) # nocov
        } else {
          if (summarize.n.per.group) {
            n_subjects <- length(unique(current_data[[subject_col]]))
            if (n_subjects < nrow(current_data)) {
              warning("Some subjects may have more than one result for ", current_parameter)
            }
            # max() because the summary table provides the N for the full row, not
            # for a single parameter.
            ret$N[row_idx] <- max(ret$N[row_idx], n_subjects, na.rm = TRUE)
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
              current_data[[result_number_col]]
            )
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
            spread <- paste0(" [", paste(spread, collapse = ", "), "]")
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
                    paste(units_to_add, collapse = ", ")
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
  }
  # If N is requested, but it is not provided, then it should be set to not
  # calculated.
  if (summarize.n.per.group) {
    if (any(mask.na.N <- is.na(ret$N))) {
      # ret$N[mask.na.N] <- not.calculated.string
      stop("Invalid subject count (please report this as a bug)") # nocov
    }
    ret$N <- as.character(ret$N)
  }
  # Extract the summarization descriptions for the caption
  summary_descriptions <-
    unlist(
      lapply(
        X = summary_instructions[names(result_data_cols)],
        FUN = `[[`,
        i = "description"
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
        collapse = ", "
      )
  }
  ret_pretty <- rename_summary_PKNCAresults(data = ret, unit_list = unit_list, pretty_names = pretty_names)
  as_summary_PKNCAresults(
    ret_pretty,
    caption = paste(
      names(simplified_summary_descriptions),
      simplified_summary_descriptions,
      sep = ": ",
      collapse = "; "
    )
  )
}

rename_summary_PKNCAresults <- function(data, unit_list, pretty_names) {
  units_to_use <-
    stats::setNames(rep(NA_character_, ncol(data)), names(data))
  if (!is.null(unit_list)) {
    # add the units to the column header, if applicable
    for (nm in names(unit_list)) {
      if (length(unit_list[[nm]]) == 1) {
        units_to_use[nm] <- unit_list[[nm]]
      }
    }
  }
  pretty_names_to_use <-
    stats::setNames(rep(NA_character_, ncol(data)), names(data))
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
    first_part <- stats::na.omit(c(pretty_names_to_use[current_col], current_col))[1]
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
    caption = caption,
    class = c("summary_PKNCAresults", "data.frame")
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
  print.data.frame(x, row.names = FALSE, ...)
  cat(paste0("\nCaption: ", attr(x, "caption"), "\n"), fill = TRUE)
  invisible(x)
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
  if (!(name %in% names(summary_instructions))) {
    stop(name, " is not in the summarization instructions from PKNCA.set.summary")
  }
  roundingInstructions <- summary_instructions[[name]]$rounding
  if (is.function(roundingInstructions)) {
    ret <- roundingInstructions(x)
  } else if (is.list(roundingInstructions)) {
    if (length(roundingInstructions) != 1) {
      stop("Cannot interpret rounding instructions for ", name, " (please report this as a bug)") # nocov
    }
    if ("signif" == names(roundingInstructions)) {
      ret <- signifString(x, roundingInstructions$signif)
    } else if ("round" == names(roundingInstructions)) {
      ret <- roundString(x, roundingInstructions$round)
    } else {
      stop("Invalid rounding instruction list name for ", name, " (please report this as a bug)") # nocov
    }
  }
  if (!is.character(ret)) {
    ret <- as.character(ret)
  }
  ret
}

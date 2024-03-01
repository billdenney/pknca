#' Summarize PKNCA results
#'
#' @details Excluded results will not be included in the summary.
#'
#' @param object The results to summarize
#' @param drop_group Which group(s) should be dropped from the formula?
#' @param drop_param Which parameters should be excluded from the summary?
#' @param not_requested A character string to use when a parameter summary was
#'   not requested for a parameter within an interval.
#' @param not_calculated A character string to use when a parameter summary was
#'   requested, but the point estimate AND spread calculations (if applicable)
#'   returned `NA`.
#' @param summarize_n Should a column for `N` be added (`TRUE` or `FALSE`)?
#'   `NA` means to automatically detect adding `N` if the data has a subject
#'   column indicated.  Note that `N` is maximum number of parameter results for
#'   any parameter; if no parameters are requested for a group, then `N` will be
#'   `NA`.
#' @param pretty_names Should pretty names (easier to understand in a report) be
#'   used?  `TRUE` is yes, `FALSE` is no, and `NULL` is yes if units are used
#'   and no if units are not used.
#' @param ... Ignored.
#' @param drop.group,summarize.n.per.group,not.requested.string,not.calculated.string
#' Deprecated use `drop_group`, `not_requested`, `not_calculated`, or
#' `summarize_n`, instead
#' @returns A data frame of NCA parameter results summarized according to the
#'   summarization settings.
#' @seealso [PKNCA.set.summary()], [print.summary_PKNCAresults()]
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
#' summary(results_obj_automatic, not_requested = "NA")
#' @export
#' @importFrom lifecycle deprecated
summary.PKNCAresults <- function(object, ...,
                                 drop_group = object$data$conc$columns$subject,
                                 drop_param = character(),
                                 summarize_n = NA,
                                 not_requested = ".",
                                 not_calculated = "NC",

                                 drop.group = deprecated(),
                                 summarize.n.per.group = deprecated(),
                                 not.requested.string = deprecated(),
                                 not.calculated.string = deprecated(),
                                 pretty_names = NULL) {
  # Process inputs ####

  ## Deprecated inputs ####
  # nocov start
  if (lifecycle::is_present(drop.group)) {
    lifecycle::deprecate_warn(
      when = "0.11.0",
      what = "PKNCA::summary.PKNCAresults(drop.group = )",
      with = "PKNCA::summary.PKNCAresults(drop_group = )"
    )
    drop_group <- drop.group
  }
  if (lifecycle::is_present(summarize.n.per.group)) {
    lifecycle::deprecate_warn(
      when = "0.11.0",
      what = "PKNCA::summary.PKNCAresults(summarize.n.per.group = )",
      with = "PKNCA::summary.PKNCAresults(summarize_n = )"
    )
    summarize_n <- summarize.n.per.group
  }
  if (lifecycle::is_present(not.requested.string)) {
    lifecycle::deprecate_warn(
      when = "0.11.0",
      what = "PKNCA::summary.PKNCAresults(not.requested.string = )",
      with = "PKNCA::summary.PKNCAresults(not_requested = )"
    )
    not_requested <- not.requested.string
  }
  if (lifecycle::is_present(not.calculated.string)) {
    lifecycle::deprecate_warn(
      when = "0.11.0",
      what = "PKNCA::summary.PKNCAresults(not.calculated.string = )",
      with = "PKNCA::summary.PKNCAresults(not_calculated = )"
    )
    not_calculated <- not.calculated.string
  }
  # nocov end

  ## Simple inputs ####
  group_cols <- get_summary_PKNCAresults_drop_group(object = object, drop_group = drop_group)
  subject_col <- object$data$conc$columns$subject
  has_subject_col <- length(subject_col) > 0
  if (is.na(summarize_n)) {
    summarize_n <- has_subject_col
  } else if (summarize_n & !has_subject_col) {
    warning("summarize_n was requested, but no subject column exists")
    summarize_n <- FALSE
  }

  # Preparation ####

  # Set excluded rows to NA, give the cleaned data.frame
  raw_results <- summarize_PKNCAresults_clean_exclude(object)

  # Find any parameters that request any summaries, and exclude ones that are
  # not requested
  parameter_cols <-
    setdiff(
      intersect(
        names(object$data$intervals),
        names(get.interval.cols())
      ),
      c(c("start", "end"), drop_param)
    )

  # Extract columns that have been requested by the user for summary in any
  # intervals
  result_data_cols_list <-
    lapply(
      X = object$data$intervals[, parameter_cols, drop = FALSE],
      FUN = any
    )
  # Then, filter them the the ones that have any "TRUE" values
  result_data_cols_list <- result_data_cols_list[unlist(result_data_cols_list)]

  # Prepare for unit management
  unit_col <- get_summary_PKNCAresults_result_unit_col(object = object)
  unit_list <- get_summary_PKNCAresults_result_unit_list(data = raw_results, unit_col = unit_col)
  if (is.null(pretty_names)) {
    pretty_names <- !is.null(unit_list)
  }

  ## Prepare the result data.frame ####
  # Populate the groups on the left and the results on the right.
  result_values <- as.data.frame(result_data_cols_list)
  # If no other value is filled in, then the default is that it was not
  # requested.
  result_values[, names(result_values)] <- not_requested
  # Rows that will have results
  result_groups <- unique(raw_results[, group_cols, drop = FALSE])
  # result_n will be the same as result_groups if summarize_n is FALSE and have
  # an "N" column added if summarize_n is TRUE.
  result_n <-
    get_summary_PKNCAresults_count_N(
      data = raw_results,
      result_group = result_groups,
      subject_col = subject_col,
      summarize_n = summarize_n,
      not_calculated = not_calculated
    )

  # Calculation ####
  ret <-
    summarize_PKNCAresults_object(
      data = raw_results,
      result_group = result_n,
      subject_col = subject_col,
      result_value_template = result_values,
      result_units = unit_list,
      intervals = object$data$intervals,
      not_calculated = not_calculated
    )

  # Decoration ####
  caption <-
    get_summary_PKNCAresults_caption(
      param_names = names(result_values),
      pretty_names = pretty_names,
      footnote_N = "N" %in% names(ret),
      footnote_n = attr(ret, "footnote_n", exact = TRUE)
    )
  attr(ret, "footnote_n") <- NULL
  ret_pretty <- rename_summary_PKNCAresults(data = ret, unit_list = unit_list, pretty_names = pretty_names)
  as_summary_PKNCAresults(
    ret_pretty,
    caption = caption
  )
}

# A helper function for summary.PKNCAresults to find the groups to drop
get_summary_PKNCAresults_drop_group <- function(object, drop_group) {
  all_group_cols <- getGroups(object)
  if (any(c("start", "end") %in% drop_group)) {
    warning("drop.group including start or end may result in incorrect groupings (such as inaccurate comparison of intervals).  Drop these with care.")
  }
  ret <-
    unique(
      setdiff(
        c("start", "end", object$data$options$keep_interval_cols, names(all_group_cols)),
        drop_group
      )
    )
  ret
}

# Get the column name with the results to use for summarization
get_summary_PKNCAresults_result_number_col <- function(object) {
  if (is.data.frame(object)) {
    data <- object
  } else {
    data <- object$result
  }
  intersect(c("PPSTRES", "PPORRES"), names(data))[1]
}

# Get the column name with the result unitss to use for summarization
get_summary_PKNCAresults_result_unit_col <- function(object) {
  if (is.data.frame(object)) {
    data <- object
  } else {
    data <- object$result
  }
  # This will return NULL if neither column is present
  ret <- intersect(c("PPSTRESU", "PPORRESU"), names(data))[1]
  if (is.na(ret)) {
    ret <- NULL
  }
  ret
}

# Get the list of units for each parameter to use for summarization
get_summary_PKNCAresults_result_unit_list <- function(data, unit_col) {
  unit_list <- NULL
  if (!is.null(unit_col)) {
    all_params <- unique(data$PPTESTCD)
    unit_list <-
      stats::setNames(
        as.list(rep(NA_character_, length(all_params))),
        nm = all_params
      )
    for (nm in names(unit_list)) {
      # Get all of the units for a given parameter that will be summarized
      unit_list[[nm]] <- unique(data[[unit_col]][data$PPTESTCD %in% nm])
    }
  }
  unit_list
}

# Count the number of subjects in each group, always return a data.frame with
# the right number of rows so that it can go into cbind
get_summary_PKNCAresults_count_N <- function(data, result_group, subject_col, summarize_n, not_calculated) {
  # subject_col %in% names(data) handles sparse data where the subject column
  # may be imputed rather than present in the data.  And, adding the `any()`
  # outside that ensures that single-subject data (with no subject column) also
  # returns FALSE.
  if (summarize_n && any(subject_col %in% names(data))) {
    # R CMD Check hack
    N <- NULL
    ret <-
      data |>
      dplyr::grouped_df(vars = names(result_group)) |>
      dplyr::summarize(
        N = length(unique(.data[[subject_col]]))
      ) |>
      dplyr::ungroup()
    # Reorder the return value to be in the same order as the original groups
    key_col <- paste0(max(names(ret)), "X")
    # This could fail if there are no groups.  But also, `summarize_n` should be
    # FALSE in that case.  So, it shouldn't be a problem.
    group_key <- do.call(paste, result_group)
    ret[[key_col]] <-
      factor(
        do.call(paste, ret[, names(result_group)]),
        levels = group_key,
        ordered = TRUE
      )
    ret <- ret[order(ret[[key_col]]), ]
    ret[[key_col]] <- NULL

    ret$N <- as.character(ret$N)
    if (any(is.na(ret$N))) {
      # If N is requested, but it is not provided, then it should be set to not
      # calculated.
      ret$N[is.na(ret$N)] <- not_calculated
    }
  } else {
    ret <- result_group
  }
  ret
}

# Provide a clean caption for summarized parameters
get_summary_PKNCAresults_caption <- function(param_names, pretty_names, footnote_N, footnote_n) {
  # Extract the summarization descriptions for the caption
  summary_descriptions <-
    unlist(
      lapply(
        X = PKNCA.set.summary()[param_names],
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
  # Remove duplicate descriptions and then find all the ones that are duplicates
  simplified_summary_descriptions <- summary_descriptions[!duplicated(summary_descriptions)]
  for (idx in seq_along(simplified_summary_descriptions)) {
    names(simplified_summary_descriptions)[idx] <-
      paste(names(summary_descriptions)[summary_descriptions %in% simplified_summary_descriptions[idx]],
            collapse = ", "
      )
  }
  ret <-
    paste(
      names(simplified_summary_descriptions),
      simplified_summary_descriptions,
      sep = ": ",
      collapse = "; "
    )
  if (footnote_N) {
    ret <- c(ret, "N: number of subjects")
  }
  if (footnote_n) {
    ret <- c(ret, "n: number of measurements included in summary")
  }
  paste(ret, collapse = "; ")
}

#' Clean up the exclusions in the object
#'
#' @param object The PKNCAresults object to clean
#' @returns The results data.frame with exclusions set to `NA`
#' @noRd
summarize_PKNCAresults_clean_exclude <- function(object) {
  result_number_col <- get_summary_PKNCAresults_result_number_col(object)
  exclude_col <- object$columns$exclude
  # Ensure that the exclude_col is NA instead of "" for subsequent processing.
  raw_results <- object$result
  raw_results[[exclude_col]] <- normalize_exclude(raw_results[[exclude_col]])
  raw_results[[result_number_col]][!is.na(raw_results[[exclude_col]])] <- NA
  raw_results
}

# A helper function for summary.PKNCAresults to summarize everything
summarize_PKNCAresults_object <- function(data, result_group, subject_col, result_value_template, result_units, intervals, not_calculated) {
  ret_values_list <- list()
  footnote_n <- FALSE
  for (idx in seq_len(nrow(result_group))) {
    ret_idx <-
      summarize_PKNCAresults_group(
        data = data,
        current_group = result_group[idx,],
        subject_col = subject_col,
        result_value_template = result_value_template,
        result_units = result_units,
        intervals = intervals,
        not_calculated = not_calculated
      )
    ret_values_list <- append(ret_values_list, list(ret_idx))
    if (attr(ret_idx, "footnote_n", exact = TRUE)) {
      footnote_n <- TRUE
    }
  }
  ret <- cbind(result_group, dplyr::bind_rows(ret_values_list))
  attr(ret, "footnote_n") <- footnote_n
  ret
}

# A helper function for summary.PKNCAresults to summarize one group
summarize_PKNCAresults_group <- function(data, current_group, subject_col, result_value_template, result_units, intervals, not_calculated) {
  ret <- result_value_template

  current_data <- dplyr::inner_join(data, current_group, by = intersect(names(data), names(current_group)))
  if (nrow(current_data) == 0) {
    # I don't think that a user can get here
    warning("No results to summarize for result row, please report a bug") # nocov
    return(ret) # nocov
  }
  current_interval <- dplyr::inner_join(intervals, current_group, by = intersect(names(intervals), names(current_group)))
  current_param_prep <-
    vapply(
      X = current_interval[, setdiff(names(get.interval.cols()), c("start", "end"))],
      FUN = any,
      FUN.VALUE = TRUE
    )
  current_param_all <-
    intersect(
      names(current_param_prep[current_param_prep]),
      # This ensures that parameters that were dropped with drop_param
      # previously are not summarized
      names(ret)
    )

  footnote_n <- FALSE
  for (current_param in current_param_all) {
    current_summary <-
      summarize_PKNCAresults_parameter(
        data = current_data,
        subject_col = subject_col,
        parameter = current_param,
        include_units = length(result_units[[current_param]]) > 1,
        not_calculated = not_calculated
      )
    # summarize N, if requested and there is a value for calculation
    if (("N" %in% names(current_group)) && (current_summary != not_calculated)) {
      N_group <- as.integer(current_group$N)
      n_summary <- attr(current_summary, "n", exact = TRUE)
      if (N_group != n_summary) {
        current_summary <- sprintf("%s, n=%d", current_summary, n_summary)
        footnote_n <- TRUE
      }
    }
    # use `as.character()` to ensure that the output stays a data.frame rather
    # than a complicated list.
    ret[[current_param]] <- as.character(current_summary)
  }
  attr(ret, "footnote_n") <- footnote_n
  ret
}

# A helper function for summary.PKNCAresults to summarize one parameter within a
# group
summarize_PKNCAresults_parameter <- function(data, parameter, subject_col, include_units, not_calculated) {
  current_data <- data[data$PPTESTCD %in% parameter, , drop = FALSE]
  number_col <- get_summary_PKNCAresults_result_number_col(data)
  unit_col <- get_summary_PKNCAresults_result_unit_col(data)

  units <- NULL
  if (!is.null(unit_col)) {
    units <- unique(current_data[[unit_col]])
    if (length(units) > 1) {
      stop(
        "Multiple units cannot be summarized together.  For ",
        parameter, ", trying to combine: ",
        paste(units, collapse = ", ")
      )
    }
  }

  if (length(subject_col) == 1) {
    N <- length(unique(current_data[[subject_col]]))
    if (any(duplicated(current_data[[subject_col]]))) {
      warning("Some subjects may have more than one result for ", parameter)
    }
  } else {
    N <- NULL
  }
  n <- sum(!is.na(current_data[[number_col]]))

  current_summary_instructions <- PKNCA.set.summary()[[parameter]]
  if (is.null(current_summary_instructions)) {
    stop("No summary function is set for parameter ", parameter, ".  Please set it with PKNCA.set.summary and report this as a bug in PKNCA.") # nocov
  }

  point <- current_summary_instructions$point(current_data[[number_col]])
  # We could count only the number of measurements included in the point
  # estimate which may differ when zeros are excluded, but that would likely be
  # less clear.  So, we are counting all points with available measurements.
  # point_n <- attr(point, which = "n", exact = TRUE)
  # if (!is.null(point_n)) {
  #   n <- point_n
  # }
  point_txt <- roundingSummarize(point, parameter)
  na_point <- is.na(point)
  result_txt <- point_txt

  spread <- NULL
  spread_txt <- NULL
  na_spread <- TRUE
  if ("spread" %in% names(current_summary_instructions) && n > 1) {
    spread <- current_summary_instructions$spread(current_data[[number_col]])
    na_spread <- all(is.na(spread))
    if (na_spread) {
      # The spread couldn't be calculated, so show that
      spread_txt <- not_calculated
    } else {
      # Round the spread
      spread_txt <- roundingSummarize(spread, parameter)
    }
    # Collapse the spread into a usable form if it is longer than one
    # (e.g. a range or a confidence interval) and put brackets around
    # it.
    spread_txt <- paste0(" [", paste(spread_txt, collapse = ", "), "]")
    result_txt <- paste0(point_txt, spread_txt)
  }

  if (na_point & na_spread) {
    result_txt <- not_calculated
  } else if (include_units) {
    result_txt <- paste(result_txt, units)
  }

  structure(
    result_txt,
    parameter = parameter,
    point = point,
    spread = spread,
    N = N,
    n = n,
    units = units,
    class = "summarize_PKNCAresults_parameter"
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
#' @param ... passed to print.data.frame (`row.names` is always set to `FALSE`)
#' @returns `x` invisibly
#' @seealso [summary.PKNCAresults()]
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
#'   [PKNCA.set.summary()])
#' @returns A string of the rounded value
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

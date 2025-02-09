#' Add specified imputation methods to the intervals in a PKNCAdata object.
#'
#' @param data A PKNCAdata object containing the intervals and data components, or a data frame of intervals.
#' @param target_impute A character string specifying the imputation method to be added.
#' @param after Numeric value specifying the position after which the imputation method should be added (optional). First is 0, last Inf. If missing, the imputation method is added at the end (Inf).
#' @param target_params A character vector specifying the parameters to be targeted (optional). If missing, all TRUE in the intervals are taken.
#' @param target_groups A data frame specifying the intervals to be targeted (optional). If missing, all relevant groups are considered.
#' @param impute_column A character string specifying the name of the impute column (optional). If missing, the default name "impute" is used.
#' @param allow_duplication A boolean specifying whether to allow creating duplicates of the target_impute in the impute column (optional). Default is TRUE.
#' @param new_rows_after_original A boolean specifying whether the new rows should be added after the original rows (optional). Default is TRUE.
#' @return A modified PKNCAdata object with the specified imputation methods added to the targeted intervals.
#' @examples
#' d_conc <- data.frame(
#'   conc = c(1, 0.6, 0.2, 0.1, 0.9, 0.4, 1.2, 0.8, 0.3, 0.2, 1.1, 0.5),
#'   time = rep(0:5, 2),
#'   ID = rep(1:2, each = 6),
#'   analyte = rep(c("Analyte1", "Analyte2"), each = 6),
#'   include_hl = c(FALSE, NA, TRUE, TRUE, TRUE, TRUE, FALSE, NA, TRUE, TRUE, TRUE, TRUE)
#' )
#'
#' d_dose <- data.frame(
#'   dose = c(100, 200),
#'   time = c(0, 0),
#'   ID = c(1, 2)
#' )
#'
#' o_conc <- PKNCAconc(d_conc, conc ~ time | ID / analyte)
#' o_dose <- PKNCAdose(d_dose, dose ~ time | ID)
#'
#' intervals <- data.frame(
#'   start = c(0, 0, 0),
#'   end = c(3, 5, Inf),
#'   half.life = c(TRUE, FALSE, TRUE),
#'   cmax = c(TRUE, TRUE, TRUE),
#'   impute = c("start_conc0,start_predose", "start_predose", "start_conc0"),
#'   analyte = c("Analyte1", "Analyte2", "Analyte1")
#' )
#'
#' o_data <- PKNCAdata(o_conc, o_dose, intervals = intervals)
#'
#' # Apply interval_add_impute function
#' o_data <- interval_add_impute(o_data, target_impute = "start_conc0", target_params = c("half.life"), target_groups = data.frame(analyte = "Analyte1"))
#' @export
interval_add_impute <- function(data, target_impute, after = Inf, target_params = NULL, target_groups = NULL, impute_column = NULL, allow_duplication = TRUE, new_rows_after_original = TRUE) {
  if (missing(data) || missing(target_impute)) {
    stop("Both 'data' and 'target_impute' must be provided.")
  }
  if (!inherits(data, "PKNCAdata") && !is.data.frame(data)) {
    stop("The 'data' object must be a PKNCAdata object or a data frame.")
  }
  if (!is.character(target_impute)) {
    stop("'target_impute' must be a character string.")
  }
  UseMethod("interval_add_impute")
}

#' @export
interval_add_impute.PKNCAdata <- function(data, target_impute, after = Inf, target_params = NULL, target_groups = NULL, impute_column = NULL, allow_duplication = TRUE, new_rows_after_original = TRUE) {
  if (is.null(impute_column) && !is.na(data$impute)) {
    impute_column <- data$impute
  }
  data$intervals <- interval_add_impute(data$intervals, target_impute, after, target_params, target_groups, impute_column, allow_duplication, new_rows_after_original)
  data
}

#' @export
interval_add_impute.data.frame <- function(data, target_impute, after = Inf, target_params = NULL, target_groups = NULL, impute_column = NULL, allow_duplication = TRUE, new_rows_after_original = TRUE) {
  
  # Add an index column to preserve the original order
  data <- dplyr::mutate(data, index = dplyr::row_number())
  
  # Get all parameter column names in the data frame
  all_param_options <- setdiff(names(PKNCA::get.interval.cols()), c("start", "end"))
  logical_cols <- names(which(colSums(data[sapply(data, is.logical)]) > 1))
  param_cols <- intersect(logical_cols, all_param_options)
  
  # Handle target_params
  if (is.null(target_params)) {
    # Take all logical columns in data that are known parameters
    target_params <- param_cols
  } else {
    # Check that all target_params are logical columns in data and known parameters
    missing_params <- setdiff(target_params, logical_cols)
    if (length(missing_params) > 0) {
      stop("The following target_params are not interval columns and/or known PKNCA parameters: ", paste(missing_params, collapse = ", "))
      target_params <- intersect(target_params, param_cols)
    }
  }
  
  # Determine the name of the impute column
  impute_col <- if (!is.null(impute_column)) {
    if (!impute_column %in% colnames(data)) {
      stop("The 'data' object does not contain the specified impute column.")
    }
    impute_column
  } else if ("impute" %in% colnames(data)) {
    "impute"
  } else {
    data$impute <- NA_character_
    "impute"
  }
  
  # Identify the targeted intervals based on the groups
  if (!is.null(target_groups)) {
    target_intervals <- dplyr::inner_join(data, target_groups, by = names(target_groups))
  } else {
    target_intervals <- data
  }
  
  # Identify the targeted intervals based on the parameters
  target_intervals <- target_intervals %>%
    dplyr::filter(rowSums(dplyr::across(dplyr::any_of(target_params), ~ . == TRUE)) > 0)
  
  # Add the imputation method to the targeted intervals
  new_intervals_with_impute <- target_intervals %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(param_cols), ~FALSE)) %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(target_params), ~TRUE)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(!!impute_col := {
      impute_methods <- unlist(strsplit(ifelse(is.na(.data[[impute_col]]), "", .data[[impute_col]]), ","))
      if (!allow_duplication && target_impute %in% impute_methods) {
        # If duplication is not allowed, do not add the impute method if it already exists
        .data[[impute_col]]
      } else {
        # Add the impute method after the specified position
        impute_methods <- append(impute_methods, target_impute, after)
        paste(impute_methods[impute_methods != ""], collapse = ",")
      }
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(index = if (new_rows_after_original) index + 0.5 else index + max(index)) %>%
    as.data.frame()
  
  # Eliminate from the old intervals the target parameters
  old_intervals_without_impute <- target_intervals %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(target_params), ~FALSE))
  
  # Make parameters FALSE in original intervals and join the new ones
  data <- data %>%
    dplyr::anti_join(target_intervals, by = names(data)) %>%
    dplyr::bind_rows(old_intervals_without_impute, new_intervals_with_impute) %>%
    dplyr::filter(rowSums(dplyr::across(dplyr::any_of(param_cols), as.numeric)) > 0) %>%
    dplyr::arrange(index) %>%
    dplyr::select(-index)
  
  data
}


#' Remove specified imputation methods from the intervals in a PKNCAdata object.
#'
#' @inheritParams interval_add_impute
#' @param target_impute A character string specifying the imputation method to be removed.
#' @return A modified PKNCAdata object with the specified imputation methods removed from the targeted intervals.
#' @examples
#' d_conc <- data.frame(
#'   conc = c(1, 0.6, 0.2, 0.1, 0.9, 0.4, 1.2, 0.8, 0.3, 0.2, 1.1, 0.5),
#'   time = rep(0:5, 2),
#'   ID = rep(1:2, each = 6),
#'   analyte = rep(c("Analyte1", "Analyte2"), each = 6),
#'   include_hl = c(FALSE, NA, TRUE, TRUE, TRUE, TRUE, FALSE, NA, TRUE, TRUE, TRUE, TRUE)
#' )
#'
#' d_dose <- data.frame(
#'   dose = c(100, 200),
#'   time = c(0, 0),
#'   ID = c(1, 2)
#' )
#'
#' o_conc <- PKNCAconc(d_conc, conc ~ time | ID / analyte)
#' o_dose <- PKNCAdose(d_dose, dose ~ time | ID)
#'
#' intervals <- data.frame(
#'   start = c(0, 0, 0),
#'   end = c(3, 5, Inf),
#'   half.life = c(TRUE, FALSE, TRUE),
#'   cmax = c(TRUE, TRUE, TRUE),
#'   impute = c("start_conc0,start_predose", "start_predose", "start_conc0"),
#'   analyte = c("Analyte1", "Analyte2", "Analyte1")
#' )
#'
#' o_data <- PKNCAdata(o_conc, o_dose, intervals = intervals)
#'
#' # Apply interval_remove_impute function
#' o_data <- interval_remove_impute(data = o_data, target_impute = "start_conc0", target_params = c("half.life"), target_groups = data.frame(analyte = "Analyte1"))
#' @export
interval_remove_impute <- function(data, target_impute, target_params = NULL, target_groups = NULL, impute_column = NULL, new_rows_after_original = TRUE) {
  if (missing(data) || missing(target_impute)) {
    stop("Both 'data' and 'target_impute' must be provided.")
  }
  if (!inherits(data, "PKNCAdata") && !is.data.frame(data)) {
    stop("The 'data' object must be a PKNCAdata object or a data frame.")
  }
  if (!is.character(target_impute)) {
    stop("'target_impute' must be a character string.")
  }
  UseMethod("interval_remove_impute")
}

#' @export
interval_remove_impute.PKNCAdata <- function(data, target_impute, target_params = NULL, target_groups = NULL, impute_column = NULL, new_rows_after_original = TRUE) {
  if (is.null(impute_column) && !is.na(data$impute)) {
    impute_column <- data$impute
  }
  data$intervals <- interval_remove_impute(data$intervals, target_impute, target_params, target_groups, impute_column, new_rows_after_original)
  data
}

#' @export
interval_remove_impute.data.frame <- function(data, target_impute, target_params = NULL, target_groups = NULL, impute_column = NULL, new_rows_after_original = TRUE) {
  
  # Add an index column to preserve the original order
  data <- dplyr::mutate(data, index = dplyr::row_number())
  
  # Get all parameter column names in the data frame
  all_param_options <- setdiff(names(get.interval.cols()), c("start", "end"))
  param_cols <- intersect(names(data), all_param_options)
  
  # Handle target_params
  if (is.null(target_params)) {
    # Take all parameter columns present with at least one TRUE value
    target_params <- param_cols[colSums(data[param_cols]) > 0]
  } else {
    # Check that all target_params are present in the data frame
    missing_params <- setdiff(target_params, param_cols)
    if (length(missing_params) > 0) {
      stop("The following target_params are not interval columns and/or known PKNCA parameters: ", paste(missing_params, collapse = ", "))
      target_params <- intersect(target_params, param_cols)
    }
  }
  
  # Determine the name of the impute column
  impute_col <- if (!is.null(impute_column)) {
    if (!impute_column %in% colnames(data)) {
      stop("The 'data' object does not contain the specified impute column.")
    }
    impute_column
  } else if ("impute" %in% colnames(data)) {
    "impute"
  } else {
    warning("No default impute column identified. No impute methods to remove. If there is an impute column, please specify it in argument 'impute_column'")
    return(data %>% dplyr::select(-index))
  }
  
  # Identify the targeted intervals based on the groups
  if (!is.null(target_groups)) {
    target_intervals <- dplyr::inner_join(data, target_groups, by = names(target_groups))
  } else {
    target_intervals <- data
  }
  
  # Identify the targeted intervals based on the impute method and parameters
  target_intervals <- target_intervals %>%
    dplyr::filter(rowSums(dplyr::across(dplyr::any_of(target_params), ~ . == TRUE)) > 0) %>%
    dplyr::filter(grepl(
      pattern = paste0(".*(", paste0(target_impute, collapse = ")|("), ").*"),
      .data[[impute_col]]
    ))
  
  # Create the new version intervals only for the target parameters
  new_intervals_without_impute <- target_intervals %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(param_cols), ~FALSE)) %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(target_params), ~TRUE)) %>%
    dplyr::rowwise() %>%
    # Eliminate the target impute method from the impute column
    dplyr::mutate(!!impute_col := paste0(setdiff(unlist(strsplit(.data[[impute_col]], ",")), target_impute),
                                         collapse = ","
    )) %>%
    dplyr::mutate(!!impute_col := ifelse(.data[[impute_col]] == "", NA_character_, .data[[impute_col]])) %>%
    dplyr::ungroup() %>%
    # Make sure the class of the impute_col remains the same
    dplyr::mutate(!!impute_col := as.character(.data[[impute_col]])) %>%
    as.data.frame()
  
  # Eliminate from the old intervals the target parameters
  old_intervals_with_impute <- target_intervals %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(target_params), ~FALSE)) %>%
    dplyr::mutate(index = if (new_rows_after_original) index + 0.5 else index + max(index))
  
  # Make parameters FALSE in original intervals and join the new ones
  data <- data %>%
    dplyr::anti_join(target_intervals, by = names(data)) %>%
    dplyr::bind_rows(old_intervals_with_impute, new_intervals_without_impute) %>%
    dplyr::filter(rowSums(dplyr::across(dplyr::any_of(param_cols), as.numeric)) > 0) %>%
    dplyr::arrange(index) %>%
    dplyr::select(-index)
  
  data
}

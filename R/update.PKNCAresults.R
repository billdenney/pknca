#' Update existing PKNCAresults with new data
#'
#' The only thing that can change in the update is the concentration or dose
#' data. All other items (column definitions, etc) must remain the same. If
#' options are not set in `data`, then the default PKNCA options will be
#' considered identical.
#'
#' If more than the allowed settings change, then a full recalculation will
#' occur.
#'
#' This function is typically used with Shiny apps which may repeat analyses
#' with small changes (e.g. point exclusion).
#'
#' @param object The PKNCAresults data
#' @param data The new PKNCAdata object
#' @param ... Ignored
#' @returns The PKNCAresults with updated new data
#' @export
update.PKNCAresults <- function(object, data, ...) {
  # Ensure that the right types of inputs are given
  assert_PKNCAresults(object)
  assert_PKNCAdata(data)
  if (length(data$options) == 0) {
    data$options <- PKNCA.options()
  }
  if (identical(as_PKNCAdata(object), data)) {
    message("No changes detected in data")
    return(object)
  }
  if (!identical(strip_source_data(as_PKNCAdata(object)), strip_source_data(data))) {
    warning("Full recalculation: changes detected in data other than source concentration or dose data")
    return(pk.nca(data))
  }
  # detect changed groups
  groups_changed <- find_changed_group(old = as_PKNCAdata(object), new = data)
  conc_changed <- filter_changed(as.data.frame(as_PKNCAconc(data)), changed = groups_changed)
  dose_changed <- filter_changed(as.data.frame(as_PKNCAdose(data)), changed = groups_changed)
  # insert the changed data into the old data as new data!
  data_new <- data
  data_new$conc$data <- conc_changed
  data_new$dose$data <- dose_changed
  result_new <- pk.nca(data_new)
  result_new_df <- as.data.frame(result_new)
  result_old_df <- as.data.frame(object)
  drop_old_groups <- unique(getGroups(result_new))
  result_old_df_keep <- dplyr::anti_join(result_old_df, drop_old_groups, by = names(drop_old_groups))
  ret <- object
  ret$data <- rbind(result_old_df_keep, result_new_df)
  addProvenance(ret, replace = TRUE)
}

# remove the original data.frames from the source data to enable comparison for
# updates
strip_source_data <- function(data) {
  assert_PKNCAdata(data)
  ret <- data
  # These are no longer valid PKNCAconc or PKNCAdose objects
  ret$conc$data <- NULL
  ret$dose$data <- NULL
  ret
}

#' Find subject identifiers that have changes
#' @param old,new Two PKNCAconc, PKNCAdose, or PKNCAdata objects (must be the
#'   same class)
#' @returns A data.frame of groups that have changed (PKNCAconc or PKNCAdose) or
#'   a list of data.frames (PKNCAdata)
#' @noRd
find_changed_group <- function(old, new) {
  stopifnot(all(class(old) == class(new)))
  if (inherits(old, "PKNCAdata")) {
    # Find subjects that changed (for PKNCAdata by going into conc and dose)
    list(
      conc = find_changed_group(old = as_PKNCAconc(old), new = as_PKNCAconc(new)),
      dose = find_changed_group(old = as_PKNCAdose(old), new = as_PKNCAdose(new))
    )
  } else {
    # Find subjects that changed (for PKNCAconc or PKNCAdose)
    group_col <- unlist(old$columns$groups, use.names = FALSE)
    d_nest_old <- tidyr::nest(old$data, data_old = !tidyr::all_of(group_col))
    d_nest_new <- tidyr::nest(new$data, data_new = !tidyr::all_of(group_col))
    d_nest_combo <- dplyr::full_join(d_nest_old, d_nest_new, by = group_col)
    mask_changed_id <-
      vapply(
        X = seq_len(nrow(d_nest_combo)),
        FUN = function(idx) !identical(d_nest_combo$data_old[[idx]], d_nest_combo$data_new[[idx]]),
        FUN.VALUE = TRUE
      )
    # Find the unique set of groups that changed
    unique(d_nest_combo[mask_changed_id, group_col])
  }
}

filter_changed <- function(data, changed) {
  rowid_col <- paste0(max(names(data)), "X")
  tracking <- data
  tracking[[rowid_col]] <- seq_len(nrow(data))
  # Find rows that are of interest from concentration or dose
  tracking_conc <- filter_changed_inner_join(tracking, changed$conc)
  tracking_dose <- filter_changed_inner_join(tracking, changed$dose)
  # filter to rows of interest
  ret <- tracking[tracking[[rowid_col]] %in% c(tracking_conc[[rowid_col]], tracking_dose[[rowid_col]]), ]
  ret[, names(data)]
}

filter_changed_inner_join <- function(data, changed) {
  if (nrow(changed) == 0) {
    # Return a zero-row data.frame if nothing changed
    data[0,]
  } else if (length(intersect(names(data), names(changed))) == 0) {
    # Return all the data if there is not an interesction in column names
    data
  } else {
    dplyr::inner_join(data, changed, by = intersect(names(data), names(changed)))
  }
}

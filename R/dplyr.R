#' @importFrom dplyr inner_join
#' @export
dplyr::inner_join
#' @importFrom dplyr left_join
#' @export
dplyr::left_join
#' @importFrom dplyr right_join
#' @export
dplyr::right_join
#' @importFrom dplyr full_join
#' @export
dplyr::full_join

#' @importFrom dplyr filter
#' @export
dplyr::filter

#' @importFrom dplyr group_by
#' @export
dplyr::group_by

#' @importFrom dplyr ungroup
#' @export
dplyr::ungroup

#' @importFrom dplyr mutate
#' @export
dplyr::mutate

join_maker_PKNCA <- function(join_fun) {
  function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ..., keep = FALSE) { # nocov
    dataname <- getDataName(x)
    x[[dataname]] <- join_fun(x=x[[dataname]], y=y, by = by, copy = copy, suffix = suffix, ..., keep = keep)
    x
  }
}
filter_PKNCA <- function(.data, ..., .preserve=FALSE) {
  dataname <- getDataName(.data)
  .data[[dataname]] <- dplyr::filter(.data[[dataname]], ..., .preserve=.preserve)
  .data
}
mutate_PKNCA <- function(.data, ...) {
  dataname <- getDataName(.data)
  .data[[dataname]] <- dplyr::mutate(.data[[dataname]], ...)
  .data
}
group_by_PKNCA <- function(.data, ..., .add = FALSE, .drop = group_by_drop_default(.data)) {
  dataname <- getDataName(.data)
  .data[[dataname]] <- dplyr::group_by(.data[[dataname]], ..., .add = FALSE, .drop = group_by_drop_default(.data))
  .data
}
ungroup_PKNCA <- function(x, ...) {
  dataname <- getDataName(x)
  x[[dataname]] <- dplyr::ungroup(x[[dataname]], ...)
  x
}

#' dplyr joins for PKNCA
#' 
#' @inheritParams dplyr::inner_join
#' @family dplyr verbs
#' @export
inner_join.PKNCAresults <- join_maker_PKNCA(dplyr::inner_join)
#' @rdname inner_join.PKNCAresults
#' @export
left_join.PKNCAresults <- join_maker_PKNCA(dplyr::left_join)
#' @rdname inner_join.PKNCAresults
#' @export
right_join.PKNCAresults <- join_maker_PKNCA(dplyr::right_join)
#' @rdname inner_join.PKNCAresults
#' @export
full_join.PKNCAresults <- join_maker_PKNCA(dplyr::full_join)

#' @rdname inner_join.PKNCAresults
#' @export
inner_join.PKNCAconc <- join_maker_PKNCA(dplyr::inner_join)
#' @rdname inner_join.PKNCAresults
#' @export
left_join.PKNCAconc <- join_maker_PKNCA(dplyr::left_join)
#' @rdname inner_join.PKNCAresults
#' @export
right_join.PKNCAconc <- join_maker_PKNCA(dplyr::right_join)
#' @rdname inner_join.PKNCAresults
#' @export
full_join.PKNCAconc <- join_maker_PKNCA(dplyr::full_join)

#' @rdname inner_join.PKNCAresults
#' @export
inner_join.PKNCAdose <- join_maker_PKNCA(dplyr::inner_join)
#' @rdname inner_join.PKNCAresults
#' @export
left_join.PKNCAdose <- join_maker_PKNCA(dplyr::left_join)
#' @rdname inner_join.PKNCAresults
#' @export
right_join.PKNCAdose <- join_maker_PKNCA(dplyr::right_join)
#' @rdname inner_join.PKNCAresults
#' @export
full_join.PKNCAdose <- join_maker_PKNCA(dplyr::full_join)

#' dplyr filtering for PKNCA
#'
#' @inheritParams dplyr::filter
#' @family dplyr verbs
#' @export
filter.PKNCAresults <- filter_PKNCA
#' @rdname filter.PKNCAresults
#' @export
filter.PKNCAconc <- filter_PKNCA
#' @rdname filter.PKNCAresults
#' @export
filter.PKNCAdose <- filter_PKNCA

#' dplyr mutate-based modification for PKNCA
#' 
#' @inheritParams dplyr::mutate
#' @family dplyr verbs
#' @export
mutate.PKNCAresults <- mutate_PKNCA
#' @rdname mutate.PKNCAresults
#' @export
mutate.PKNCAconc <- mutate_PKNCA
#' @rdname mutate.PKNCAresults
#' @export
mutate.PKNCAdose <- mutate_PKNCA

#' dplyr grouping for PKNCA
#' 
#' @inheritParams dplyr::group_by
#' @family dplyr verbs
#' @export
group_by.PKNCAresults <- group_by_PKNCA
#' @rdname group_by.PKNCAresults
#' @export
group_by.PKNCAconc <- group_by_PKNCA
#' @rdname group_by.PKNCAresults
#' @export
group_by.PKNCAdose <- group_by_PKNCA
#' @rdname group_by.PKNCAresults
#' @inheritParams dplyr::ungroup
#' @export
ungroup.PKNCAresults <- ungroup_PKNCA
#' @rdname group_by.PKNCAresults
#' @export
ungroup.PKNCAconc <- ungroup_PKNCA
#' @rdname group_by.PKNCAresults
#' @export
ungroup.PKNCAdose <- ungroup_PKNCA

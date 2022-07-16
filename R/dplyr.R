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

#' dplyr filtering for PKNCA
#'
#' @inheritParams dplyr::filter
#' @family dplyr verbs
#' @export
filter.PKNCAresults <- filter_PKNCA
#' @rdname filter.PKNCAresults
#' @export
filter.PKNCAconc <- filter_PKNCA

#' dplyr mutate-based modification for PKNCA
#' 
#' @inheritParams dplyr::mutate
#' @family dplyr verbs
#' @export
mutate.PKNCAresults <- mutate_PKNCA
#' @rdname mutate.PKNCAresults
#' @export
mutate.PKNCAconc <- mutate_PKNCA

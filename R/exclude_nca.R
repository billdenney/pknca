#' Exclude NCA parameters based on examining the parameter set.
#'
#' @param min.span.ratio The minimum acceptable span ratio (uses
#'   `PKNCA.options("min.span.ratio")` if not provided).
#' @param max.aucinf.pext The maximum acceptable percent AUC extrapolation (uses
#'   `PKNCA.options("max.aucinf.pext")` if not provided).
#' @param min.hl.r.squared The minimum acceptable r-squared value for half-life
#'   (uses `PKNCA.options("min.hl.r.squared")` if not provided).
#' @examples
#' my_conc <- PKNCAconc(data.frame(conc=1.1^(3:0),
#'                                 time=0:3,
#'                                 subject=1),
#'                      conc~time|subject)
#' my_data <- PKNCAdata(my_conc,
#'                      intervals=data.frame(start=0, end=Inf,
#'                                           aucinf.obs=TRUE,
#'                                           aucpext.obs=TRUE))
#' my_result <- pk.nca(my_data)
#' my_result_excluded <- exclude(my_result,
#'                               FUN=exclude_nca_max.aucinf.pext())
#' as.data.frame(my_result_excluded)
#' @name exclude_nca
#' @family Result exclusions
NULL

#' @describeIn exclude_nca Exclude based on span.ratio
#' @export
exclude_nca_span.ratio <- function(min.span.ratio) {
  affected_parameters <- get.parameter.deps("half.life")
  missing_min.span.ratio <- missing(min.span.ratio)
  function(x, ...) {
    if (missing_min.span.ratio) {
      min.span.ratio <- PKNCA.options("min.span.ratio")
    }
    ret <- rep(NA_character_, nrow(x))
    if (!is.na(min.span.ratio)) {
      idx_span.ratio <- which(x$PPTESTCD %in% "span.ratio")
      if (length(idx_span.ratio) == 0) {
        # Do nothing, it wasn't calculated
      } else if (length(idx_span.ratio) == 1) {
        current_span.ratio <- x$PPORRES[idx_span.ratio]
        drop_span_ratio <-
          !is.na(current_span.ratio) &
          current_span.ratio < min.span.ratio
        if (drop_span_ratio) {
          ret[x$PPTESTCD %in% affected_parameters] <-
            sprintf("Span ratio < %g", min.span.ratio)
        }
      } else if (length(idx_span.ratio) > 1) { # nocov
        stop("Should not see more than one span.ratio (please report this as a bug)") # nocov
      }
    }
    ret
  }
}

#' @describeIn exclude_nca Exclude based on AUC percent extrapolated (both
#'   observed and predicted)
#' @export
exclude_nca_max.aucinf.pext <-  function(max.aucinf.pext) {
  affected_parameters <-
    list(obs=get.parameter.deps("aucinf.obs"),
         pred=get.parameter.deps("aucinf.pred"))
  missing_max.aucinf.pext <- missing(max.aucinf.pext)
  function(x, ...) {
    if (missing_max.aucinf.pext) {
      max.aucinf.pext <- PKNCA.options("max.aucinf.pext")
    }
    ret <- rep(NA_character_, nrow(x))
    if (!is.na(max.aucinf.pext)) {
      for (ext_type in c("obs", "pred")) {
        idx_pext <- which(x$PPTESTCD %in% paste0("aucpext.", ext_type))
        if (length(idx_pext) == 0) {
          # Do nothing, it wasn't calculated
        } else if (length(idx_pext) == 1) {
          current_pext <- x$PPORRES[idx_pext]
          drop_pext <-
            !is.na(current_pext) &
            current_pext > max.aucinf.pext
          if (drop_pext) {
            ret[x$PPTESTCD %in% affected_parameters[[ext_type]]] <-
              sprintf("AUC percent extrapolated > %g", max.aucinf.pext)
          }
        } else if (length(idx_pext) > 1) { # nocov
          stop("Should not see more than one aucpext.", ext_type, " (please report this as a bug)") # nocov
        }
      }
    }
    ret
  }
}

#' @describeIn exclude_nca Exclude AUC measurements based on count of
#'   concentrations measured and not below the lower limit of quantification
#' @param min_count Minimum number of measured concentrations
#' @param exclude_param_pattern Character vector of regular expression patterns
#'   to exclude
#' @export
exclude_nca_count_conc_measured <-  function(min_count, exclude_param_pattern = c("^aucall", "^aucinf", "^aucint", "^auciv", "^auclast", "^aumc", "^sparse_auc")) {
  all_parameters <- names(PKNCA::get.interval.cols())
  affected_parameters_base <-
    sort(unique(unlist(
      lapply(
        X = exclude_param_pattern,
        FUN = grep,
        x = all_parameters,
        value = TRUE
      )
    )))
  affected_parameters <-
    sort(unique(unlist(
      lapply(
        X = affected_parameters_base,
        FUN = get.parameter.deps
      )
    )))
  force(min_count)
  function(x, ...) {
    ret <- rep(NA_character_, nrow(x))
    if (!is.na(min_count)) {
      idx_count <- which(x$PPTESTCD %in% "count_conc_measured")
      if (length(idx_count) == 0) {
        # Do nothing, it wasn't calculated
      } else if (length(idx_count) == 1) {
        current_count <- x$PPORRES[idx_count]
        drop_count <-
          !is.na(current_count) &
          current_count < min_count
        if (drop_count) {
          ret[x$PPTESTCD %in% affected_parameters] <-
            sprintf("Number of measured concentrations is < %g", min_count)
        }
      } else { # nocov
        stop("Should not see more than one count_conc_measured (please report this as a bug)") # nocov
      }
    }
    ret
  }
}

#' @describeIn exclude_nca Exclude based on half-life r-squared
#' @export
exclude_nca_min.hl.r.squared <- function(min.hl.r.squared) {
  affected_parameters <- get.parameter.deps("half.life")
  missing_min.hl.r.squared <- missing(min.hl.r.squared)
  function(x, ...) {
    if (missing_min.hl.r.squared) {
      min.hl.r.squared <- PKNCA.options("min.hl.r.squared")
    }
    ret <- rep(NA_character_, nrow(x))
    if (!is.na(min.hl.r.squared)) {
      idx_r.squared <- which(x$PPTESTCD %in% "r.squared")
      if (length(idx_r.squared) == 0) {
        # Do nothing, it wasn't calculated
      } else if (length(idx_r.squared) == 1) {
        current_r.squared <- x$PPORRES[idx_r.squared]
        drop_r.squared <-
          !is.na(current_r.squared) &
          current_r.squared < min.hl.r.squared
        if (drop_r.squared) {
          ret[x$PPTESTCD %in% affected_parameters] <-
            sprintf("Half-life r-squared < %g", min.hl.r.squared)
        }
      } else if (length(idx_r.squared) > 1) { # nocov
        stop("Should not see more than one r.squared (please report this as a bug)") # nocov
      }
    }
    ret
  }
}

#' @describeIn exclude_nca Exclude based on implausibly early Tmax (often used
#'   for extravascular dosing with a Tmax value of 0)
#' @param tmax_early The time for Tmax which is considered too early to be a
#'   valid NCA result
#' @export
exclude_nca_tmax_early <- function(tmax_early = 0) {
  force(tmax_early)
  function(x, ...) {
    ret <- rep(NA_character_, nrow(x))
    idx_tmax <- which(x$PPTESTCD %in% "tmax")
    if (length(idx_tmax) == 1) {
      current_tmax <- x$PPORRES[idx_tmax]
      drop <- !is.na(current_tmax) & current_tmax <= tmax_early
      if (drop) {
        ret <- rep(sprintf("Tmax is <=%g (likely missed dose, insufficient PK samples, or PK sample swap)", tmax_early), nrow(x))
      }
    } else if (length(idx_tmax) != 0) {
      stop("Should not see more than one tmax (please report this as a bug)") # nocov
    }
    ret
  }
}

#' @describeIn exclude_nca Exclude based on implausibly early Tmax (special case
#'   for `tmax_early = 0`)
#' @export
exclude_nca_tmax_0 <- function() {
  exclude_nca_tmax_early()
}

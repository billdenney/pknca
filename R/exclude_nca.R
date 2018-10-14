#' Exclude NCA parameters based on examining the parameter set.
#' 
#' @param min.span.ratio The minimum acceptable span ratio (uses 
#'   \code{PKNCA.options("min.span.ratio")} if not provided).
#' @param max.aucinf.pext The maximum acceptable percent AUC
#'   extrapolation (uses \code{PKNCA.options("max.aucinf.pext")} if not
#'   provided).
#' @param min.hl.r.squared The minimum acceptable r-squared value for
#'   half-life (uses \code{PKNCA.options("min.hl.r.squared")} if not
#'   provided).
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

#' @describeIn exclude_nca Exclude based on AUC percent extrapolated
#'   (both observed and predicted)
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

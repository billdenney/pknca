#' Assess the AIC for all models in a list of models
#'
#' @param x the list of models
#' @param \dots parameters passed to the underlying AIC function
#' (typically the parameter k)
#' @param assess.best determine which model is the best (by lowest
#' AIC)
#' @return a data frame with row names matching the names of the list
#' \code{x} and columns for degrees of freedom (\code{df}) and
#' \code{AIC}.  If \code{assess.best} is true, then there will be
#' another column \code{isBest}.
#' @export
AIC.list <- function(x, ..., assess.best=TRUE) {
  ret <- data.frame()
  for (i in 1:length(x)) {
    if (identical(x[[i]], NA)) {
      tmp <- data.frame(df=NA, AIC=NA)
    } else {
      tmp <- AIC(x[[1]], x[[i]], ...)
    }
    ret <- rbind(ret, tmp[2,])
  }
  rownames(ret) <- names(x)
  if (assess.best) {
    ret$isBest <- ""
    ret$isBest[ret$AIC %in% min(ret$AIC, na.rm=TRUE)] <- "Best Model"
  }
  ret
}

get.best.model <- function(x, ...)
  x[AIC(x, ...)$isBest %in% "Best Model"][[1]]

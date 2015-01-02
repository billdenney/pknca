#' Assess the AIC for all models in a list of models
#'
#' @param object the list of models
#' @param \dots parameters passed to the underlying AIC function
#' (typically the parameter k)
#' @param assess.best determine which model is the best (by lowest
#' AIC)
#' @seealso \code{\link{get.best.model}}
#' @return a data frame with row names matching the names of the list
#' \code{x} and columns for degrees of freedom (\code{df}) and
#' \code{AIC}.  If \code{assess.best} is true, then there will be
#' another column \code{isBest}.
#' @export
AIC.list <- function(object, ..., assess.best=TRUE) {
  ret <- data.frame()
  for (i in 1:length(object)) {
    if (identical(object[[i]], NA)) {
      tmp <- data.frame(df=NA, AIC=NA)
    } else {
      tmp <- AIC(object[[1]], object[[i]], ...)
    }
    ret <- rbind(ret, tmp[2,])
  }
  rownames(ret) <- names(object)
  if (assess.best) {
    ret$isBest <- ""
    ret$isBest[ret$AIC %in% min(ret$AIC, na.rm=TRUE)] <- "Best Model"
  }
  ret
}

#' Extract the best model from a list of models using AIC.list.
#'
#' @param object the list of models
#' @param \dots Parameters passed to AIC.list
#' @seealso \code{\link{AIC.list}}
#' @return The model which is assessed as best.  If more than one are
#' equal, the first is chosen.
#' @export
get.best.model <- function(object, ...)
  object[AIC(object, ...)$isBest %in% "Best Model"][[1]]

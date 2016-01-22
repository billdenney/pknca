#' Compute noncompartmental pharmacokinetics
#'
#' This package computes pharmacokinetic (PK) noncompartmental
#' analysis (NCA) parameters.
#'
#' Typically, you will start with the \code{\link{pk.nca}}
#' function and then use the \code{\link{summary}} function to
#' generate tables for reporting.
#'
#' @docType package
#' @name PKNCA
#' @importFrom digest digest
#' @importFrom doBy renameCol splitBy summaryBy
#' @importFrom graphics plot
#' @importFrom lattice xyplot
#' @importFrom nlme fixef getData getGroups gnls intervals lme nlme ranef
#' @importFrom parallel mclapply
#' @importFrom plyr rbind.fill
#' @importFrom stats AIC as.formula coef confint formula glm logLik lm median model.frame na.exclude na.omit predict sd update update.formula
#' @importFrom utils head sessionInfo
NULL

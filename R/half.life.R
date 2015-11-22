#' Compute the half-life and associated parameters
#'
#' The half-life is calculated by computing the best fit line for all
#' available sets of points.  The best one is chosen by the following
#' rules in order:
#'
#' \itemize{
#' \item{At least \code{min.hl.points} points included}
#' \item{A \code{lambda.z} > 0}
#' \item{The best adjusted r-squared (within \code{adj.r.squared.factor})}
#' \item{The one with the most points included}
#' }
#'
#' @param conc Concentration measured
#' @param time Time of concentration measurement
#' @param tmax Time of maximum concentration (will be calculated and
#' included in the return data frame if not given)
#' @param tlast Time of last concentration above the limit of
#' quantification (will be calculated and included in the return data
#' frame if not given)
#' @param options List of changes to the default
#' \code{\link{PKNCA.options}} for calculations.
#' @param min.hl.points The minimum number of points that must be
#' included to calculate the half-life
#' @param adj.r.squared.factor The allowance in adjusted r-squared for
#' adding another point.
#' @param conc.blq See \code{\link{clean.conc.blq}}
#' @param conc.na See \code{\link{clean.conc.na}}
#' @param check Run \code{\link{check.conc.time}},
#' \code{\link{clean.conc.blq}}, and \code{\link{clean.conc.na}}?
#' @param first.tmax See \code{\link{pk.calc.tmax}}.
#' @param allow.tmax.in.half.life Allow the concentration point for
#' tmax to be included in the half-life slope calculation.
#' @return A data frame with one row and columns for
#' \describe{
#'   \item{tmax}{Time of maximum observed concentration (only included if not given as an input)}
#'   \item{tlast}{Time of last observed concentration above the LOQ (only included if not given as an input)}
#'   \item{r.squared}{coefficient of determination}
#'   \item{adj.r.squared}{adjusted coefficient of determination}
#'   \item{lambda.z}{elimination rate}
#'   \item{lambda.z.time.first}{first time for half-life calculation}
#'   \item{lambda.z.n.points}{number of points in half-life calculation}
#'   \item{clast.pred}{Concentration at tlast as predicted by the half-life
#'       line}
#'   \item{half.life}{half-life}
#'   \item{span.ratio}{span ratio [ratio of half-life to time used for half-life
#'       calculation}
#' }
#' @references
#'
#' Gabrielsson J, Weiner D.  "Section 2.8.4 Strategies for estimation
#' of lambda-z."  Pharmacokinetic & Pharmacodynamic Data Analysis:
#' Concepts and Applications, 4th Edition.  Stockholm, Sweden: Swedish
#' Pharmaceutical Press, 2000.  167-9.
#' @export
pk.calc.half.life <- function(conc, time, tmax, tlast,
                              options=list(),
                              min.hl.points=PKNCA.choose.option("min.hl.points", options),
                              adj.r.squared.factor=PKNCA.choose.option("adj.r.squared.factor", options),
                              conc.blq=PKNCA.choose.option("conc.blq", options),
                              conc.na=PKNCA.choose.option("conc.na", options),
                              first.tmax=PKNCA.choose.option("first.tmax", options),
                              allow.tmax.in.half.life=PKNCA.choose.option("allow.tmax.in.half.life", options),
                              check=TRUE) {
  ## Check inputs
  min.hl.points <-
    PKNCA.options(min.hl.points=min.hl.points, check=TRUE)
  adj.r.squared.factor <-
    PKNCA.options(adj.r.squared.factor=adj.r.squared.factor,
                  check=TRUE)
  if (check) {
    check.conc.time(conc, time)
    data <- clean.conc.blq(conc, time, conc.blq=conc.blq, conc.na=conc.na)
  } else {
    data <- data.frame(conc, time)
  }
  ## Prepare the return values
  ret <- data.frame(
    ## Terminal elimination slope
    lambda.z=NA,
    ## R-squared of terminal elimination slope
    r.squared=NA,
    ## Adjusted r-squared of terminal elimination slope
    adj.r.squared=NA,
    ## First time point used in the slope estimation
    ## (for plotting later)
    lambda.z.time.first=NA,
    ## Number of points in the half-life estimate
    lambda.z.n.points=NA,
    ## Concentration at Tlast predicted by the half-life
    clast.pred=NA,
    ## Half-life
    half.life=NA,
    ## T1/2 span range
    span.ratio=NA)
  if (missing(tmax)) {
    ret$tmax <- pk.calc.tmax(data$conc, data$time,
                             first.tmax=first.tmax, check=FALSE)
  } else {
    ret$tmax <- tmax
  }
  if (missing(tlast)) {
    ret$tlast <- pk.calc.tlast(data$conc, data$time, check=FALSE)
  } else {
    ret$tlast <- tlast
  }
  ## Data frame to use for computation of half-life
  if (allow.tmax.in.half.life) {
    dfK <- data[data$time >= ret$tmax & data$conc != 0,]
  } else {
    dfK <- data[data$time > ret$tmax & data$conc != 0,]
  }
  dfK$logDV <- log(dfK$conc)
  if (nrow(dfK) >= min.hl.points) {
    ## If we have enough data to estimate a slope, then
    EQP <- data.frame(r.squared=0,
                      adj.r.squared=0,
                      PROB=NA,
                      clast.pred=NA,
                      lambda.z=NA,
                      lambda.z.n.points=NA,
                      lambda.z.time.first=dfK$time,
                      DV=dfK$conc)
    EQP <- EQP[order(-EQP$lambda.z.time.first),]
    for(i in min.hl.points:nrow(EQP)) {
      ## Fit the terminal slopes until the adjusted r-squared value
      ## is not improving (or it only gets worse by a small factor).
      DF2 <- data.frame(Y=log(EQP$DV[1:i]),X=EQP$lambda.z.time.first[1:i])
      fit <- stats::lm(Y~X, data=DF2, na.action=stats::na.exclude)
      sfit <- summary(fit)
      EQP$r.squared[i] <- sfit$r.squared
      EQP$adj.r.squared[i] <- adj.r.squared(sfit$r.squared, i)
      EQP$PROB[i] <- sfit$coefficients["X", "Pr(>|t|)"]
      EQP$lambda.z[i] <- -stats::coef(fit)["X"]
      EQP$clast.pred[i] <- exp(stats::predict(fit, newdata=data.frame(X=ret$tlast)))
      EQP$lambda.z.n.points[i] <- i
    }
    ## Find the best model
    mask.best <-
      (EQP$adj.r.squared > (max(EQP$adj.r.squared) - adj.r.squared.factor) &
       EQP$lambda.z > 0)
    ## Missing values are not the best
    mask.best[is.na(mask.best)] <- FALSE
    if (sum(mask.best) > 1) {
      ## If more than one models qualify, choose the one with the
      ## most points used.
      mask.best <- (mask.best &
                    EQP$lambda.z.n.points == max(EQP$lambda.z.n.points[mask.best]))
    }
    ## If the half-life fit, set all associated parameters
    if (any(mask.best)) {
      ## Put in all the computed values
      replacements <- c("lambda.z", "r.squared", "adj.r.squared", "lambda.z.time.first",
                        "lambda.z.n.points", "clast.pred")
      ret[,replacements] <- EQP[mask.best,replacements]
      ## Compute the half-life and span ratio
      ret$half.life <- log(2)/ret$lambda.z
      ret$span.ratio <- (ret$tlast - ret$lambda.z.time.first)/ret$half.life
    }
  } else {
    warning(sprintf(
      "Too few points for half-life calculation (min.hl.points=%g with only %g points)",
      min.hl.points, nrow(dfK)))
  }
  ## Drop the inputs of tmax and tlast, if given.
  if (!missing(tmax))
    ret$tmax <- NULL
  if (!missing(tlast))
    ret$tlast <- NULL
  ret
}

## Add the column to the interval specification
add.interval.col("half.life",
                 FUN="pk.calc.half.life",
                 values=c(FALSE, TRUE),
                 desc="The (terminal) half-life",
                 depends=c("tmax", "tlast"))
PKNCA.set.summary("half.life", business.mean, business.sd)
add.interval.col("r.squared",
                 FUN=NA,
                 values=c(FALSE, TRUE),
                 desc="The r^2 value of the half-life calculation",
                 depends=c("half.life"))
add.interval.col("adj.r.squared",
                 FUN=NA,
                 values=c(FALSE, TRUE),
                 desc="The adjusted r^2 value of the half-life calculation",
                 depends=c("half.life"))
add.interval.col("lambda.z",
                 FUN=NA,
                 values=c(FALSE, TRUE),
                 desc="The elimination rate of the terminal half-life",
                 depends=c("half.life"))
PKNCA.set.summary("lambda.z", business.geomean, business.geocv)
add.interval.col("lambda.z.time.first",
                 FUN=NA,
                 values=c(FALSE, TRUE),
                 desc="The first time point used for the calculation of half-life",
                 depends=c("half.life"))
add.interval.col("lambda.z.n.points",
                 FUN=NA,
                 values=c(FALSE, TRUE),
                 desc="The number of points used for the calculation of half-life",
                 depends=c("half.life"))
PKNCA.set.summary("lambda.z.n.points", business.median, business.range)
add.interval.col("clast.pred",
                 FUN=NA,
                 values=c(FALSE, TRUE),
                 desc="The concentration at Tlast as predicted by the half-life",
                 depends=c("half.life"))
PKNCA.set.summary("clast.pred", business.geomean, business.geocv)
add.interval.col("span.ratio",
                 FUN=NA,
                 values=c(FALSE, TRUE),
                 desc="The ratio of the half-life to the duration used for half-life calculation",
                 depends=c("half.life"))
PKNCA.set.summary("span.ratio", business.geomean, business.geocv)

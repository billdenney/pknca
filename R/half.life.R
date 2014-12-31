#' Compute the half-life and associated parameters
#'
#' The half-life is calculated by computing the best fit line for all
#' available sets of points.  The best one is chosen by the following
#' rules in order:
#'
#' \describe{
#' \item At least \code{min.hl.points} points included
#' \item A \code{lambda.z} > 0
#' \item The best adjusted r-squared (within \code{adj.r.squared.factor})
#' \item The one with the most points included
#' }
#'
#' @param conc Concentration measured
#' @param time Time of concentration measurement
#' @param min.hl.points The minimum number of points that must be
#' included to calculate the half-life
#' @param adj.r.squared.factor The allowance in adjusted r-squared for
#' adding another point.
#' @return A data frame with one row and columns for
#' \describe{
#'   \item{tmax}{Time of maximum observed concentration}
#'   \item{tlast}{Time of last observed concentration above the LOQ}
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
pk.calc.half.life <- function(conc, time,
                              min.hl.points=PKNCA.options("min.hl.points"),
                              adj.r.squared.factor=PKNCA.options("adj.r.squared.factor"),
                              conc.blq=PKNCA.options("conc.blq"),
                              conc.na=PKNCA.options("conc.na"),
                              use.first=PKNCA.options("first.tmax"),
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
    tmax=pk.calc.tmax(data$conc, data$time,
      use.first=use.first, check=FALSE),
    tlast=pk.calc.tlast(data$conc, data$time, check=FALSE),
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
  ## Data frame to use for computation of half-life
  dfK <- data[data$time >= ret$tmax & data$conc != 0,]
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
      fit <- lm(Y~X, data=DF2, na.action=na.exclude)
      sfit <- summary(fit)
      EQP$r.squared[i] <- sfit$r.squared
      EQP$adj.r.squared[i] <- adj.r.squared(sfit$r.squared, i)
      EQP$PROB[i] <- sfit$coefficients["X", "Pr(>|t|)"]
      EQP$lambda.z[i] <- -coef(fit)["X"]
      EQP$clast.pred[i] <- exp(predict(fit, newdata=data.frame(X=ret$tlast)))
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
  ret
}

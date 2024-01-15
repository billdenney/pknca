#' Compute the half-life and associated parameters
#'
#' The terminal elimination half-life is estimated from the final points in the
#' concentration-time curve using semi-log regression (`log(conc)~time`)
#' with automated selection of the points for calculation (unless
#' `manually.selected.points` is `TRUE`).
#'
#' See the "Half-Life Calculation" vignette for more details on the calculation
#' methods used.
#'
#' @details If `manually.selected.points` is `FALSE` (default), the
#' half-life is calculated by computing the best fit line for all points at or
#' after tmax (based on the value of `allow.tmax.in.half.life`).  The best
#' half-life is chosen by the following rules in order:
#'
#' \itemize{
#'  \item{At least `min.hl.points` points included}
#'  \item{A `lambda.z` > 0 and at the same time the best adjusted r-squared
#'  (within `adj.r.squared.factor`)}
#'  \item{The one with the most points included}
#' }
#'
#' If `manually.selected.points` is `TRUE`, the `conc`
#' and `time` data are used as-is without any form of selection for
#' the best-fit half-life.
#'
#' @inheritParams assert_conc_time
#' @inheritParams choose_interval_method
#' @inheritParams PKNCA.choose.option
#' @param tmax Time of maximum concentration (will be calculated and
#'   included in the return data frame if not given)
#' @param tlast Time of last concentration above the limit of
#'   quantification (will be calculated and included in the return data
#'   frame if not given)
#' @param manually.selected.points Have the input points (`conc` and
#'   `time`) been manually selected?  The impact of setting this to
#'   `TRUE` is that no selection for the best points will be done.  When
#'   `TRUE`, this option causes the options of `adj.r.squared.factor`,
#'   `min.hl.points`, and `allow.tmax.in.half.life` to be ignored.
#' @param min.hl.points The minimum number of points that must be
#'   included to calculate the half-life
#' @param adj.r.squared.factor The allowance in adjusted r-squared for
#'   adding another point.
#' @param conc.blq See [clean.conc.blq()]
#' @param conc.na See [clean.conc.na()]
#' @param check Run [assert_conc_time()],
#'   [clean.conc.blq()], and [clean.conc.na()]?
#' @param first.tmax See [pk.calc.tmax()].
#' @param allow.tmax.in.half.life Allow the concentration point for tmax
#'   to be included in the half-life slope calculation.
#' @return A data frame with one row and columns for
#'  \describe{
#'   \item{tmax}{Time of maximum observed concentration (only included
#'     if not given as an input)}
#'   \item{tlast}{Time of last observed concentration above the LOQ (only
#'     included if not given as an input)}
#'   \item{r.squared}{coefficient of determination}
#'   \item{adj.r.squared}{adjusted coefficient of determination}
#'   \item{lambda.z}{elimination rate}
#'   \item{lambda.z.time.first}{first time for half-life calculation}
#'   \item{lambda.z.n.points}{number of points in half-life calculation}
#'   \item{clast.pred}{Concentration at tlast as predicted by the half-life
#'     line}
#'   \item{half.life}{half-life}
#'   \item{span.ratio}{span ratio [ratio of half-life to time used for
#'     half-life calculation}
#'  }
#' @references
#'
#' Gabrielsson J, Weiner D.  "Section 2.8.4 Strategies for estimation of
#' lambda-z."  Pharmacokinetic & Pharmacodynamic Data Analysis: Concepts
#' and Applications, 4th Edition.  Stockholm, Sweden: Swedish
#' Pharmaceutical Press, 2000.  167-9.
#' @family NCA parameter calculations
#' @export
pk.calc.half.life <- function(conc, time, tmax, tlast,
                              manually.selected.points=FALSE,
                              options=list(),
                              min.hl.points=NULL,
                              adj.r.squared.factor=NULL,
                              conc.blq=NULL,
                              conc.na=NULL,
                              first.tmax=NULL,
                              allow.tmax.in.half.life=NULL,
                              check=TRUE) {
  # Check inputs
  min.hl.points <-
    PKNCA.choose.option(
      name="min.hl.points", value=min.hl.points, options=options
    )
  adj.r.squared.factor <-
    PKNCA.choose.option(
      name="adj.r.squared.factor", value=adj.r.squared.factor, options=options
    )
  conc.blq <-
    PKNCA.choose.option(
      name="conc.blq", value=conc.blq, options=options
    )
  conc.na <-
    PKNCA.choose.option(
      name="conc.na", value=conc.na, options=options
    )
  first.tmax <-
    PKNCA.choose.option(
      name="first.tmax", value=first.tmax, options=options
    )
  allow.tmax.in.half.life <-
    PKNCA.choose.option(
      name="allow.tmax.in.half.life", value=allow.tmax.in.half.life, options=options
    )
  if (check) {
    assert_conc_time(conc = conc, time = time)
    data <- clean.conc.blq(conc, time, conc.blq=conc.blq, conc.na=conc.na)
  } else {
    data <- data.frame(conc, time)
  }
  # if (inherits(data$conc, "units")) {
  #   conc_units <- units(data$conc)
  # } else {
    conc_units <- NULL
  #}
  data$log_conc <- log(data$conc)
  # as.numeric() to handle units objects
  data <- data[as.numeric(data$conc) > 0,]
  # Prepare the return values
  ret <- data.frame(
    # Terminal elimination slope
    lambda.z=NA,
    # R-squared of terminal elimination slope
    r.squared=NA,
    # Adjusted r-squared of terminal elimination slope
    adj.r.squared=NA,
    # First time point used in the slope estimation (for plotting later)
    lambda.z.time.first=NA,
    # Number of points in the half-life estimate
    lambda.z.n.points=NA,
    # Concentration at Tlast predicted by the half-life
    clast.pred=NA,
    # Half-life
    half.life=NA,
    # T1/2 span ratio
    span.ratio=NA)
  ret_replacements <-
    c("lambda.z", "r.squared", "adj.r.squared", "lambda.z.time.first",
      "lambda.z.n.points", "clast.pred", "half.life", "span.ratio")
  if (missing(tmax)) {
    ret$tmax <-
      pk.calc.tmax(
        data$conc, data$time, first.tmax=first.tmax, check=FALSE
      )
  } else {
    ret$tmax <- tmax
  }
  if (missing(tlast)) {
    ret$tlast <- pk.calc.tlast(data$conc, data$time, check=FALSE)
  } else {
    ret$tlast <- tlast
  }
  # Data frame to use for computation of half-life
  if (allow.tmax.in.half.life) {
    # as.numeric is for units handling
    dfK <- data[as.numeric(data$time) >= as.numeric(ret$tmax), ]
  } else {
    # as.numeric is for units handling
    dfK <- data[as.numeric(data$time) > as.numeric(ret$tmax), ]
  }
  if (manually.selected.points) {
    if (nrow(data) > 0) {
      fit <- fit_half_life(data=data, tlast=ret$tlast, conc_units=conc_units)
      ret[,ret_replacements] <- fit[,ret_replacements]
    } else {
      warning("No data to manually fit for half-life (all concentrations may be 0 or excluded)")
      ret <-
        structure(
          ret,
          exclude="No data to manually fit for half-life (all concentrations may be 0 or excluded)"
        )
    }
  } else if (nrow(dfK) >= min.hl.points) {
    # If we have enough data to estimate a slope, then
    half_lives_for_selection <-
      data.frame(
        r.squared=-Inf,
        adj.r.squared=-Inf,
        clast.pred=NA_real_,
        lambda.z=-Inf,
        lambda.z.n.points=NA_integer_,
        lambda.z.time.first=dfK$time,
        log_conc=dfK$log_conc,
        span.ratio=NA_real_,
        half.life=NA_real_
      )
    half_lives_for_selection <-
      half_lives_for_selection[order(-half_lives_for_selection$lambda.z.time.first), ]
    for(i in min.hl.points:nrow(half_lives_for_selection)) {
      # Fit the terminal slopes until the adjusted r-squared value
      # is not improving (or it only gets worse by a small factor).
      fit <-
        fit_half_life(
          data=
            data.frame(
              # pass in the conc so that we can use its units, if applicable
              log_conc=half_lives_for_selection$log_conc[1:i],
              time=half_lives_for_selection$lambda.z.time.first[1:i]
            ),
          tlast=ret$tlast,
          conc_units=conc_units
        )
      half_lives_for_selection[i,names(fit)] <- fit
    }
    # Find the best model
    mask_best <-
      half_lives_for_selection$lambda.z > 0 &
      if (min.hl.points == 2 & nrow(half_lives_for_selection) == 2) {
        rlang::warn(
          message = "2 points used for half-life calculation",
          class = "pknca_halflife_2points"
        )
        TRUE
      } else {
        half_lives_for_selection$adj.r.squared >
          (max(half_lives_for_selection$adj.r.squared, na.rm=TRUE) - adj.r.squared.factor)
      }
    # Missing values are not the best
    mask_best[is.na(mask_best)] <- FALSE
    if (sum(mask_best) > 1) {
      # If more than one models qualify, choose the one with the
      # most points used.
      mask_best <-
        (mask_best &
           half_lives_for_selection$lambda.z.n.points == max(half_lives_for_selection$lambda.z.n.points[mask_best]))
    }
    # If the half-life fit, set all associated parameters
    if (any(mask_best)) {
      # Put in all the computed values
      ret[,ret_replacements] <- half_lives_for_selection[mask_best, ret_replacements]
    }
  } else {
    attr(ret, "exclude") <-
      sprintf(
        "Too few points for half-life calculation (min.hl.points=%g with only %g points)",
        min.hl.points, nrow(dfK)
      )
    rlang::warn(
      message = attr(ret, "exclude"),
      class = "pknca_halflife_too_few_points"
    )
  }
  # Drop the inputs of tmax and tlast, if given.
  if (!missing(tmax))
    ret$tmax <- NULL
  if (!missing(tlast))
    ret$tlast <- NULL
  ret
}

#' Perform the half-life fit given the data.  The function simply fits
#' the data without any validation.  No selection of points or any other
#' components are done.
#'
#' @param data The data to fit.  Must have two columns named "log_conc"
#'   and "time"
#' @param tlast The time of last observed concentration above the limit
#'   of quantification.
#' @param conc_units NULL or the units to set for concentration measures
#' @return A data.frame with one row and columns named "r.squared",
#'   "adj.r.squared", "PROB", "lambda.z", "clast.pred",
#'   "lambda.z.n.points", "half.life", "span.ratio"
#' @seealso [pk.calc.half.life()]
fit_half_life <- function(data, tlast, conc_units) {
  fit <- stats::.lm.fit(x=cbind(1, data$time), y=data$log_conc)
  # unit handling
  # if (inherits(tlast, "units")) {
  #   time_units <- units(tlast)
  # } else if (inherits(tlast, "mixed_units")) {
  #   time_units <- units(units::as_units(tlast))
  # } else {
    time_units <- NULL
  # }
  # if (!is.null(time_units)) {
  #   inverse_time_units <- time_units
  #   inverse_time_units$numerator <- time_units$denominator
  #   inverse_time_units$denominator <- time_units$numerator
  # } else {
    inverse_time_units <- NULL
  # }

  # as.numeric is so that it works for units objects
  r_squared <- 1 - as.numeric(sum(fit$residuals^2))/as.numeric(sum((data$log_conc - mean(data$log_conc))^2))
  clast_pred <- exp(sum(fit$coefficients*c(1, as.numeric(tlast))))
  # if (!is.null(conc_units)) {
  #   clast_pred <- units::set_units(clast_pred, conc_units, mode="standard")
  # }
  lambda_z <- -fit$coefficients[2]
  # if (!is.null(inverse_time_units)) {
  #   lambda_z <- units::set_units(lambda_z, inverse_time_units, mode="standard")
  # }
  ret <-
    data.frame(
      r.squared=r_squared,
      adj.r.squared=adj.r.squared(r_squared, nrow(data)),
      lambda.z=lambda_z,
      clast.pred=clast_pred,
      lambda.z.time.first=min(data$time, na.rm=TRUE),
      lambda.z.n.points=nrow(data)
    )
  ret$half.life <- log(2)/ret$lambda.z
  ret$span.ratio <- (max(data$time) - min(data$time))/ret$half.life
  ret
}

# Add the column to the interval specification
add.interval.col("half.life",
                 FUN="pk.calc.half.life",
                 values=c(FALSE, TRUE),
                 unit_type="time",
                 pretty_name="Half-life",
                 desc="The (terminal) half-life",
                 depends=c("tmax", "tlast"))
PKNCA.set.summary(
  name="half.life",
  description="arithmetic mean and standard deviation",
  point=business.mean,
  spread=business.sd
)
add.interval.col("r.squared",
                 FUN=NA,
                 values=c(FALSE, TRUE),
                 unit_type="unitless",
                 pretty_name="$r^2$",
                 desc="The r^2 value of the half-life calculation",
                 depends="half.life")
PKNCA.set.summary(
  name="r.squared",
  description="arithmetic mean and standard deviation",
  point=business.mean,
  spread=business.sd
)
add.interval.col("adj.r.squared",
                 FUN=NA,
                 values=c(FALSE, TRUE),
                 unit_type="unitless",
                 pretty_name="$r^2_{adj}$",
                 desc="The adjusted r^2 value of the half-life calculation",
                 depends="half.life")
PKNCA.set.summary(
  name="adj.r.squared",
  description="arithmetic mean and standard deviation",
  point=business.mean,
  spread=business.sd
)
add.interval.col("lambda.z",
                 FUN=NA,
                 values=c(FALSE, TRUE),
                 unit_type="inverse_time",
                 pretty_name="$\\lambda_z$",
                 desc="The elimination rate of the terminal half-life",
                 depends="half.life")
PKNCA.set.summary(
  name="lambda.z",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("lambda.z.time.first",
                 FUN=NA,
                 values=c(FALSE, TRUE),
                 unit_type="time",
                 pretty_name="First time for $\\lambda_z$",
                 desc="The first time point used for the calculation of half-life",
                 depends="half.life")
PKNCA.set.summary(
  name="lambda.z.time.first",
  description="median and range",
  point=business.median,
  spread=business.range
)
add.interval.col("lambda.z.n.points",
                 FUN=NA,
                 values=c(FALSE, TRUE),
                 unit_type="count",
                 pretty_name="Number of points used for lambda_z",
                 desc="The number of points used for the calculation of half-life",
                 depends="half.life")
PKNCA.set.summary(
  name="lambda.z.n.points",
  description="median and range",
  point=business.median,
  spread=business.range
)
add.interval.col("clast.pred",
                 FUN=NA,
                 values=c(FALSE, TRUE),
                 unit_type="conc",
                 pretty_name="Clast,pred",
                 desc="The concentration at Tlast as predicted by the half-life",
                 depends="half.life")
PKNCA.set.summary(
  name="clast.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("span.ratio",
                 FUN=NA,
                 values=c(FALSE, TRUE),
                 unit_type="fraction",
                 pretty_name="Span ratio",
                 desc="The ratio of the half-life to the duration used for half-life calculation",
                 depends="half.life")
PKNCA.set.summary(
  name="span.ratio",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

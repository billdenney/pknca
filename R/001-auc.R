#' Check the formatting of an AUC specification data frame.
#'
#' AUC specifications are data frames with four columns: \code{start},
#' \code{end}, \code{end.text}, \code{half.life}.
#'
#' The columns are defined as:
#' \describe{
#'   \item{\code{start}}{The starting time as a number}
#'   \item{\code{end}}{The ending time as either NA or a number (including Inf)}
#'   \item{\code{auc.type}}{What type of AUC should be computed: 'AUCinf',
#'     'AUClast', 'AUCall', or NA (for half-life only)}
#'   \item{\code{half.life}}{Compute the half-life for the interval (logical variable)}
#' }
#'
#' They are interpreted with the following rules:
#'
#' \itemize{
#'   \item The \code{start} time must always be given.
#'   \item The \code{start} must be before the \code{end} (if \code{end} is
#'         given)
#'   \item If the \code{end} time is given, then the other
#'         specifications (i.e. \code{auc.type} and
#'         \code{half.life} will only be done in that interval.
#'   \item If the \code{end} time is NA, then
#'   \itemize{
#'     \item No AUCs will be calculated.
#'     \item It is an error to set \code{auc.type} to anything other than NA.
#'     \item It is an error to set \code{half.life} to FALSE.
#'     \item \code{half.life} is computed from the \code{start} to the last
#'           available measurement.
#'   }
#' }
#'
#' @param x The data frame specifying what AUCs to calculate
#' @return x The potentially updated data frame with the AUC
#' specification.
#' @export
check.auc.specification <- function(x) {
  if (!is.data.frame(x)) {
    ## Just a warning and let as.data.frame make it an error if
    ## it can't be coerced.
    warning("AUC specification must be a data.frame")
    x <- as.data.frame(x, stringsAsFactors=FALSE)
  }
  if (nrow(x) == 0) {
    ## Return it as is-- nothing is requested
    return(x)
  }
  ## If end, last, all, or half.life is missing, add it as all NA.
  ## Inappropriate combinations will be checked later.
  added.cols <- list(end=NA, auc.type="AUClast", half.life=FALSE)
  for (new.col in names(added.cols))
    if (!(new.col %in% names(x))) {
      warning(sprintf("AUC specification column '%s' is missing.  Adding with all %s",
                      new.col, as.character(added.cols[[new.col]])))
      x[,new.col] <- added.cols[[new.col]]
    }
  required.cols <- c("start", "end", "auc.type", "half.life")
  if (!all(required.cols %in% names(x)))
    stop(paste("AUC specification must have columns for",
               paste(required.cols, collapse=", "),
               "Column(s) missing:",
               paste(setdiff(required.cols, names(x)))))
  ## If there are additional columns, remove them
  addl.cols <- setdiff(names(x), required.cols)
  if (length(addl.cols) > 0) {
    message(paste("AUC specification columns to specify the group(s):",
                  paste(addl.cols, collapse=", ")))
  }
  ## Ensure that all columns have the right edit(s)
  if (!is.logical(x$half.life)) {
    warning("AUC specification for 'half.life' must be a logical vector, attempting conversion")
    x$half.life <- check.conversion(x$half.life, as.logical)
  }
  if (!is.numeric(x$start)) {
    warning("AUC specification for 'start' time must be numeric, attempting conversion")
    x$start <- check.conversion(x$start, as.numeric)
  }
  if (!all(is.na(x$end) | is.numeric(x$end))) {
    warning("AUC specification for 'end' time must be numeric or NA, attempting conversion")
    x$end <- check.conversion(x$end, as.numeric)
  }
  ## Ensure that only valid combinations are provided
  if (any(is.na(x$start)))
    stop("AUC specification may not have NA for the starting time")
  if (!(is.character(x$auc.type) | is.factor(x$auc.type) | all(is.na(x$auc.type)))) {
    stop("auc.type must be either a character or factor (or NA)")
  }
  if (any(is.na(x$end) & !is.na(x$auc.type)))
    stop("AUC specification may not have NA for the end time and request an auc.type")
  if (any(is.na(x$end) & !x$half.life))
    stop("AUC specification may not have NA for the end time and not request half.life")
  if (any(!is.na(x$end) & x$end <= x$start))
    stop("AUC specification end must be after the start when end is given")
  x
}

#' A compute the Area Under the (Moment) Curve
#'
#' Compute the area under the curve (AUC) and the area under the
#' moment curve (AUMC) for pharmacokinetic (PK) data.  AUC and AUMC
#' are used for many purposes when analyzing PK in drug development.
#'
#' \code{pk.calc.auc.last} is simply a shortcut setting the
#' \code{interval} parameter to \code{c(0, "last")}.
#'
#' Extrapolation beyond Clast occurs using the half-life and
#' Clast,obs; Clast,pred is not yet supported.
#'
#' @param conc Concentration measured
#' @param time Time of concentration measurement (must be
#' monotonically increasing and the same length as the concentration
#' data)
#' @param interval Numeric vector of two numbers for the start and
#' end time of integration
#' @param auc.type The type of AUC to compute.  Choices are 'AUCinf',
#' 'AUClast', and 'AUCall'.
#' @param lambda.z The elimination rate (in units of inverse time) for
#' extrapolation
#' @param options List of changes to the default
#' \code{\link{PKNCA.options}} for calculations.
#' @param method The method for integration (either 'lin up/log down'
#' or 'linear')
#' @param conc.blq How to handle BLQ values in between the first and
#' last above LOQ concentrations. (See \code{\link{clean.conc.blq}}
#' for usage instructions.)
#' @param conc.na How to handle missing concentration values.  (See
#' \code{\link{clean.conc.na}} for usage instructions.)
#' @param check Run \code{\link{check.conc.time}},
#' \code{\link{clean.conc.blq}}, and \code{\link{clean.conc.na}}?
#' @param fun.linear The function to use for integration of the linear
#' part of the curve (not required for AUC or AUMC functions)
#' @param fun.log The function to use for integration of the
#' logarithmic part of the curve (if log integration is used; not
#' required for AUC or AUMC functions)
#' @param fun.inf The function to use for extrapolation from the final
#' measurement to infinite time (not required for AUC or AUMC
#' functions.
#' @param ... For functions other than \code{pk.calc.auxc}, these
#' values are passed to \code{pk.calc.auxc}
#' @return A numeric value for the AU(M)C
#' @aliases pk.calc.auc pk.calc.aumc pk.calc.auc.last
#' @seealso \code{\link{pk.calc.auc.all}},
#' \code{\link{pk.calc.auc.last}}, \code{\link{clean.conc.blq}}
#' @references
#'
#' Gabrielsson J, Weiner D.  "Section 2.8.1 Computation methods -
#' Linear trapezoidal rule."  Pharmacokinetic & Pharmacodynamic Data
#' Analysis: Concepts and Applications, 4th Edition.  Stockholm,
#' Sweden: Swedish Pharmaceutical Press, 2000.  162-4.
#'
#' Gabrielsson J, Weiner D.  "Section 2.8.3 Computation methods -
#' Log-linear trapezoidal rule."  Pharmacokinetic & Pharmacodynamic
#' Data Analysis: Concepts and Applications, 4th Edition.  Stockholm,
#' Sweden: Swedish Pharmaceutical Press, 2000.  164-7.
#' @examples
#' myconc <- c(0, 1, 2, 1, 0.5, 0.25, 0)
#' mytime <- c(0, 1, 2, 3, 4,   5,    6)
#' pk.calc.auc(myconc, mytime, interval=c(0, 6))
#' pk.calc.auc(myconc, mytime, interval=c(0, Inf))
#' @export
pk.calc.auxc <- function(conc, time, interval=c(0, Inf),
                         lambda.z=NA,
                         auc.type="AUClast",
                         options=list(),
                         method=PKNCA.choose.option("auc.method", options),
                         conc.blq=PKNCA.choose.option("conc.blq", options),
                         conc.na=PKNCA.choose.option("conc.na", options),
                         check=TRUE,
                         fun.linear, fun.log, fun.inf) {
  ## Check the inputs
  method <- PKNCA.options(auc.method=method, check=TRUE)
  if (check) {
    check.conc.time(conc, time)
    data <- clean.conc.blq(conc, time,
                           conc.blq=conc.blq,
                           conc.na=conc.na,
                           check=FALSE)
  } else {
    data <- data.frame(conc, time)
  }
  if (nrow(data) == 0) {
    ## All the data were missing
    return(NA)
  } else if (all(data$conc %in% c(0, NA))) {
    ## All the data were missing or 0
    return(0)
  }
  auc.type <- match.arg(auc.type, c("AUClast", "AUCinf", "AUCall"))
  if (interval[1] >= interval[2])
    stop("The AUC interval must be increasing")
  if (auc.type %in% "AUCinf" &
        is.finite(interval[2]))
    warning("Requesting AUCinf when the end of the interval is not Inf")
  ##############################
  ## Subset the data to the range of interest
  interval.start <- interval[1]
  interval.end <- interval[2]
  ## Find the first time point
  if (interval.start < min(data$time)) {
    warning(sprintf(
      "Requesting an AUC range starting (%g) before the first measurement (%g) is not allowed",
      interval.start, min(data$time)))
    return(NA)
  } else if (interval.start > max(data$time)) {
    ## Give this as a warning, but allow it to continue
    warning(sprintf("AUC start time (%g) is after the maximum observed time (%g)",
                    interval.start, max(data$time)))
  }
  ## Ensure that we have clean concentration and time data.  This
  ## means that we need to make sure that we have our starting point.
  ## Interpolation ensures that (and will give the same answer if it
  ## already exists in the right form).
  conc.start <-
    interp.extrap.conc(data$conc, data$time,
                       time.out=interval.start,
                       lambda.z=lambda.z,
                       interp.method=method,
                       extrap.method=auc.type,
                       check=FALSE)
  ## Add that concentration and time to the vectors removing the
  ## original version if it was there.
  data <- rbind(data.frame(conc=conc.start,
                           time=interval.start),
                data[!(data$time %in% interval.start),])
  ## * either have our ending point or the ending point is Inf
  if (is.finite(interval.end)) {
    conc.end <-
      interp.extrap.conc(data$conc, data$time,
                         interval.end,
                         interp.method=method,
                         extrap.method=auc.type,
                         check=FALSE)
    ## !mask.end because we may be replacing an entry that is a 0.
    data <- rbind(data[!(data$time %in% interval.end),],
                  data.frame(conc=conc.end,
                             time=interval.end))
  }
  ## Subset the conc and time data to the interval of interest
  data <- subset(data, (interval.start <= time &
                          time <= interval.end))
  ## Set the overall tlast
  tlast <- pk.calc.tlast(data$conc, data$time, check=FALSE)
  if (is.na(tlast)) {
    ## All concentrations are BLQ (note that this has to be checked
    ## after full subsetting and interpolation to ensure that it is
    ## still true)
    if (all(data$conc %in% 0)) {
      ret <- 0
    } else {
      stop("Unknown error with NA tlast but non-BLQ concentrations")
    }
  } else {
    ## ############################
    ## Compute the AUxC
    ## Compute it in linear space from the start to Tlast
    if (auc.type %in% "AUCall" &
          tlast != max(data$time)) {
      ## Include the first point after tlast if it exists and we are
      ## looking for AUCall
      idx.1 <- 1:sum(data$time <= tlast)
    } else {
      idx.1 <- 1:(sum(data$time <= tlast) - 1)
    }
    idx.2 <- idx.1 + 1
    ret <- fun.linear(data$conc[idx.1], data$conc[idx.2],
                      data$time[idx.1], data$time[idx.2])
    if (method %in% "lin up/log down") {
      ## Compute log down if required (and if the later point is not 0)
      mask.down <- (data$conc[idx.2] < data$conc[idx.1] &
                      data$conc[idx.2] != 0)
      idx.1 <- idx.1[mask.down]
      idx.2 <- idx.2[mask.down]
      ret[mask.down] <- fun.log(data$conc[idx.1], data$conc[idx.2],
                                data$time[idx.1], data$time[idx.2])
    } else if (!(method %in% "linear")) {
      stop("Invalid AUC integration method")
    }
    if (auc.type %in% "AUCinf") {
      ## Extrapolate to infinity using C.last.obs
      ## FIXME: should allow for C.last.pred, too.
      ret[length(ret)+1] <- fun.inf(data$conc, data$time, lambda.z)
    }
    ret <- sum(ret)
  }
  ret
}

fun.auc.linear <- function(conc.1, conc.2, time.1, time.2)
  (time.2-time.1)*(conc.2+conc.1)/2
fun.auc.log <- function(conc.1, conc.2, time.1, time.2)
  (time.2-time.1)*(conc.2-conc.1)/log(conc.2/conc.1)
fun.auc.inf <- function(conc, time, lambda.z)
  pk.calc.clast.obs(conc, time, check=FALSE)/lambda.z

#' @describeIn pk.calc.auxc Compute the area under the curve
#' @export
pk.calc.auc <- function(...)
  pk.calc.auxc(...,
               fun.linear=fun.auc.linear,
               fun.log=fun.auc.log,
               fun.inf=fun.auc.inf)

## Note that lambda.z is set to NA for both auc.last and auc.all
## because all interpolation should happen within given points.
## lambda.z should not be used, and if it is used, that should be
## caught as an error.
#' @describeIn pk.calc.auxc Compute the AUClast.
#' @export
pk.calc.auc.last <- function(...)
  pk.calc.auc(...,
              auc.type="AUClast",
              lambda.z=NA)

#' @describeIn pk.calc.auxc Compute the AUCall.
#' @export
pk.calc.auc.all <- function(...)
  pk.calc.auc(...,
              auc.type="AUCall",
              lambda.z=NA)

#' @describeIn pk.calc.auxc Compute the area under the moment curve
#' @export
pk.calc.aumc <- function(...)
  pk.calc.auxc(...,
    fun.linear=function(conc.1, conc.2, time.1, time.2) {
      (time.2-time.1)*(conc.2*time.2+conc.1*time.1)/2
    },
    fun.log=function(conc.1, conc.2, time.1, time.2) {
      ((time.2-time.1)*(conc.2*time.2-conc.1*time.1)/log(conc.2/conc.1)-
       (time.2-time.1)^2*(conc.2-conc.1)/(log(conc.2/conc.1)^2))
    },
    fun.inf=function(conc.last, time.last, lambda.z) {
      (conc.last*time.last/lambda.z) + conc.last/(lambda.z^2)
    })

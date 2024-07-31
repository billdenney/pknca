#' Generate a sparse_pk object
#'
#' @inheritParams assert_conc_time
#' @param subject Subject identifiers (may be any class; may not be null)
#' @returns A sparse_pk object which is a list of lists.  The inner lists have
#'   elements named: "time", The time of measurement; "conc", The concentration
#'   measured; "subject", The subject identifiers.  The object will usually be
#'   modified by future functions to add more named elements to the inner list.
#' @family Sparse Methods
#' @export
as_sparse_pk <- function(conc, time, subject) {
  if (is.data.frame(conc) && missing(time) && missing(subject)) {
    time <- conc$time
    subject <- conc$subject
    conc <- conc$conc
  }
  assert_conc_time(conc = conc, time = time, any_missing_conc = FALSE, sorted_time = FALSE)
  checkmate::check_vector(subject, any.missing=FALSE, len=length(conc), null.ok=FALSE)
  unique_times <- sort(unique(time))
  ret <- list()
  for (current_time in unique_times) {
    current_mask <- time %in% current_time
    ret <-
      append(
        ret,
        list(list(
          time=current_time,
          conc=conc[current_mask],
          subject=subject[current_mask]
        ))
      )
  }
  class(ret) <- "sparse_pk"
  ret
}

#' Set or get a sparse_pk object attribute
#'
#' @param sparse_pk A sparse_pk object from [as_sparse_pk()]
#' @param ... Either a character string (to get that value) or a named vector
#'   the same length as `sparse_pk` to set the value.
#' @returns Either the attribute value or an updated `sparse_pk` object
#' @keywords Internal
sparse_pk_attribute <- function(sparse_pk, ...) {
  args <- list(...)
  stopifnot(length(args) == 1)
  if (is.null(names(args))) {
    vapply(X=sparse_pk, FUN="[[", args[[1]], FUN.VALUE = 1)
  } else {
    stopifnot(length(args[[1]]) == length(sparse_pk))
    for (idx in seq_along(sparse_pk)) {
      sparse_pk[[idx]][names(args)[1]] <- args[[1]][idx]
    }
    sparse_pk
  }
}

#' Calculate the weight for sparse AUC calculation with the linear-trapezoidal
#' rule
#'
#' The weight is used as the \eqn{w_i}{w_i} parameter in [pk.calc.sparse_auc()]
#'
#' \deqn{w_i = \frac{\delta_{time,i-1,i} + \delta_{time,i,i+1}}{2}}{w_i = (d_time[i-1,i] + d_time[i,i+1])/2}
#' \deqn{\delta_{time,i,i+1} = t_{i+1} - t_i}{d_time = t_[i+1] - t_i, and zero if i < 1 or i > K}
#'
#' Where:
#'
#' \describe{
#'   \item{\eqn{w_i}{w_i}}{is the weight at time i}
#'   \item{\eqn{\delta_{time,i-1,i}}{d_time[i-1,i]} and \eqn{\delta_{time,i,i+1}}{d_time[i,i+1]}}{are the changes between time i-1 and i or i and i+1 (zero outside of the time range)}
#'   \item{\eqn{t_i}{t_i}}{is the time at time i}
#' }
#'
#' @inheritParams sparse_pk_attribute
#' @returns A numeric vector of weights for sparse AUC calculations the same
#'   length as `sparse_pk`
#' @family Sparse Methods
#' @export
sparse_auc_weight_linear <- function(sparse_pk) {
  times <- vapply(X=sparse_pk, FUN="[[", "time", FUN.VALUE = 1)
  half_diff_times <- diff(times)/2
  weights <- c(0, half_diff_times) + c(half_diff_times, 0)
  sparse_pk_attribute(sparse_pk=sparse_pk, weight=weights)
}

#' Calculate the mean concentration at all time points for use in sparse NCA
#' calculations
#'
#' Choices for the method of calculation (the argument `sparse_mean_method`)
#' are:
#'
#' \describe{
#'   \item{"arithmetic mean"}{Arithmetic mean (ignoring number of BLQ samples)}
#'   \item{"arithmetic mean, <=50% BLQ"}{If >= 50% of the measurements are BLQ, zero.  Otherwise, the arithmetic mean of all samples (including the BLQ as zero).}
#' }
#'
#' @inheritParams sparse_pk_attribute
#' @param sparse_mean_method The method used to calculate the sparse mean (see
#'   details)
#' @returns A vector the same length as `sparse_pk` with the mean concentration
#'   at each of those times.
#' @family Sparse Methods
#' @export
sparse_mean <- function(sparse_pk, sparse_mean_method=c("arithmetic mean, <=50% BLQ", "arithmetic mean")) {
  sparse_mean_method <- match.arg(sparse_mean_method)
  ret <-
    vapply(
      X = sparse_pk,
      FUN = function(current_time) mean(current_time$conc),
      FUN.VALUE = 1
    )
  if (sparse_mean_method == "arithmetic mean, <=50% BLQ") {
    numerator <-
      vapply(
        X=sparse_pk,
        FUN=function(current_time) sum(current_time$conc == 0),
        FUN.VALUE = 1
      )
    denominator <-
      vapply(
        X = sparse_pk,
        FUN = function(current_time) length(current_time$conc),
        FUN.VALUE = 1
      )
    frac_blq <- numerator/denominator
    ret[frac_blq > 0.5] <- 0
  } else if (sparse_mean_method == "arithmetic mean") {
    # do nothing
  } else {
    stop("Invalid sparse_mean_method: ", sparse_mean_method) # nocov
  }
  sparse_pk <- sparse_pk_attribute(sparse_pk, mean=ret)
  sparse_pk <- sparse_pk_attribute(sparse_pk, mean_method=rep(sparse_mean_method, length(ret)))
  sparse_pk
}

#' Calculate the variance for the AUC of sparsely sampled PK
#'
#' Equation 7.vii in Nedelman and Jia, 1998 is used for this calculation:
#'
#' \deqn{var\left(\hat{AUC}\right) = \sum\limits_{i=0}^m\left(\frac{w_i^2 s_i^2}{r_i}\right) + 2\sum\limits_{i<j}\left(\frac{w_i w_j r_{ij} s_{ij}}{r_i r_j}\right)}{var(AUC) = sum_(i=0)^(m) ((w_i^2 * s_i^2)/(r_i) + + 2*sum_(i<j)((w_i * w_j * r_ij * s_ij)/(r_i * r_j))}
#'
#' The degrees of freedom are calculated as described in equation 6 of the same
#' paper.
#'
#' @inheritParams sparse_pk_attribute
#' @references
#' Nedelman JR, Jia X. An extension of Satterthwaite’s approximation applied to
#' pharmacokinetics. Journal of Biopharmaceutical Statistics. 1998;8(2):317-328.
#' doi:10.1080/10543409808835241
#' @export
var_sparse_auc <- function(sparse_pk) {
  covariance <- cov_holder(sparse_pk)
  var_auc <- 0
  weights <- sparse_pk_attribute(sparse_pk, "weight")
  # number of subjects at a given time point
  n <- rep(0, length(sparse_pk))
  df <- 0
  for (idx1 in seq_along(sparse_pk)) {
    n_idx1 <- length(unique(sparse_pk[[idx1]]$subject))
    n[idx1] <- n_idx1
    var_auc <-
      var_auc +
      weights[idx1]^2*covariance[idx1, idx1]/n_idx1
    for (idx2 in seq_len(idx1 - 1)) {
      n_idx2 <- length(unique(sparse_pk[[idx2]]$subject))
      n_both <- length(unique(intersect(sparse_pk[[idx1]]$subject, sparse_pk[[idx2]]$subject)))
      var_auc <-
        var_auc +
        2*weights[idx1]*weights[idx2]*n_both*covariance[idx1, idx2]/(n_idx1*n_idx2)
    }
  }
  # df based on equation 6 of Nedelman and Jia 1998
  # df_e <- sum(diag(covariance))
  # df_v <- 2*sum(diag(covariance %*% covariance))
  # df <- 2*df_e^2/df_v
  # df based on equation 6a of Nedelman et al 1995
  df <-
    sum(weights^2 * diag(covariance)/n)^2 /
    sum(weights^4 * diag(covariance)^2/(n^2*(n-1)))
  if (sum(covariance[lower.tri(covariance)] != 0) > 0) {
    rlang::warn(
      message = "Cannot yet calculate sparse degrees of freedom for multiple samples per subject",
      class = "pknca_sparse_df_multi"
    )
    df <- NA_real_
  }
  attr(var_auc, "df") <- df
  var_auc
}

#' Calculate the covariance for two time points with sparse sampling
#'
#' The calculation follows equation A3 in Holder 2001 (see references below):
#'
#' \deqn{\hat{\sigma}_{ij} = \sum\limits_{k=1}^{r_{ij}}{\frac{\left(x_{ik} - \bar{x}_i\right)\left(x_{jk} - \bar{x}_j\right)}{\left(r_{ij} - 1\right) + \left(1 - \frac{r_{ij}}{r_i}\right)\left(1 - \frac{r_{ij}}{r_j}\right)}}}{sigma_ij = sum_(k=1)^(r_ij)((x_ik-xbar_i)(x_jk-xbar_j)/((r_ij-1)+(1-r_ij/r_i)*(1-r_ij/r_j)))}
#'
#' If \eqn{r_{ij} = 0}{r_ij = 0}, then \eqn{\hat{\sigma}_{ij}}{sigma_ij} is
#' defined as zero (rather than dividing by zero).
#'
#' Where:
#' \describe{
#'   \item{\eqn{\hat{\sigma}_{ij}}{sigma_ij}}{The covariance of times i and j}
#'   \item{\eqn{r_i}{r_i} and \eqn{r_j}{r_j}}{The number of subjects (usually animals) at times i and j, respectively}
#'   \item{\eqn{r_{ij}{r_ij}}}{The number of subjects (usually animals) at both times i and j}
#'   \item{\eqn{x_{ik}}{x_ik} and \eqn{x_{jk}}{x_jk}}{The concentration measured for animal k at times i and j, respectively}
#'   \item{\eqn{\bar{x}_i}{xbar_i} and \eqn{\bar{x}_j}{xbar_j}}{The mean of the concentrations at times i and j, respectively}
#' }
#'
#' The Cauchy-Schwartz inequality is enforced for covariances to keep
#' correlation coefficients between -1 and 1, inclusive, as described in
#' equations 8 and 9 of Nedelman and Jia 1998.
#'
#' @inheritParams sparse_pk_attribute
#' @returns A matrix with one row and one column for each element of
#'   `sparse_pk_attribute`.  The covariances are on the off diagonals, and for
#'   simplicity of use, it also calculates the variance on the diagonal
#'   elements.
#' @keywords Internal
#' @references
#' Holder DJ. Comments on Nedelman and Jia’s Extension of Satterthwaite’s
#' Approximation Applied to Pharmacokinetics. Journal of Biopharmaceutical
#' Statistics. 2001;11(1-2):75-79. doi:10.1081/BIP-100104199
#'
#' Nedelman JR, Jia X. An extension of Satterthwaite’s approximation applied to
#' pharmacokinetics. Journal of Biopharmaceutical Statistics. 1998;8(2):317-328.
#' doi:10.1080/10543409808835241
#' @export
cov_holder <- function(sparse_pk) {
  ret <-
    matrix(
      data=0,
      nrow=length(sparse_pk),
      ncol=length(sparse_pk)
    )

  time_means <- sparse_pk_attribute(sparse_pk, "mean")

  for (idx1 in seq_along(sparse_pk)) {
    # Variance on the diagonal
    ret[idx1, idx1] <- stats::var(sparse_pk[[idx1]]$conc)
    for (idx2 in seq_len(idx1 - 1)) {
      subject_idx1 <- sparse_pk[[idx1]]$subject
      subject_idx2 <- sparse_pk[[idx2]]$subject
      subject_both <- intersect(subject_idx1, subject_idx2)
      if (length(subject_both) > 1) {
        # Holder covariance on the off-diagonals when there is more than one
        # subject in both times
        cov_ij <- 0
        for (current_subject in subject_both) {
          cov_ij <-
            cov_ij +
            (sparse_pk[[idx1]]$conc[sparse_pk[[idx1]]$subject %in% current_subject] - time_means[[idx1]]) *
            (sparse_pk[[idx2]]$conc[sparse_pk[[idx2]]$subject %in% current_subject] - time_means[[idx2]])
        }
        # Apply the common denominator
        cov_ij <-
          cov_ij /
          (
            (length(subject_both) - 1) + (1 - length(subject_both)/length(subject_idx1))*(1 - length(subject_both)/length(subject_idx2))
          )
        # Enforce the Cauchy-Schwartz inequality
        cov_cs <- sqrt(ret[idx1, idx1] * ret[idx2, idx2])
        if (abs(cov_ij) > cov_cs) {
          cov_ij <- sign(cov_ij)*cov_cs
        }
        # The matrix is symmetric
        ret[idx1, idx2] <- ret[idx2, idx1] <- cov_ij
      }
    }
  }
  ret
}

#' Extract the mean concentration-time profile as a data.frame
#'
#' @inheritParams sparse_pk_attribute
#' @return A data.frame with names of "conc" and "time"
#' @keywords Internal
sparse_to_dense_pk <- function(sparse_pk) {
  data.frame(
    conc=sparse_pk_attribute(sparse_pk, "mean"),
    time=sparse_pk_attribute(sparse_pk, "time")
  )
}

#' Calculate AUC and related parameters using sparse NCA methods
#'
#' The AUC is calculated as:
#'
#' \deqn{AUC=\sum\limits_{i} w_i \bar{C}_i}{AUC = sum(w_i * Cbar_i)}
#'
#' Where:
#'
#' \describe{
#'   \item{\eqn{AUC}{AUC}}{is the estimated area under the concentration-time curve}
#'   \item{\eqn{w_i}{w_i}}{is the weight applied to the concentration at time i (related to the time which it affects, see [sparse_auc_weight_linear()])}
#'   \item{\eqn{\bar{C}_i}{Cbar_i}}{is the average concentration at time i}
#' }
#' @inheritParams pk.calc.auc
#' @inheritParams as_sparse_pk
#' @family Sparse Methods
#' @export
pk.calc.sparse_auc <- function(conc, time, subject,
                               method=NULL,
                               auc.type="AUClast",
                               ...,
                               options=list()) {
  sparse_pk <- as_sparse_pk(conc=conc, time=time, subject=subject)
  sparse_pk_wt <- sparse_auc_weight_linear(sparse_pk)
  sparse_pk_mean <- sparse_mean(sparse_pk=sparse_pk_wt, sparse_mean_method="arithmetic mean, <=50% BLQ")
  auc <-
    pk.calc.auc(
      conc=sparse_pk_attribute(sparse_pk_mean, "mean"),
      time=sparse_pk_attribute(sparse_pk_mean, "time"),
      auc.type=auc.type,
      method="linear"
    )
  var_auc <- var_sparse_auc(sparse_pk_mean)
  data.frame(
    sparse_auc=auc,
    # as.numeric() drops the "df" attribute
    sparse_auc_se=sqrt(as.numeric(var_auc)),
    sparse_auc_df=attr(var_auc, "df")
  )
}

#' @describeIn pk.calc.sparse_auc Compute the AUClast for sparse PK
#' @export
pk.calc.sparse_auclast <- function(conc, time, subject, ..., options=list()) {
  if ("auc.type" %in% names(list(...))) {
    rlang::abort(
      message = "auc.type cannot be changed when calling pk.calc.sparse_auclast, please use pk.calc.sparse_auc",
      class = "pknca_sparse_auclast_change_auclast"
    )
  }
  ret <-
    pk.calc.sparse_auc(
      conc=conc, time=time, subject=subject, ...,
      options=options,
      auc.type="AUClast",
      lambda.z=NA
    )
  names(ret)[names(ret) == "sparse_auc"] <- "sparse_auclast"
  ret
}

add.interval.col(
  "sparse_auclast",
  sparse=TRUE,
  FUN="pk.calc.sparse_auclast",
  values=c(FALSE, TRUE),
  unit_type="auc",
  pretty_name="Sparse AUClast",
  desc="For sparse PK sampling, the area under the concentration time curve from the beginning of the interval to the last concentration above the limit of quantification"
)
PKNCA.set.summary(
  name="sparse_auclast",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

add.interval.col(
  "sparse_auc_se",
  FUN=NA,
  values=c(FALSE, TRUE),
  unit_type="auc",
  pretty_name="Sparse AUClast standard error",
  desc="For sparse PK sampling, the standard error of the area under the concentration time curve from the beginning of the interval to the last concentration above the limit of quantification",
  depends="sparse_auclast"
)
PKNCA.set.summary(
  name="sparse_auc_se",
  description="arithmetic mean and standard deviation",
  point=business.mean,
  spread=business.sd
)

add.interval.col(
  "sparse_auc_df",
  FUN=NA,
  values=c(FALSE, TRUE),
  unit_type="count",
  pretty_name="Sparse AUClast degrees of freedom",
  desc="For sparse PK sampling, the standard error degrees of freedom of the area under the concentration time curve from the beginning of the interval to the last concentration above the limit of quantification",
  depends="sparse_auclast"
)
PKNCA.set.summary(
  name="sparse_auc_df",
  description="arithmetic mean and standard deviation",
  point=business.mean,
  spread=business.sd
)

#' Is a PKNCA object used for sparse PK?
#'
#' @param object The object to see if it includes sparse PK
#' @returns `TRUE` if sparse and `FALSE` if dense (not sparse)
#' @export
is_sparse_pk <- function(object) {
  UseMethod("is_sparse_pk")
}

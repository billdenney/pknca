#' Calculate Satterwaithe correlation for sparse AUC using Bailer's method to
#' provide AUC, SE, and degrees of freedom.
#' 
#' @param m number of times
#' @param a number of animals
#' @param D design matrix
#' @param w vector of trapezoidal weights
#' @param y data structured as D
#' @return auc = estimated AUC
#' se.auc = estimated standard error of estimated AUC
#' df = estimated degrees of freedom
#' width = half-width of the 95% confidence interval
#' lower = lower bound of the 95% confidence interval
#' upper = upper bound of the 95% confidence interval
#' @references The function here is copied from the appendix of
#' Nedelman JR, Jia X. An extension of Satterthwaite’s approximation applied to
#' pharmacokinetics. Journal of Biopharmaceutical Statistics. 1998;8(2):317-328.
#' doi:10.1080/10543409808835241
corsatt <- function(m, a, D, w, y) {
  # Create identity matrices and matrices of ones
  Ia <- diag(rep(1, a))
  Im <- diag(rep(1, m))
  Ima <- diag(rep(1, m*a))
  Ja <- matrix(1, ncol=a, nrow=a)
  Jma <- matrix(1, nrow=m, ncol=a)
  
  R <- D %*% t(D)
  
  Delta <- diag(c(D), nrow=m*a, ncol=m*a)
  DeltaStar <- diag(c(D)/kronecker(rep(1, a), diag(R)))
  
  A <- diag(w/diag(R)) %*% (R/(R - 1)) %*% diag(w/diag(R))
  
  P <- kronecker(Ja, Im) %*% DeltaStar
  
  M <- t(Ima - P) %*% Delta %*% kronecker(Ia, A) %*% Delta %*% (Ima - P)
  
  # Convert zeros where animals are not observed to NA
  y[D == 0] <- NA
  
  # Compute vector of means
  y_mean <- matrix(rowMeans(x=y, na.rm=TRUE), nrow=m, ncol=1)
  
  # Compute AUC from the data using the linear trapezoidal rule
  auc <- t(w) %*% y_mean
  
  # Estimate covariance matrix from the data
  sigest <- matrix(NA_real_, nrow=m, ncol=m)
  y_dif <- y - kronecker(matrix(1, nrow=1, ncol=a), y_mean)
  # first compute variances (diagonals of cov matrix)
  for (j1 in seq_len(nrow(y))) {
    r <- R[j1, j1]
    sigest[j1, j1] <- mean(y_dif[j1,] * y_dif[j1,], na.rm=TRUE) * r/(r - 1)
  }
  # now compute above the diagonal
  for (j1 in 1:(nrow(y) - 1)) {
    for (j2 in (j1 + 1):nrow(y)) {
      r <- R[j1, j2]
      cov <- mean(y_dif[j1,] * y_dif[j2,], na.rm=TRUE) * r/(r - 1)
      # prepare to enforce the Cauchy-Schwartz inequality
      cs_cov <- sign(cov)*sqrt(sigest[j1, j1]*sigest[j2, j2])
      if (is.na(cov)) {
        sigest[j1, j2] <- 0
      } else if (abs(cov) <= abs(cs_cov)) {
        sigest[j1, j2] <- cov
      } else {
        # Must enforce the Cauchy-Schwartz inequality
        sigest[j1, j2] <- cs_cov
      }
      # Now fill in below the diagonal
      sigest[j2, j1] <- sigest[j1, j2]
    }
  }
  
  # Compute df and confidence interval information
  EV <- M %*% kronecker(Ia, sigest)
  E <- sum(diag(EV))
  V <- 2*sum(diag(EV %*% EV))
  df <- 2*E^2/V
  
  data.frame(
    auc=auc,
    se_auc=sqrt(E),
    df=df
  )
}

# Example 1a from Yeh, C. Estimation and significant tests of area under the
# curve derived from incomplete sampling. ASA Proceedings of the
# Biopharmaceutical Section, 1990, pp. 74-81.

# Number of times is
mm <- 9

# Number of animals is
aa <- 9

# D matrix
DD <- matrix(rep(c(1, 0, 0,
                   0, 1, 0,
                   1, 0, 0,
                   0, 1, 0,
                   1, 0, 0,
                   0, 1, 0,
                   0, 0, 1,
                   0, 0, 1,
                   0, 0, 1), each=3),
             nrow=mm, ncol=aa, byrow=TRUE)

time <- c(0, 0.5, 1, 2, 4, 6, 8, 12, 24)
# AUC weights
ww_prep <- time[-1] - time[-length(time)]
ww <- c(0, ww_prep/2) + c(ww_prep/2, 0)
conc_raw <-
  c(
    0, 0, 0,
    4, 1.3, 3.2,
    4.69, 2.07, 6.45,
    6.68, 3.83, 6.08,
    4.69, 4.06, 6.45, # Accurate relative to the source paper; duplicated data relative to two rows ago suggests possible transcription error
    8.13, 9.54, 6.29,
    9.36, 13, 5.48,
    5.18, 5.18, 2.79,
    1.06, 2.15, 0.827
  )
conc <- t(DD)
conc[t(DD) == 1] <- conc_raw
conc <- t(conc)

# According to source paper, should return
# data.frame(auc=110.1015, se_auc=14.78964, df=2.144318)
corsatt(mm, aa, DD, ww, conc)

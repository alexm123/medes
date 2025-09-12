#' Upsilon from covariance matrix only
#'
#' Compute unadjusted and adjusted upsilon mediation measures using only a
#' covariance matrix (of x, m, y) and the sample size n. Optionally, obtain
#' bootstrap confidence intervals via parametric bootstrap (MVN simulation).
#'
#' @param S A 3x3 covariance matrix whose row/col names include x, m, y.
#'          Order can be arbitrary but names must be present.
#' @param x,m,y Character names of the variables in S.
#' @param n Integer, sample size used to compute S (number of observations).
#' @param R Number of bootstrap draws (default 1000).
#' @param ci_level Confidence level for intervals (default 0.95).
#' @param do_bootstrap Logical; if TRUE, run parametric bootstrap (default TRUE).
#' @param stat One of "adjusted", "unadjusted", or "both" (default "both").
#' @return A named list with point estimates and (optionally) CIs.
#' @examples
#' set.seed(42)
#' # Fake covariance with var(x)=1, var(m)=1, var(y)=1.5, some covariances
#' S <- matrix(c(
#'   1.00, 0.30, 0.25,
#'   0.30, 1.00, 0.40,
#'   0.25, 0.40, 1.50
#' ), 3, 3, byrow = TRUE)
#' dimnames(S) <- list(c("x","m","y"), c("x","m","y"))
#' upsilons_cov(S, x="x", m="m", y="y", n=500, R=500, do_bootstrap=FALSE)
upsilons_cov <- function(S, x, m, y,
                         n,
                         R = 1000,
                         ci_level = 0.95,
                         do_bootstrap = TRUE,
                         stat = c("adjusted", "unadjusted", "both")) {

  stat <- match.arg(stat)

  # -------- Input checks --------
  if (!is.matrix(S) || nrow(S) != 3 || ncol(S) != 3)
    stop("`S` must be a 3x3 covariance matrix for variables x, m, y.")

  if (is.null(rownames(S)) || is.null(colnames(S)))
    stop("`S` must have row and column names that include x, m, y.")

  needed <- c(x, m, y)
  if (!all(needed %in% rownames(S)) || !all(needed %in% colnames(S)))
    stop("Row/col names of `S` must include: ", paste(needed, collapse=", "))

  if (!isTRUE(all.equal(rownames(S), colnames(S))))
    stop("Row and column names of `S` must match in the same order.")

  if (!is.numeric(n) || length(n) != 1 || n < 3)
    stop("`n` must be a single integer >= 3 (sample size).")

  if (!is.numeric(R) || R < 1) stop("`R` must be a positive integer.")
  if (!is.numeric(ci_level) || ci_level <= 0 || ci_level >= 1)
    stop("`ci_level` must be in (0, 1).")

  # Reorder S so that its rows/cols are [x, m, y]
  S <- S[needed, needed]

  # Basic symbols
  Sxx <- S[x, x]
  Smm <- S[m, m]
  Syy <- S[y, y]
  Sxm <- S[x, m]
  Sxy <- S[x, y]
  Smy <- S[m, y]

  # -------- Coefficients via covariance algebra (OLS) --------
  # a-path: m ~ x
  a_hat <- Sxm / Sxx

  # b-path (and c'): y ~ x + m
  # Z = [x, m]
  S_ZZ <- matrix(c(Sxx, Sxm, Sxm, Smm), 2, 2)
  S_Zy <- c(Sxy, Smy)
  S_ZZ_inv <- solve(S_ZZ)
  beta <- S_ZZ_inv %*% S_Zy
  # names: c("b_x","b_m"), we want b_m
  b_hat <- as.numeric(beta[2])

  # -------- Variance of a_hat and b_hat from S, n --------
  # For simple regression (m ~ x):
  # R2_mx = Sxm^2 / (Sxx * Smm)
  R2_mx <- (Sxm^2) / (Sxx * Smm)
  sigma2_m_given_x <- Smm * (1 - R2_mx)
  var_a <- sigma2_m_given_x / (n * Sxx)

  # For multiple regression (y ~ x + m):
  # R2 = S_Zy' S_ZZ^{-1} S_Zy / Syy
  R2_y_Z <- as.numeric(t(S_Zy) %*% S_ZZ_inv %*% S_Zy) / Syy
  sigma2_y_given_Z <- Syy * (1 - R2_y_Z)
  # Var(beta) = (sigma^2 / n) * S_ZZ^{-1}
  var_beta <- (sigma2_y_given_Z / n) * S_ZZ_inv
  var_b <- as.numeric(var_beta[2, 2])

  # -------- Upsilon (OLS-based) --------
  upsilon_ols <- (a_hat^2 * b_hat^2) * (Sxx / Syy)
  adj_upsilon_ols <- ((a_hat^2 - var_a) * (b_hat^2 - var_b)) * (Sxx / Syy)

  out <- list()
  if (stat %in% c("unadjusted", "both")) out$upsilon_ols <- as.numeric(upsilon_ols)
  if (stat %in% c("adjusted",   "both")) out$`adj.upsilon_ols` <- as.numeric(adj_upsilon_ols)

  # -------- Parametric bootstrap (optional) --------
  if (isTRUE(do_bootstrap)) {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is required for parametric bootstrap. Please install it or set do_bootstrap=FALSE.")
    }

    alpha <- 1 - ci_level

    boot_fun <- function() {
      # Simulate centered data ~ MVN(0, S) with sample size n
      Xsim <- MASS::mvrnorm(n = n, mu = c(0, 0, 0), Sigma = S)
      colnames(Xsim) <- c(x, m, y)  # already ordered

      # Sample covariance from the simulated data
      Sc <- stats::cov(Xsim)

      # Recompute a, b, upsilons using the *sample* covariance
      Sxx_b <- Sc[x, x]; Smm_b <- Sc[m, m]; Syy_b <- Sc[y, y]
      Sxm_b <- Sc[x, m]; Sxy_b <- Sc[x, y]; Smy_b <- Sc[m, y]

      a_b <- Sxm_b / Sxx_b
      S_ZZ_b <- matrix(c(Sxx_b, Sxm_b, Sxm_b, Smm_b), 2, 2)
      S_Zy_b <- c(Sxy_b, Smy_b)
      S_ZZ_inv_b <- tryCatch(solve(S_ZZ_b), error = function(e) return(NA))
      if (is.na(S_ZZ_inv_b)[1]) return(c(NA_real_, NA_real_))

      beta_b <- S_ZZ_inv_b %*% S_Zy_b
      b_b <- as.numeric(beta_b[2])

      # Variances for adjustment at this bootstrap draw
      R2_mx_b <- (Sxm_b^2) / (Sxx_b * Smm_b)
      sigma2_m_given_x_b <- Smm_b * (1 - R2_mx_b)
      var_a_b <- sigma2_m_given_x_b / (n * Sxx_b)

      R2_y_Z_b <- as.numeric(t(S_Zy_b) %*% S_ZZ_inv_b %*% S_Zy_b) / Syy_b
      sigma2_y_given_Z_b <- Syy_b * (1 - R2_y_Z_b)
      var_beta_b <- (sigma2_y_given_Z_b / n) * S_ZZ_inv_b
      var_b_b <- as.numeric(var_beta_b[2, 2])

      ups_b  <- (a_b^2 * b_b^2) * (Sxx_b / Syy_b)
      adj_b  <- ((a_b^2 - var_a_b) * (b_b^2 - var_b_b)) * (Sxx_b / Syy_b)
      c(ups_b, adj_b)
    }

    boot_mat <- replicate(R, boot_fun())
    boot_mat <- t(boot_mat)

    if (stat %in% c("unadjusted", "both")) {
      out$ci.upsilon_ols <- stats::quantile(boot_mat[, 1],
                                            probs = c(alpha/2, 1 - alpha/2),
                                            na.rm = TRUE)
    }
    if (stat %in% c("adjusted", "both")) {
      out$`ci.adj.upsilon_ols` <- stats::quantile(boot_mat[, 2],
                                                  probs = c(alpha/2, 1 - alpha/2),
                                                  na.rm = TRUE)
    }
  }

  return(out)
}

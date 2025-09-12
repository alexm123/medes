implied_cov_mediation <- function(tau_prime, alpha, beta, check_psd = TRUE) {
  stopifnot(length(tau_prime)==1, length(alpha)==1, length(beta)==1)
  if (abs(alpha) > 1) stop("|alpha| must be <= 1 (since Var(X)=Var(M)=1).")

  ey_var <- 1 - (tau_prime^2 + beta^2 + 2*tau_prime*beta*alpha)
  if (ey_var < -1e-10) {
    stop("Invalid paths: implied Var(e_y) < 0. Ensure tau'^2 + beta^2 + 2*tau'*beta*alpha <= 1.")
  }

  Sigma <- matrix(
    c(1,                     alpha,                  tau_prime + beta*alpha,
      alpha,                 1,                      beta + alpha*tau_prime,
      tau_prime + beta*alpha, beta + alpha*tau_prime, 1),
    nrow = 3, byrow = TRUE
  )
  colnames(Sigma) <- rownames(Sigma) <- c("x","m","y")

  psd_ok <- TRUE
  if (check_psd) {
    vals <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
    psd_ok <- all(vals >= -1e-8)
    if (!psd_ok) stop("Implied covariance matrix is not positive semidefinite; check inputs.")
  }

  list(Sigma = Sigma, ey_var = max(ey_var, 0), psd_ok = psd_ok)
}

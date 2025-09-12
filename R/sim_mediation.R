sim_mediation <- function(tau_prime, alpha, beta, n = NULL, empirical = FALSE) {
  stopifnot(is.numeric(tau_prime), length(tau_prime) == 1,
            is.numeric(alpha),     length(alpha)     == 1,
            is.numeric(beta),      length(beta)      == 1)

  # Residual variance for Y under standardized SEM
  ey_var <- 1 - (tau_prime^2 + beta^2 + 2 * tau_prime * beta * alpha)
  if (ey_var < -1e-10) {
    stop("Invalid path values: implied residual variance for Y is negative.\n",
         "Ensure tau_prime^2 + beta^2 + 2*tau_prime*beta*alpha <= 1.")
  }
  if (ey_var < 0 && ey_var > -1e-10) ey_var <- 0

  # Implied standardized covariance matrix
  Sigma <- matrix(
    c(1,                    alpha,                 tau_prime + beta * alpha,
      alpha,                1,                     beta + alpha * tau_prime,
      tau_prime + beta*alpha, beta + alpha*tau_prime, 1),
    nrow = 3, byrow = TRUE
  )
  dimnames(Sigma) <- list(c("x","m","y"), c("x","m","y"))

  ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev < -1e-8)) stop("Implied covariance matrix is not positive semidefinite; check inputs.")

  if (isTRUE(empirical)) {
    return(list(
      Sigma      = Sigma,
      alpha      = alpha,
      beta       = beta,
      tau_prime  = tau_prime
    ))
  }

  # Otherwise we need n
  if (is.null(n)) stop("n must be provided when empirical = FALSE.")
  n <- as.integer(n)
  if (n < 3) stop("n must be >= 3.")
  if (abs(alpha) > 1) stop("|alpha| must be <= 1 (since cor(X,M) = alpha).")

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for empirical = FALSE.")
  }
  dat <- as.data.frame(MASS::mvrnorm(
    n = n, mu = c(0, 0, 0), Sigma = Sigma, empirical = FALSE
  ))
  names(dat) <- c("x", "m", "y")
  dat
}

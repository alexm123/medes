#' Upsilon Mediation Effect Estimates via SEM and/or OLS (MC CIs only)
#'
#' Computes Upsilon statistics for a single-mediator model using
#' SEM (\pkg{lavaan}) and/or OLS. Confidence intervals are obtained
#' via Monte Carlo simulation (parametric), not bootstrap.
#'
#' @param data A `data.frame` containing variables.
#' @param x,m,y Character strings: names of predictor, mediator, outcome.
#' @param engine "ols", "sem", or "both" (default).
#' @param R Integer. Number of Monte Carlo draws (default 5000).
#' @param ci_level Numeric in (0,1). Confidence level for MC intervals (default 0.95).
#' @param do_monte_carlo Logical. Whether to compute MC CIs (default TRUE).
#' @param stat "adjusted" (default), "unadjusted", or "both".
#'
#' @details
#' - **SEM MC CIs:** Draw (a, b) from a bivariate normal with mean equal to the
#'   SEM estimates and covariance from the SEM parameter covariance matrix;
#'   compute Upsilon on each draw with the observed X/Y variance ratio held fixed.
#'   (No model refits; typically very fast.)
#' - **OLS MC CIs:** Parametric simulation: generate new \code{m} and \code{y}
#'   given \code{x} from the fitted OLS models with normal errors, then refit the
#'   same OLS models on the simulated data to get (a, b). This is much lighter
#'   than a nonparametric bootstrap and avoids degeneracy from resampling rows.
#'
#' The "adjusted" Upsilon uses the usual variance corrections based on the
#' fitted-model variance estimates (treated as fixed during MC for SEM; recomputed
#' each draw for OLS).
#'
#' @return A named list including some or all of:
#' \describe{
#'   \item{upsilon_sem, adj.upsilon_sem}{Point estimates (if SEM selected).}
#'   \item{ci.upsilon_sem, ci.adj.upsilon_sem}{MC CIs for SEM (if requested).}
#'   \item{upsilon_ols, adj.upsilon_ols}{Point estimates (if OLS selected).}
#'   \item{ci.upsilon_ols, ci.adj.upsilon_ols}{MC CIs for OLS (if requested).}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' d <- data.frame(x = rnorm(200), m = rnorm(200), y = rnorm(200))
#'
#' # OLS only, MC CIs
#' res_ols <- upsilons_mc(d, x="x", m="m", y="y", engine="ols", R=2000)
#'
#' # SEM only, MC CIs (requires lavaan)
#' res_sem <- upsilons_mc(d, x="x", m="m", y="y", engine="sem", R=2000)
#'
#' # Both engines (default)
#' res_both <- upsilons_mc(d, x="x", m="m", y="y", R=2000)
#'
#' # Adjusted only (default) vs unadjusted vs both
#' res_adj  <- upsilons_mc(d, x="x", m="m", y="y", stat="adjusted")
#' res_un   <- upsilons_mc(d, x="x", m="m", y="y", stat="unadjusted")
#' res_both_stats <- upsilons_mc(d, x="x", m="m", y="y", stat="both")
#' }
#' @export
upsilons_mc <- function(data, x, m, y,
                        engine = c("both", "sem", "ols"),
                        R = 5000,
                        ci_level = 0.95,
                        do_monte_carlo = TRUE,
                        stat = c("adjusted", "unadjusted", "both")) {
  engine <- match.arg(engine)
  stat   <- match.arg(stat)
  do_sem <- engine %in% c("both", "sem")
  do_ols <- engine %in% c("both", "ols")

  stopifnot(is.data.frame(data))
  vars_needed <- c(x, m, y)
  missing_vars <- setdiff(vars_needed, names(data))
  if (length(missing_vars) > 0) {
    stop("These variables are missing from `data`: ",
         paste(missing_vars, collapse = ", "))
  }
  if (do_sem && !requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required for SEM but is not installed.")
  }
  if (!is.numeric(R) || R < 1) stop("`R` must be a positive integer.")
  if (!is.numeric(ci_level) || ci_level <= 0 || ci_level >= 1)
    stop("`ci_level` must be in (0, 1).")
  alpha <- 1 - ci_level

  # Consistent listwise deletion for OLS computations
  keeps <- stats::complete.cases(data[, c(x, m, y)])
  data_cc <- data[keeps, , drop = FALSE]
  S_cc <- stats::cov(data_cc[, c(x, y, m)])
  scov_ratio <- S_cc[x, x] / S_cc[y, y]

  out <- list()

  # ---------- SEM ----------
  if (do_sem) {
    med.model <- paste0(m, " ~ ", x, "\n",
                        y, " ~ ", x, " + ", m)
    med.res <- lavaan::sem(med.model, data)

    a_name <- paste0(m, "~", x)
    b_name <- paste0(y, "~", m)

    coefs  <- lavaan::coef(med.res)
    vc_sem <- lavaan::vcov(med.res)

    a_sem <- as.numeric(coefs[a_name])
    b_sem <- as.numeric(coefs[b_name])

    upsilon_sem <- (a_sem^2 * b_sem^2) * scov_ratio
    adj_sem <- ((a_sem^2 - vc_sem[a_name, a_name]) *
                  (b_sem^2 - vc_sem[b_name, b_name]) * scov_ratio)

    if (stat %in% c("unadjusted", "both")) out$upsilon_sem <- upsilon_sem
    if (stat %in% c("adjusted", "both"))   out$`adj.upsilon_sem` <- adj_sem

    if (do_monte_carlo) {
      # Draw (a,b) ~ MVN(mean, cov) without requiring MASS
      mu <- c(a_sem, b_sem)
      Sigma <- vc_sem[c(a_name, b_name), c(a_name, b_name), drop = FALSE]
      # Stable Cholesky (fallback to eigen if needed)
      chol_safe <- function(S) {
        ok <- TRUE
        out <- try(chol(S), silent = TRUE)
        if (inherits(out, "try-error")) {
          ok <- FALSE
          ev  <- eigen((S + t(S))/2, symmetric = TRUE)
          ev$values[ev$values < .Machine$double.eps] <- .Machine$double.eps
          out <- chol(ev$vectors %*% diag(ev$values, nrow(S)) %*% t(ev$vectors))
        }
        out
      }
      C <- chol_safe(Sigma)
      Z <- matrix(stats::rnorm(2 * R), ncol = 2)
      AB <- sweep(Z %*% t(C), 2, mu, FUN = "+")
      a_draw <- AB[,1]; b_draw <- AB[,2]

      ups_draw <- (a_draw^2 * b_draw^2) * scov_ratio
      adj_draw <- ((a_draw^2 - Sigma[1,1]) * (b_draw^2 - Sigma[2,2])) * scov_ratio

      if (stat %in% c("unadjusted", "both")) {
        out$ci.upsilon_sem <- stats::quantile(ups_draw, c(alpha/2, 1 - alpha/2), na.rm = TRUE)
      }
      if (stat %in% c("adjusted", "both")) {
        out$`ci.adj.upsilon_sem` <- stats::quantile(adj_draw, c(alpha/2, 1 - alpha/2), na.rm = TRUE)
      }
    }
  }

  # ---------- OLS ----------
  if (do_ols) {
    # Point estimates
    fit_a <- stats::lm(stats::as.formula(paste(m, "~", x)), data = data_cc)
    fit_b <- stats::lm(stats::as.formula(paste(y, "~", x, "+", m)), data = data_cc)

    a_ols <- stats::coef(fit_a)[x]
    b_ols <- stats::coef(fit_b)[m]

    vc_a  <- stats::vcov(fit_a)[x, x]
    vc_b  <- stats::vcov(fit_b)[m, m]

    upsilon_ols <- (a_ols^2 * b_ols^2) * scov_ratio
    adj_ols     <- ((a_ols^2 - vc_a) * (b_ols^2 - vc_b)) * scov_ratio

    if (stat %in% c("unadjusted", "both")) out$upsilon_ols <- as.numeric(upsilon_ols)
    if (stat %in% c("adjusted", "both"))   out$`adj.upsilon_ols` <- as.numeric(adj_ols)

    if (do_monte_carlo) {
      # Parametric simulation of data given x; refit simple OLS on each draw
      n <- nrow(data_cc)
      if (n < 3) stop("Not enough complete cases after listwise deletion.")

      # Extract fitted-model pieces
      a0 <- stats::coef(fit_a)[["(Intercept)"]]
      sig2_M <- summary(fit_a)$sigma^2

      b0 <- stats::coef(fit_b)[["(Intercept)"]]
      bx <- stats::coef(fit_b)[x]
      sig2_Y <- summary(fit_b)$sigma^2

      x_vec <- data_cc[[x]]

      # Storage
      ups_store <- adj_store <- numeric(R)

      for (r in seq_len(R)) {
        m_sim <- a0 + a_ols * x_vec + stats::rnorm(n, sd = sqrt(sig2_M))
        y_sim <- b0 + bx * x_vec + b_ols * m_sim + stats::rnorm(n, sd = sqrt(sig2_Y))

        # Refit on simulated data
        fit_a_s <- stats::lm(m_sim ~ x_vec)
        fit_b_s <- stats::lm(y_sim ~ x_vec + m_sim)

        a_s <- stats::coef(fit_a_s)[["x_vec"]]
        b_s <- stats::coef(fit_b_s)[["m_sim"]]

        vc_a_s <- stats::vcov(fit_a_s)[["x_vec","x_vec"]]
        vc_b_s <- stats::vcov(fit_b_s)[["m_sim","m_sim"]]

        S_sim <- stats::cov(cbind(x = x_vec, y = y_sim, m = m_sim))
        ratio_s <- S_sim["x","x"] / S_sim["y","y"]

        ups_store[r] <- (a_s^2 * b_s^2) * ratio_s
        adj_store[r] <- ((a_s^2 - vc_a_s) * (b_s^2 - vc_b_s)) * ratio_s
      }

      if (stat %in% c("unadjusted", "both")) {
        out$ci.upsilon_ols <- stats::quantile(ups_store, c(alpha/2, 1 - alpha/2), na.rm = TRUE)
      }
      if (stat %in% c("adjusted", "both")) {
        out$`ci.adj.upsilon_ols` <- stats::quantile(adj_store, c(alpha/2, 1 - alpha/2), na.rm = TRUE)
      }
    }
  }

  out
}

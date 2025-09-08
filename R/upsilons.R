#' Upsilon Mediation Effect Estimates via SEM and/or OLS
#'
#' Computes the Upsilon statistic for mediation effects using
#' Structural Equation Modeling (SEM, via \pkg{lavaan}) and/or
#' Ordinary Least Squares (OLS). Optionally, percentile bootstrap
#' confidence intervals are provided for the chosen engine(s).
#'
#' @param data A `data.frame` containing the variables of interest.
#' @param x Character string. Name of the independent (predictor) variable.
#' @param m Character string. Name of the mediator variable.
#' @param y Character string. Name of the dependent (outcome) variable.
#' @param engine Character. One of `"ols"`, `"sem"`, or `"both"` (default).
#'   Controls which estimator(s) are computed and, if bootstrapping,
#'   which bootstrap(s) are run.
#' @param R Integer. Number of bootstrap replications. Default is 5000.
#' @param ci_level Numeric. Confidence level for bootstrap intervals
#'   (between 0 and 1). Default is 0.95.
#' @param do_bootstrap Logical. Whether to perform bootstrapping
#'   to compute confidence intervals. Default is `TRUE`.
#'
#' @details
#' The function fits a simple mediation model with one mediator using
#' SEM (\pkg{lavaan}) and/or two OLS regressions. Upsilon statistics
#' (original and adjusted) are computed from path coefficients and covariance.
#'
#' If `do_bootstrap = TRUE`:
#' - SEM bootstrap is performed via [lavaan::bootstrapLavaan()] (if `engine` includes `"sem"`).
#' - OLS bootstrap uses row-resampling (pairs bootstrap) (if `engine` includes `"ols"`).
#'
#' The returned list includes only elements for the requested `engine`.
#'
#' @return A named `list` containing some or all of:
#' \describe{
#'   \item{upsilon_sem}{Point estimate of Upsilon from SEM.}
#'   \item{adj.upsilon_sem}{Adjusted Upsilon from SEM.}
#'   \item{ci.upsilon_sem}{Bootstrap CI for SEM Upsilon (if bootstrapped).}
#'   \item{ci.adj.upsilon_sem}{Bootstrap CI for adjusted SEM Upsilon (if bootstrapped).}
#'   \item{upsilon_ols}{Point estimate of Upsilon from OLS.}
#'   \item{adj.upsilon_ols}{Adjusted Upsilon from OLS.}
#'   \item{ci.upsilon_ols}{Bootstrap CI for OLS Upsilon (if bootstrapped).}
#'   \item{ci.adj.upsilon_ols}{Bootstrap CI for adjusted OLS Upsilon (if bootstrapped).}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' d <- data.frame(x = rnorm(100), m = rnorm(100), y = rnorm(100))
#'
#' # OLS only
#' res_ols <- upsilons(d, x="x", m="m", y="y", engine="ols",
#'                             do_bootstrap=FALSE)
#'
#' # SEM only with bootstrap
#' res_sem <- upsilons(d, x="x", m="m", y="y", engine="sem",
#'                             do_bootstrap=TRUE, R=1000)
#'
#' # Both (default)
#' res_both <- upsilons(d, x="x", m="m", y="y")
#' }
#'
#' @export
upsilons <- function(data, x, m, y,
                     engine = c("both", "sem", "ols"),
                     R = 5000,
                     ci_level = 0.95,
                     do_bootstrap = TRUE) {
  engine <- match.arg(engine)
  do_sem <- engine %in% c("both", "sem")
  do_ols <- engine %in% c("both", "ols")

  # Basic input checks
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

  # Build mediation model (SEM syntax)
  med.model <- paste0(m, " ~ ", x, "\n",
                      y, " ~ ", x, " + ", m)

  out <- list()

  # ---------- SEM (lavaan) ----------
  if (do_sem) {
    med.res <- lavaan::sem(med.model, data)

    a_sem <- lavaan::coef(med.res)[paste0(m, "~", x)]
    b_sem <- lavaan::coef(med.res)[paste0(y, "~", m)]
    vc_sem <- lavaan::vcov(med.res)

    # Sample covariance (SEM's notion)
    scov_sem <- lavaan::inspectSampleCov(med.model, data)$cov

    upsilon_sem <- (a_sem^2 * b_sem^2) * scov_sem[x, x] / scov_sem[y, y]
    adj.upsilon_sem <- ((a_sem^2 - vc_sem[paste0(m, "~", x), paste0(m, "~", x)]) *
                          (b_sem^2 - vc_sem[paste0(y, "~", m), paste0(y, "~", m)]) *
                          scov_sem[x, x] / scov_sem[y, y])

    out$upsilon_sem <- as.numeric(upsilon_sem)
    out$`adj.upsilon_sem` <- as.numeric(adj.upsilon_sem)
  }

  # ---------- OLS ----------
  if (do_ols) {
    fit_a <- stats::lm(stats::as.formula(paste(m, "~", x)), data = data)
    fit_b <- stats::lm(stats::as.formula(paste(y, "~", x, "+", m)), data = data)

    a_ols <- stats::coef(fit_a)[x]
    b_ols <- stats::coef(fit_b)[m]

    vc_a_ols <- stats::vcov(fit_a)[x, x]
    vc_b_ols <- stats::vcov(fit_b)[m, m]

    scov_ols <- stats::cov(data[, c(x, y, m)], use = "pairwise.complete.obs")

    upsilon_ols <- (a_ols^2 * b_ols^2) * scov_ols[x, x] / scov_ols[y, y]
    adj.upsilon_ols <- ((a_ols^2 - vc_a_ols) *
                          (b_ols^2 - vc_b_ols) *
                          scov_ols[x, x] / scov_ols[y, y])

    out$upsilon_ols <- as.numeric(upsilon_ols)
    out$`adj.upsilon_ols` <- as.numeric(adj.upsilon_ols)
  }

  # ---------- Bootstrap ----------
  if (do_bootstrap) {
    alpha <- 1 - ci_level

    # --- SEM bootstrap via lavaan ---
    if (do_sem) {
      lavaan.med.boot.fun <- function(out){
        data.boot <- lavaan::lavInspect(out, what = "data")
        colnames(data.boot) <- lavaan::lavNames(out)

        scov.boot <- lavaan::inspectSampleCov(med.model, data.boot)$cov

        a.boot  <- lavaan::coef(out)[paste0(m, "~", x)]
        b.boot  <- lavaan::coef(out)[paste0(y, "~", m)]
        vc.boot <- lavaan::vcov(out)

        upsilon.boot <- (a.boot^2 * b.boot^2) * scov.boot[x, x] / scov.boot[y, y]

        adj.upsilon.boot <- ((a.boot^2 - vc.boot[paste0(m, "~", x),
                                                 paste0(m, "~", x)]) *
                               (b.boot^2 - vc.boot[paste0(y, "~", m),
                                                   paste0(y, "~", m)]) *
                               scov.boot[x, x] / scov.boot[y, y])

        as.vector(c(upsilon.boot, adj.upsilon.boot))
      }

      if (do_ols) {
        message("Bootstrapping (SEM + OLS) may take several minutes...")
      } else {
        message("Bootstrapping (SEM) may take several minutes...")
      }

      boot.med.res <- lavaan::bootstrapLavaan(
        med.res, R = R, type = "ordinary", FUN = lavaan.med.boot.fun
      )
      boot.sem.mat <- if (is.list(boot.med.res)) do.call(rbind, boot.med.res) else
        as.matrix(boot.med.res)

      out$ci.upsilon_sem <- stats::quantile(boot.sem.mat[, 1],
                                            probs = c(alpha/2, 1 - alpha/2),
                                            na.rm = TRUE)
      out$`ci.adj.upsilon_sem` <- stats::quantile(boot.sem.mat[, 2],
                                                  probs = c(alpha/2, 1 - alpha/2),
                                                  na.rm = TRUE)
    }

    # --- OLS bootstrap (pairs bootstrap on rows) ---
    if (do_ols) {
      if (!do_sem) message("Bootstrapping (OLS) may take several minutes...")

      n <- nrow(data)
      boot_ols <- replicate(R, {
        idx <- sample.int(n, n, replace = TRUE)
        db  <- data[idx, , drop = FALSE]

        a_val <- b_val <- NA_real_
        vc_a  <- vc_b  <- NA_real_

        # a-path
        fit_a_b <- try(stats::lm(stats::as.formula(paste(m, "~", x)), data = db),
                       silent = TRUE)
        if (!inherits(fit_a_b, "try-error")) {
          coefs <- try(stats::coef(fit_a_b), silent = TRUE)
          vcv   <- try(stats::vcov(fit_a_b), silent = TRUE)
          if (!inherits(coefs, "try-error") && !inherits(vcv, "try-error") &&
              !is.na(coefs[x]) && !is.na(vcv[x, x])) {
            a_val <- coefs[x]; vc_a <- vcv[x, x]
          }
        }

        # b-path
        fit_b_b <- try(stats::lm(stats::as.formula(paste(y, "~", x, "+", m)), data = db),
                       silent = TRUE)
        if (!inherits(fit_b_b, "try-error")) {
          coefs <- try(stats::coef(fit_b_b), silent = TRUE)
          vcv   <- try(stats::vcov(fit_b_b), silent = TRUE)
          if (!inherits(coefs, "try-error") && !inherits(vcv, "try-error") &&
              !is.na(coefs[m]) && !is.na(vcv[m, m])) {
            b_val <- coefs[m]; vc_b <- vcv[m, m]
          }
        }

        scov_b <- try(stats::cov(db[, c(x, y, m)], use = "pairwise.complete.obs"),
                      silent = TRUE)
        if (inherits(scov_b, "try-error") || anyNA(dim(scov_b))) {
          scov_b <- matrix(NA_real_, 3, 3,
                           dimnames = list(c(x, y, m), c(x, y, m)))
        }

        ups_b  <- (a_val^2 * b_val^2) * scov_b[x, x] / scov_b[y, y]
        adj_b  <- ((a_val^2 - vc_a) * (b_val^2 - vc_b)) *
          scov_b[x, x] / scov_b[y, y]
        c(ups_b, adj_b)
      })

      boot_ols <- t(boot_ols)
      out$ci.upsilon_ols <- stats::quantile(boot_ols[, 1],
                                            probs = c(alpha/2, 1 - alpha/2),
                                            na.rm = TRUE)
      out$`ci.adj.upsilon_ols` <- stats::quantile(boot_ols[, 2],
                                                  probs = c(alpha/2, 1 - alpha/2),
                                                  na.rm = TRUE)
    }
  }

  return(out)
}

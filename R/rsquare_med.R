#' Fairchild et al.'s R-squared Mediated Effect (OLS and/or SEM)
#'
#' Computes the proportion of variance in the outcome that is mediated
#' via a single mediator, following the approach of Fairchild et al.
#' Supports an OLS (correlation-based) engine and a SEM engine (via
#' \pkg{lavaan}) that uses model-implied correlations. Optionally
#' returns percentile bootstrap confidence intervals.
#'
#' If missing data are detected in \code{x}, \code{m}, or \code{y},
#' the function prints: "Hello! You have missing data. We recommend SEM,
#' which will use FIML." When SEM is requested in this case, it fits the
#' model with FIML (\code{missing = "fiml"}).
#'
#' @param data A `data.frame` containing the variables of interest.
#' @param x Character string. Name of the independent (predictor) variable.
#' @param m Character string. Name of the mediator variable.
#' @param y Character string. Name of the dependent (outcome) variable.
#' @param engine Character. One of `"ols"`, `"sem"`, or `"both"` (default).
#'   Controls which estimator(s) are computed and bootstrapped.
#' @param R Integer. Number of bootstrap replications (default = 5000).
#' @param ci_level Numeric. Confidence level for bootstrap intervals
#'   (between 0 and 1). Default is 0.95.
#' @param do_bootstrap Logical. Whether to perform bootstrapping to
#'   compute confidence intervals. Default is `TRUE`.
#'
#' @details
#' **OLS engine.** Uses observed pairwise correlations among \code{x}, \code{m}, \code{y}:
#' \deqn{R^2_{y\cdot xm} = \frac{r_{xy}^2 + r_{my}^2 - 2 r_{xy} r_{my} r_{xm}}{1 - r_{xm}^2}}
#' and \deqn{R^2_\mathrm{med} = r_{my}^2 - \left(R^2_{y\cdot xm} - r_{xy}^2\right).}
#' OLS uses listwise deletion.
#'
#' **SEM engine.** Fits the path model \code{m ~ x; y ~ x + m} in \pkg{lavaan},
#' obtains the model-implied covariance for \code{x,m,y}, converts to correlations,
#' and plugs them into the same Fairchild formula. If missing data are detected,
#' \pkg{lavaan} is run with \code{missing = "fiml"} (FIML).
#'
#' If `do_bootstrap = TRUE`:
#' - OLS uses a pairs bootstrap (row resampling) with \pkg{boot}.
#' - SEM uses [lavaan::bootstrapLavaan()] and recomputes the statistic from each
#'   bootstrapped fitted object (respecting the \code{missing} setting).
#'
#' @return
#' - If `engine = "ols"` and `do_bootstrap = FALSE`, returns a single numeric
#'   (backward-compatible).
#' - Otherwise, returns a named `list` with some or all of:
#'   \describe{
#'     \item{rsquaredmediated_ols}{Point estimate (OLS).}
#'     \item{rsquaredmediated_ci_lower_ols, rsquaredmediated_ci_upper_ols}{Percentile CI (OLS) if bootstrapped.}
#'     \item{rsquaredmediated_sem}{Point estimate (SEM).}
#'     \item{rsquaredmediated_ci_lower_sem, rsquaredmediated_ci_upper_sem}{Percentile CI (SEM) if bootstrapped.}
#'   }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' d <- data.frame(x = rnorm(200), m = rnorm(200), y = rnorm(200))
#'
#' # OLS only (numeric return, no bootstrap)
#' rsq_ols <- rsquare_med(d, x="x", m="m", y="y", engine="ols", do_bootstrap=FALSE)
#'
#' # SEM only with bootstrap
#' rsq_sem <- rsquare_med(d, x="x", m="m", y="y", engine="sem", R=1000)
#'
#' # Both engines
#' rsq_both <- rsquare_med(d, x="x", m="m", y="y", engine="both", R=1000)
#' }
#'
#' @export
rsquare_med <- function(data, x, m, y,
                        engine = c("both", "ols", "sem"),
                        R = 1000,
                        ci_level = 0.95,
                        do_bootstrap = TRUE) {
  engine <- match.arg(engine)
  do_ols <- engine %in% c("both", "ols")
  do_sem <- engine %in% c("both", "sem")

  stopifnot(is.data.frame(data))
  vars_needed <- c(x, m, y)
  missing_vars <- setdiff(vars_needed, names(data))
  if (length(missing_vars) > 0) {
    stop("These variables are missing from `data`: ",
         paste(missing_vars, collapse = ", "))
  }
  if (!is.numeric(R) || R < 1) stop("`R` must be a positive integer.")
  if (!is.numeric(ci_level) || ci_level <= 0 || ci_level >= 1)
    stop("`ci_level` must be in (0, 1).")
  if (do_sem && !requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required for SEM but is not installed.")
  }

  # Detect missingness in the analysis variables
  has_missing <- any(is.na(data[, vars_needed, drop = FALSE]))
  if (has_missing) {
    message("Hello! You have missing data. We recommend SEM, which will use FIML.")
  }

  # Helper: Fairchild mediated R^2 from correlations
  fairchild_rsquared <- function(rxm, rxy, rmy) {
    overall_r2 <- ((rxy^2 + rmy^2) - 2 * rxy * rmy * rxm) / (1 - rxm^2)
    rmy^2 - (overall_r2 - rxy^2)
  }

  out <- list()

  # ---------- OLS engine (observed correlations) ----------
  if (do_ols) {
    # Listwise deletion
    keeps <- stats::complete.cases(data[, c(x, m, y)])
    data_cc  <-  data[keeps, , drop = FALSE]
    if (nrow(data_cc) < 3) {
      stop("Not enough complete cases for OLS correlations after listwise deletion.")
    }

    rxm <- stats::cor(data_cc[[x]], data_cc[[m]])
    rxy <- stats::cor(data_cc[[x]], data_cc[[y]])
    rmy <- stats::cor(data_cc[[m]], data_cc[[y]])

    rsq_ols <- fairchild_rsquared(rxm, rxy, rmy)
    out$rsquaredmediated_ols <- as.numeric(rsq_ols)
  }

  # ---------- SEM engine (model-implied correlations) ----------
  if (do_sem) {
    med.model <- paste0(m, " ~ ", x, "\n", y, " ~ ", x, " + ", m)

    # Use FIML if missing data are present
    fit <- lavaan::sem(
      model = med.model,
      data  = data,
      missing = if (has_missing) "fiml" else "listwise",
      fixed.x = FALSE,
      meanstructure = if (has_missing) TRUE else FALSE
    )

    # Model-implied covariance among observed variables
    covlist <- lavaan::fitted(fit)$cov
    cov_imp <- if (is.list(covlist)) covlist[[1]] else covlist
    cov_imp <- cov_imp[c(x, m, y), c(x, m, y), drop = FALSE]
    cor_imp <- stats::cov2cor(cov_imp)

    rxm_sem <- cor_imp[x, m]
    rxy_sem <- cor_imp[x, y]
    rmy_sem <- cor_imp[m, y]

    rsq_sem <- fairchild_rsquared(rxm_sem, rxy_sem, rmy_sem)
    out$rsquaredmediated_sem <- as.numeric(rsq_sem)
  }

  # ---------- Bootstrap ----------
  if (do_bootstrap) {
    if (do_ols && do_sem) {
      message("Bootstrapping (SEM + OLS) may take several minutes...")
    } else if (do_ols) {
      message("Bootstrapping (OLS) may take several minutes...")
    } else if (do_sem) {
      message("Bootstrapping (SEM) may take several minutes...")
    }

    alpha <- 1 - ci_level

    # OLS bootstrap (pairs)
    if (do_ols) {
      boot_stat <- function(d, indices) {
        db <- d[indices, , drop = FALSE]
        db_cc <- db[stats::complete.cases(db[, c(x, m, y)]), , drop = FALSE]
        if (nrow(db_cc) < 3) return(NA_real_)
        rxm_b <- stats::cor(db_cc[[x]], db_cc[[m]])
        rxy_b <- stats::cor(db_cc[[x]], db_cc[[y]])
        rmy_b <- stats::cor(db_cc[[m]], db_cc[[y]])
        fairchild_rsquared(rxm_b, rxy_b, rmy_b)
      }
      boot_results <- boot::boot(data = data, statistic = boot_stat, R = R)
      # Percentile CI
      boot_ci <- boot::boot.ci(boot_results, conf = ci_level, type = "perc", index = 1)
      out$rsquaredmediated_ci_lower_ols <- unname(boot_ci$percent[4])
      out$rsquaredmediated_ci_upper_ols <- unname(boot_ci$percent[5])
    }

    # SEM bootstrap via lavaan (inherits missing="fiml" if used above)
    if (do_sem) {
      sem_boot_fun <- function(fit_obj) {
        covlist_b <- lavaan::fitted(fit_obj)$cov
        cov_imp_b <- if (is.list(covlist_b)) covlist_b[[1]] else covlist_b
        cov_imp_b <- cov_imp_b[c(x, m, y), c(x, m, y), drop = FALSE]
        cor_imp_b <- stats::cov2cor(cov_imp_b)
        fairchild_rsquared(cor_imp_b[x, m], cor_imp_b[x, y], cor_imp_b[m, y])
      }

      boot_sem <- lavaan::bootstrapLavaan(
        fit, R = R, type = "ordinary", FUN = sem_boot_fun
      )
      sem_vals <- if (is.list(boot_sem)) unlist(boot_sem) else as.numeric(boot_sem)
      qs <- stats::quantile(sem_vals, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
      out$rsquaredmediated_ci_lower_sem <- unname(qs[1])
      out$rsquaredmediated_ci_upper_sem <- unname(qs[2])
    }
  }

  # Back-compat: if user requests classic behavior (OLS only, no bootstrap), return numeric
  if (engine == "ols" && !do_bootstrap) {
    return(out$rsquaredmediated_ols)
  }

  return(out)
}

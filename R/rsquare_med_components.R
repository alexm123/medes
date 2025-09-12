#' Pairwise Correlations (OLS and/or SEM) with Optional Bootstrap CIs
#'
#' Returns the three correlations among \code{x}, \code{m}, \code{y}:
#' \eqn{r_{xm}}, \eqn{r_{xy}}, \eqn{r_{my}}. Supports an OLS
#' (observed correlations) engine and a SEM engine (via \pkg{lavaan})
#' using model-implied correlations from the path model \code{m ~ x; y ~ x + m}.
#' Optionally returns percentile bootstrap confidence intervals.
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
#' **OLS engine.** Uses observed pairwise correlations among \code{x}, \code{m}, \code{y}
#' after listwise deletion.
#'
#' **SEM engine.** Fits \code{m ~ x; y ~ x + m} in \pkg{lavaan}, obtains the
#' model-implied covariance for \code{x,m,y}, converts to correlations, and
#' returns those \eqn{r_{xm}}, \eqn{r_{xy}}, \eqn{r_{my}}.
#'
#' If `do_bootstrap = TRUE`:
#' - OLS uses a pairs bootstrap (row resampling) with \pkg{boot}, returning
#'   percentile CIs for each correlation.
#' - SEM uses [lavaan::bootstrapLavaan()] and recomputes the three correlations
#'   from each bootstrapped fitted object, returning percentile CIs.
#'
#' @return A named `list` containing some or all of:
#' \describe{
#'   \item{rxm_ols, rxy_ols, rmy_ols}{Point estimates (OLS), if requested by `engine`.}
#'   \item{rxm_ci_lower_ols, rxm_ci_upper_ols, ...}{Percentile CIs (OLS), if bootstrapped.}
#'   \item{rxm_sem, rxy_sem, rmy_sem}{Point estimates (SEM), if requested by `engine`.}
#'   \item{rxm_ci_lower_sem, rxm_ci_upper_sem, ...}{Percentile CIs (SEM), if bootstrapped.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' d <- data.frame(x = rnorm(200), m = rnorm(200), y = rnorm(200))
#'
#' # OLS only (no bootstrap)
#' out_ols <- rsquare_med_components(d, x="x", m="m", y="y",
#'                                   engine="ols", do_bootstrap=FALSE)
#'
#' # SEM only with bootstrap
#' out_sem <- rsquare_med_components(d, x="x", m="m", y="y",
#'                                   engine="sem", R=1000)
#'
#' # Both engines with bootstrap
#' out_both <- rsquare_med_components(d, x="x", m="m", y="y",
#'                                    engine="both", R=1000)
#' }
#'
#' @export
rsquare_med_components <- function(data, x, m, y,
                                   engine = c("both", "ols", "sem"),
                                   R = 5000,
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

  out <- list()

  # ---------- OLS engine (observed correlations) ----------
  if (do_ols) {
    keeps <- stats::complete.cases(data[, c(x, m, y)])
    data_cc  <-  data[keeps, , drop = FALSE]

    rxm <- stats::cor(data_cc[[x]], data_cc[[m]])
    rxy <- stats::cor(data_cc[[x]], data_cc[[y]])
    rmy <- stats::cor(data_cc[[m]], data_cc[[y]])



    out$rxm_ols <- as.numeric(rxm)
    out$rxy_ols <- as.numeric(rxy)
    out$rmy_ols <- as.numeric(rmy)
  }

  # ---------- SEM engine (model-implied correlations) ----------
  if (do_sem) {
    med.model <- paste0(m, " ~ ", x, "\n", y, " ~ ", x, " + ", m)
    fit <- lavaan::sem(med.model, data = data)

    covlist <- lavaan::fitted(fit)$cov
    cov_imp <- if (is.list(covlist)) covlist[[1]] else covlist
    cov_imp <- cov_imp[c(x, m, y), c(x, m, y), drop = FALSE]
    cor_imp <- stats::cov2cor(cov_imp)

    rxm_sem <- unname(cor_imp[x, m])
    rxy_sem <- unname(cor_imp[x, y])
    rmy_sem <- unname(cor_imp[m, y])

    out$rxm_sem <- as.numeric(rxm_sem)
    out$rxy_sem <- as.numeric(rxy_sem)
    out$rmy_sem <- as.numeric(rmy_sem)
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
      stat_fun <- function(d, indices) {
        db <- d[indices, , drop = FALSE]
        keeps_b <- stats::complete.cases(db[, c(x, m, y)])
        db_cc <- db[keeps_b, , drop = FALSE]
        if (nrow(db_cc) < 3) return(c(NA_real_, NA_real_, NA_real_))
        rxm_b <- stats::cor(db_cc[[x]], db_cc[[m]])
        rxy_b <- stats::cor(db_cc[[x]], db_cc[[y]])
        rmy_b <- stats::cor(db_cc[[m]], db_cc[[y]])
        c(rxm_b, rxy_b, rmy_b)
      }
      # Only pass the needed columns to speed up boot
      dat_small <- data[, c(x, m, y)]
      boot_results <- boot::boot(data = dat_small, statistic = stat_fun, R = R)
      # Percentile CIs per correlation
      qs <- apply(boot_results$t, 2, stats::quantile,
                  probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE, names = FALSE)
      out$rxm_ci_lower_ols <- unname(qs[1, 1]); out$rxm_ci_upper_ols <- unname(qs[2, 1])
      out$rxy_ci_lower_ols <- unname(qs[1, 2]); out$rxy_ci_upper_ols <- unname(qs[2, 2])
      out$rmy_ci_lower_ols <- unname(qs[1, 3]); out$rmy_ci_upper_ols <- unname(qs[2, 3])
    }

    # SEM bootstrap via lavaan
    if (do_sem) {
      sem_boot_fun <- function(fit_obj) {
        covlist_b <- lavaan::fitted(fit_obj)$cov
        cov_imp_b <- if (is.list(covlist_b)) covlist_b[[1]] else covlist_b
        cov_imp_b <- cov_imp_b[c(x, m, y), c(x, m, y), drop = FALSE]
        cor_imp_b <- stats::cov2cor(cov_imp_b)
        c(unname(cor_imp_b[x, m]),
          unname(cor_imp_b[x, y]),
          unname(cor_imp_b[m, y]))
      }

      boot_sem <- lavaan::bootstrapLavaan(
        lavaan::sem(paste0(m, " ~ ", x, "\n", y, " ~ ", x, " + ", m), data = data),
        R = R, type = "ordinary", FUN = sem_boot_fun
      )

      sem_mat <- if (is.matrix(boot_sem)) {
        boot_sem
      } else if (is.list(boot_sem)) {
        do.call(rbind, boot_sem)
      } else {
        as.matrix(boot_sem)
      }

      qs_sem <- apply(sem_mat, 2, stats::quantile,
                      probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE, names = FALSE)
      out$rxm_ci_lower_sem <- unname(qs_sem[1, 1]); out$rxm_ci_upper_sem <- unname(qs_sem[2, 1])
      out$rxy_ci_lower_sem <- unname(qs_sem[1, 2]); out$rxy_ci_upper_sem <- unname(qs_sem[2, 2])
      out$rmy_ci_lower_sem <- unname(qs_sem[1, 3]); out$rmy_ci_upper_sem <- unname(qs_sem[2, 3])
    }
  }

  return(out)
}

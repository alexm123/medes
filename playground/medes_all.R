#' Comprehensive Mediation Effect Sizes via SEM (data or sample moments)
#'
#' Fits the simple mediation model \code{m ~ x; y ~ x + m} using \pkg{lavaan}
#' either from raw \code{data} or from supplied sample moments
#' (\code{sample.cov}, optional \code{sample.mean}, and \code{sample.nobs}).
#' Always returns BOTH unadjusted and adjusted Upsilon statistics, plus the
#' effect-size measures described by Fairchild et al. (2009), all computed from
#' the SEM model-implied correlations so they are mutually consistent.
#' Optionally returns percentile bootstrap CIs via \code{lavaan::bootstrapLavaan()}.
#'
#' A short message is printed to the console: "Using SEM (lavaan).
#' Returning BOTH unadjusted and adjusted Upsilons."
#'
#' @param data A `data.frame` with variables \code{x}, \code{m}, \code{y}.
#'   Ignored if \code{sample.cov} is provided.
#' @param x Character. Name of the predictor.
#' @param m Character. Name of the mediator.
#' @param y Character. Name of the outcome.
#' @param R Integer. Number of bootstrap replications (default 5000).
#' @param ci_level Numeric in (0,1). Confidence level for percentile CIs (default 0.95).
#' @param do_bootstrap Logical. If `TRUE`, percentile CIs are returned (default `FALSE`).
#' @param sample.cov Optional covariance matrix for variables \code{c(x,m,y)}.
#'   Must be named and include rows/cols for \code{x,m,y}.
#' @param sample.mean Optional named mean vector for \code{c(x,m,y)}.
#'   Only used with \code{sample.cov}.
#' @param sample.nobs Optional integer sample size used with \code{sample.cov}.
#'
#' @details
#' **SEM fit.** If raw \code{data} contain missing values, the model is fit with
#' FIML (\code{missing = "fiml"}), otherwise listwise. When using sample moments,
#' \code{sample.cov} (and optionally \code{sample.mean}) with \code{sample.nobs}
#' are passed directly to \code{lavaan::sem()}.
#'
#' **Upsilon statistics.** Both unadjusted and adjusted Upsilons are returned.
#' These use the path coefficients \eqn{a} and \eqn{b}, their sampling variances,
#' and the *sample* covariance used by the fit.
#'
#' **Fairchild et al. (2009) measures** (computed from SEM model-implied
#' correlations among \code{x,m,y}):
#' \itemize{
#'   \item \code{r2_med}: second-order common (mediated) effect in \eqn{Y}.
#'   \item \code{r2_xm}: variance in \eqn{M} explained by \eqn{X} (\eqn{r_{XM}^2}).
#'   \item \code{r2_my_x}: squared partial correlation \eqn{r^2_{MY\cdot X}}.
#'   \item \code{r2_xy_m}: squared partial correlation \eqn{r^2_{XY\cdot M}}.
#'   \item \code{r2_overall_y_xm}: overall multiple \eqn{R^2} of \eqn{Y} on \eqn{X,M}.
#'   \item \code{prop_r2_mediated}: proportion of explained outcome variance that is
#'     mediated: \code{r2_med / r2_overall_y_xm}.
#' }
#' Standardized path coefficients (\code{a_std}, \code{b_std}, \code{tauprime_std})
#' are also returned for convenience. Note that \code{r2_med} can be negative
#' (suppression), as discussed by Fairchild et al.
#'
#' **Bootstrapping.** If \code{do_bootstrap = TRUE}, percentile CIs are returned
#' for: \code{upsilon_sem}, \code{adj.upsilon_sem}, \code{r2_med},
#' \code{r2_overall_y_xm}, \code{prop_r2_mediated}, \code{r2_xm},
#' \code{r2_my_x}, and \code{r2_xy_m}.
#'
#' @return A named `list` containing:
#' \describe{
#'   \item{upsilon_sem}{Upsilon (unadjusted).}
#'   \item{adj.upsilon_sem}{Upsilon (adjusted).}
#'   \item{ci.upsilon_sem}{Percentile CI for \code{upsilon_sem} (if bootstrapped).}
#'   \item{ci.adj.upsilon_sem}{Percentile CI for \code{adj.upsilon_sem} (if bootstrapped).}
#'   \item{r2_med}{Fairchild mediated \eqn{R^2}.}
#'   \item{r2_xm}{\eqn{r_{XM}^2}.}
#'   \item{r2_my_x}{\eqn{r^2_{MY\cdot X}}.}
#'   \item{r2_xy_m}{\eqn{r^2_{XY\cdot M}}.}
#'   \item{r2_overall_y_xm}{Overall \eqn{R^2} of \eqn{Y \sim X + M}.}
#'   \item{prop_r2_mediated}{\code{r2_med / r2_overall_y_xm}.}
#'   \item{a_std, b_std, tauprime_std}{Standardized path coefficients.}
#'   \item{ci.*}{Percentile CIs corresponding to the above (if bootstrapped).}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' d <- data.frame(x = rnorm(300), m = rnorm(300), y = rnorm(300))
#'
#' # Raw data, point estimates only
#' out <- medes_all(d, x="x", m="m", y="y", do_bootstrap = FALSE)
#'
#' # With bootstrap CIs
#' out_ci <- medes_all(d, "x","m","y", R = 1000, do_bootstrap = TRUE)
#'
#' # From sample moments
#' S  <- cov(d[,c("x","m","y")])
#' Mu <- colMeans(d[,c("x","m","y")])
#' out_mom <- medes_all(
#'   data = NULL, x="x", m="m", y="y",
#'   sample.cov = S, sample.mean = Mu, sample.nobs = nrow(d),
#'   do_bootstrap = FALSE
#' )
#' }
#'
#' @importFrom stats vcov cov2cor quantile
#' @importFrom lavaan sem fitted standardizedSolution lavInspect
#' @export
medes_all <- function(data, x, m, y,
                      R = 5000,
                      ci_level = 0.95,
                      do_bootstrap = FALSE) {

  out <- list()
  vars_needed <- c(x, m, y)

  # Build mediation model (SEM syntax)
  med.model <- paste0(m, " ~ ", x, "\n",
                      y, " ~ ", x, " + ", m)

  # Detect missingness in the analysis variables
  has_missing <- any(is.na(data[, vars_needed, drop = FALSE]))

  fit <- lavaan::sem(
    model = med.model,
    data  = data,
    missing = if (has_missing) "fiml" else "listwise",
    fixed.x = FALSE,
    meanstructure = if (has_missing) TRUE else FALSE
  )

  ############### Compute upsilons #############################
  a <- lavaan::coef(fit)[paste0(m, "~", x)]
  b <- lavaan::coef(fit)[paste0(y, "~", m)]
  vc <- lavaan::vcov(fit)

  # Sample covariance (SEM's notion)
  scov <- lavaan::inspectSampleCov(med.model, data)$cov

  upsilon <- (a^2 * b^2) * scov[x, x] / scov[y, y]
  adj.upsilon <- ((a^2 - vc[paste0(m, "~", x), paste0(m, "~", x)]) *
                        (b^2 - vc[paste0(y, "~", m), paste0(y, "~", m)]) *
                        scov[x, x] / scov[y, y])

  out$upsilon <- as.numeric(upsilon)
  out$adj.upsilon <- as.numeric(adj.upsilon)
  #######################################################
  ########### Compute R2 and Components ##################

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

  rxm <- cor_imp[x, m]
  rxy <- cor_imp[x, y]
  rmy <- cor_imp[m, y]

  overall_r2 <- ((rxy^2 + rmy^2) - 2 * rxy * rmy * rxm) / (1 - rxm^2)
  rsquared_mediated <- rmy^2 - (overall_r2 - rxy^2)
  out$rsquared_mediated <- as.numeric(rsquared_mediated)
  out$overall_r2 <- as.numeric(overall_r2)
  out$rxm<- as.numeric(rxm)
  out$rxy<- as.numeric(rxy)
  out$rmy<- as.numeric(rmy)

  return(out)
}

medes_all(d, x = "x", m = "m", y = "y")
upsilons(d, x = "x", m = "m", y = "y", do_bootstrap = FALSE)
rsquare_med(d, x = "x", m = "m", y = "y", do_bootstrap = FALSE)

rsquare_med(d, x = "x", m = "m", y = "y", do_bootstrap = FALSE)








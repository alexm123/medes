#' Fairchild et al.'s R-squared Mediated Effect
#'
#' Computes the proportion of variance in the outcome that is mediated
#' via a single mediator, following the approach of Fairchild et al.
#' Optionally returns bootstrap confidence intervals.
#'
#' @param data A `data.frame` containing the variables of interest.
#' @param x Character string. Name of the independent (predictor) variable.
#' @param m Character string. Name of the mediator variable.
#' @param y Character string. Name of the dependent (outcome) variable.
#' @param R Integer. Number of bootstrap replications (default = 5000).
#' @param ci_level Numeric. Confidence level for bootstrap intervals
#'   (between 0 and 1). Default is 0.95.
#' @param do_bootstrap Logical. Whether to perform bootstrapping to
#'   compute confidence intervals. Default is `TRUE`.
#'
#' @details
#' The mediated R-squared is calculated from pairwise correlations among
#' the predictor, mediator, and outcome variables. When bootstrapping is
#' enabled, the function resamples rows with replacement and recomputes
#' the statistic, returning percentile confidence intervals.
#'
#' @return
#' If `do_bootstrap = FALSE`, a single numeric value:
#' \item{rsquaredmediated}{The estimated mediated R-squared.}
#'
#' If `do_bootstrap = TRUE`, a list with elements:
#' \describe{
#'   \item{rsquaredmediated_ci_lower}{Lower percentile CI bound.}
#'   \item{rsquaredmediated}{Point estimate of mediated R-squared.}
#'   \item{rsquaredmediated_ci_upper}{Upper percentile CI bound.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' d <- data.frame(
#'   x = rnorm(100),
#'   m = rnorm(100),
#'   y = rnorm(100)
#' )
#'
#' # Without bootstrap
#' rsquare_med(d, x = "x", m = "m", y = "y", do_bootstrap = FALSE)
#'
#' # With bootstrap (small R for speed)
#' rsquare_med(d, x = "x", m = "m", y = "y", R = 100, do_bootstrap = TRUE)
#' }
#'
#' @export
rsquare_med <- function(data, x, m, y,
                        R = 5000,
                        ci_level = 0.95,
                        do_bootstrap = TRUE) {
  # Compute correlations among the variables
  rxm <- cor(data[[x]], data[[m]])
  rxy <- cor(data[[x]], data[[y]])
  rmy <- cor(data[[m]], data[[y]])

  overallrsquared <- (((rxy^2) + (rmy^2)) - (2 * rxy * rmy * rxm)) / (1 - rxm^2)

  rsquaredmediated <- (rmy^2) - (overallrsquared - (rxy^2))

  # Bootstrapping if requested
  if (do_bootstrap) {
    boot_stat <- function(d, indices) {
      d_boot <- d[indices, ]
      rxm_b <- cor(d_boot[[x]], d_boot[[m]])
      rxy_b <- cor(d_boot[[x]], d_boot[[y]])
      rmy_b <- cor(d_boot[[m]], d_boot[[y]])
      rxmsquared_b <- rxm_b^2
      overallrsquared_b <- (((rxy_b^2) + (rmy_b^2)) - (2 * rxy_b * rmy_b * rxm_b)) / (1 - rxmsquared_b)

      rsquaredmediated_b <- (rmy_b^2) - (overallrsquared_b - (rxy_b^2))

      return(rsquaredmediated_b)
    }

    boot_results <- boot::boot(data = data, statistic = boot_stat, R = R)
    boot_ci_rsquaredmediated <- boot::boot.ci(boot_results,
                                              conf = ci_level, type = "perc", index = 1)

    rsquaredmediated_ci_lower <- boot_ci_rsquaredmediated$percent[4]
    rsquaredmediated_ci_upper <- boot_ci_rsquaredmediated$percent[5]

    results = list(
      rsquaredmediated_ci_lower = rsquaredmediated_ci_lower,
      rsquaredmediated = rsquaredmediated,
      rsquaredmediated_ci_upper = rsquaredmediated_ci_upper
    )
    return(results)
  }

  return(rsquaredmediated = rsquaredmediated)
}

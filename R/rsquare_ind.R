#' # Function to compute de Heus' R squared Indirect
#'
#' Computes the de Heus R^2 for an indirect (mediation) pathway with optional
#' percentile bootstrap CIs.
#'
#' @param data A data.frame containing the variables.
#' @param x,m,y Character names of the X (predictor), M (mediator), and Y (outcome) columns in `data`.
#' @param R Integer, bootstrap replications. Default 5000.
#' @param ci_level Confidence level for the percentile CI. Default 0.95.
#' @param do_bootstrap Logical; if `TRUE`, returns bootstrap percentile CI.
#'
#' @return If `do_bootstrap = FALSE`, a single numeric `rsq_ind`.
#' If `TRUE`, a list with `rsq_ind_lower`, `rsq_ind`, and `rsq_ind_upper`.
#'
#' @examples
#' \dontrun{
#'   set.seed(42)
#'   n <- 200
#'   x <- rnorm(n, mean = 50, sd = 10)
#'   m  <- 0.5*x + rnorm(n, mean = 0, sd = 5)
#'   y <- 0.3*x + 0.6*m + rnorm(n, mean = 0, sd = 5)
#'   d <- data.frame(x, m, y)
#'   rsq_ind(d, "x", "m", "y", R = 1000)
#' }
#'
#' @export
#' @importFrom stats cor lm as.formula

rsq_indirect <- function(data, x, m, y,
                        R = 5000,
                        ci_level = 0.95,
                        do_bootstrap = TRUE) {

  rxm <- cor(data[[x]], data[[m]])

  # Regression: y ~ x + m (for tau_prime and beta)
  model2 <- lm(as.formula(paste(y, "~", x, "+", m)), data = data)
  beta_YX_M <- QuantPsyc::lm.beta(model2)[[x]]

  r2_y_x.m <- (1 - rxm^2) * beta_YX_M
  rsq_ind <- r2_y_x.m * (rxm)^2

  if (do_bootstrap) {
    boot_stat <- function(d, indices) {
      d_boot <- d[indices, ]
      rxm_b <- cor(d_boot[[x]], d_boot[[m]])
      model2_b <- lm(as.formula(paste(y, "~", x, "+", m)), data = d_boot)

      beta_YX_M_b <- QuantPsyc::lm.beta(model2_b)[[x]]

      r2_y_x.m_b <- (1 - rxm_b^2) * beta_YX_M_b
      rsq_ind_b <- r2_y_x.m_b * (rxm_b)^2

      return(rsq_ind_b)
    }

    boot_results <- boot::boot(data = data, statistic = boot_stat, R = R)
    boot_ci_rsq_ind <- boot::boot.ci(boot_results,
                                              conf = ci_level,
                                              type = "perc",
                                              index = 1)

    rsq_ind_lower <- boot_ci_rsq_ind$percent[4]
    rsq_ind_upper <- boot_ci_rsq_ind$percent[5]

    results = list(rsq_ind_lower = rsq_ind_lower,
                   rsq_ind = rsq_ind,
                   rsq_ind_upper = rsq_ind_upper)
    return(results)
  }

  return(rsq_ind = rsq_ind)
}

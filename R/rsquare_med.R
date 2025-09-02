# Function to compute Fairchild et al.'s R^2 mediated

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

  results = list(rsquaredmediated_ci_lower = rsquaredmediated_ci_lower,
                 rsquaredmediated = rsquaredmediated,
                 rsquaredmediated_ci_upper = rsquaredmediated_ci_upper)
  return(results)
  }

  return(rsquaredmediated = rsquaredmediated)
}







upsilons_ols <- function(data, x, m, y,
                         R = 5000,
                         ci_level = 0.95,
                         do_bootstrap = TRUE) {
  # Build formulas
  f_m <- as.formula(paste(m, "~", x))
  f_y <- as.formula(paste(y, "~", x, "+", m))

  # Use a common, complete-case frame for all three variables to keep everything aligned
  vars <- unique(c(x, m, y))
  df <- stats::model.frame(as.formula(paste(vars[1], "~ 1")), data = data, na.action = na.omit)
  if (length(vars) > 1) {
    for (v in vars[-1]) {
      df[[v]] <- data[[v]]
    }
    df <- stats::na.omit(df[, vars, drop = FALSE])
  }

  # Refit models on the aligned data
  mod_m <- stats::lm(f_m, data = df)
  mod_y <- stats::lm(f_y, data = df)

  # Coefficients a (m ~ x) and b (y ~ m)
  a <- stats::coef(mod_m)[x]
  b <- stats::coef(mod_y)[m]

  # Variances of a and b from OLS vcov
  vc_m <- stats::vcov(mod_m)
  vc_y <- stats::vcov(mod_y)

  # Defensive indexing in case of name quirks
  idx_a <- which(names(stats::coef(mod_m)) == x)
  idx_b <- which(names(stats::coef(mod_y)) == m)

  var_a <- vc_m[idx_a, idx_a]
  var_b <- vc_y[idx_b, idx_b]

  # Sample covariance matrix for x and y on the same rows used above
  scov <- stats::cov(df[, c(x, y)], use = "everything")

  # Upsilon and adjusted upsilon
  upsilon <- (a^2 * b^2) * scov[x, x] / scov[y, y]
  adj_upsilon <- ((a^2 - var_a) *
                    (b^2 - var_b) *
                    scov[x, x] / scov[y, y])

  if (do_bootstrap) {
    alpha <- 1 - ci_level
    n <- nrow(df)

    boot_once <- function() {
      idx <- sample.int(n, size = n, replace = TRUE)
      d <- df[idx, , drop = FALSE]

      # Refit on bootstrap sample
      bm <- stats::lm(f_m, data = d)
      by <- stats::lm(f_y, data = d)

      a_b <- suppressWarnings(stats::coef(bm)[x])
      b_b <- suppressWarnings(stats::coef(by)[m])

      # Handle potential failures or singular fits
      if (any(!is.finite(c(a_b, b_b)))) return(c(NA_real_, NA_real_))

      vcm <- try(stats::vcov(bm), silent = TRUE)
      vcy <- try(stats::vcov(by), silent = TRUE)
      if (inherits(vcm, "try-error") || inherits(vcy, "try-error")) return(c(NA_real_, NA_real_))

      ia <- which(names(stats::coef(bm)) == x)
      ib <- which(names(stats::coef(by)) == m)
      if (length(ia) != 1L || length(ib) != 1L) return(c(NA_real_, NA_real_))

      var_a_b <- vcm[ia, ia]
      var_b_b <- vcy[ib, ib]

      sc <- stats::cov(d[, c(x, y)], use = "everything")
      if (!all(is.finite(c(sc[x, x], sc[y, y])))) return(c(NA_real_, NA_real_))

      ups_b <- (a_b^2 * b_b^2) * sc[x, x] / sc[y, y]
      adj_b <- ((a_b^2 - var_a_b) *
                  (b_b^2 - var_b_b) *
                  sc[x, x] / sc[y, y])

      c(ups_b, adj_b)
    }

    # Run bootstrap (matrix with 2 columns: upsilon, adj.upsilon)
    boot_mat <- replicate(R, boot_once())
    boot_mat <- t(boot_mat)

    # Determine which replications converged
    success <- is.finite(boot_mat[, 1]) & is.finite(boot_mat[, 2])
    R_success <- sum(success)
    R_failed  <- R - R_success
    fail_rate <- R_failed / R

    # Compute CIs using only successful reps
    ups_s <- boot_mat[success, 1]
    adj_s <- boot_mat[success, 2]

    ci_ups <- stats::quantile(ups_s, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
    ci_adj <- stats::quantile(adj_s, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)

    return(list(
      upsilon = as.numeric(upsilon),
      adj.upsilon = as.numeric(adj_upsilon),

      ci.upsilon = ci_ups,
      ci.adj.upsilon = ci_adj,

      # ADDED: report bootstrap convergence counts
      boot.converged = R_success,
      boot.failed = R_failed
    ))
  }

  list(
    upsilon = as.numeric(upsilon),
    adj.upsilon = as.numeric(adj_upsilon)
  )
}



















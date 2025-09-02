upsilons_sem_ols <- function(data, x, m, y,
                     R = 5000,
                     ci_level = 0.95,
                     do_bootstrap = TRUE) {
  # Build mediation model (SEM)
  med.model <- paste0(m, " ~ ", x, "\n",
                      y, " ~ ", x, " + ", m)

  # ---------- SEM (lavaan) ----------
  # Fit SEM
  med.res <- lavaan::sem(med.model, data)

  # Path coefficients
  a_sem <- lavaan::coef(med.res)[paste0(m, "~", x)]
  b_sem <- lavaan::coef(med.res)[paste0(y, "~", m)]

  # Variance-covariance of parameters
  vc_sem <- lavaan::vcov(med.res)

  # Sample covariance (SEM's notion of sample cov)
  scov_sem <- lavaan::inspectSampleCov(med.model, data)$cov

  # Upsilon (SEM) and adjusted Upsilon (SEM)
  upsilon_sem <- (a_sem^2 * b_sem^2) * scov_sem[x, x] / scov_sem[y, y]

  adj.upsilon_sem <- ((a_sem^2 - vc_sem[paste0(m, "~", x), paste0(m, "~", x)]) *
                        (b_sem^2 - vc_sem[paste0(y, "~", m), paste0(y, "~", m)]) *
                        scov_sem[x, x] / scov_sem[y, y])

  # ---------- OLS ----------
  # Regressions
  fit_a <- stats::lm(stats::as.formula(paste(m, "~", x)), data = data)
  fit_b <- stats::lm(stats::as.formula(paste(y, "~", x, "+", m)), data = data)

  # Path coefficients (OLS)
  a_ols <- stats::coef(fit_a)[x]
  b_ols <- stats::coef(fit_b)[m]

  # Variance-covariance of coefficients (OLS)
  vc_a_ols <- stats::vcov(fit_a)[x, x]
  vc_b_ols <- stats::vcov(fit_b)[m, m]

  # Plain sample covariance from the data frame
  scov_ols <- stats::cov(data[, c(x, y, m)], use = "pairwise.complete.obs")

  # Upsilon (OLS) and adjusted Upsilon (OLS)
  upsilon_ols <- (a_ols^2 * b_ols^2) * scov_ols[x, x] / scov_ols[y, y]

  adj.upsilon_ols <- ((a_ols^2 - vc_a_ols) *
                        (b_ols^2 - vc_b_ols) *
                        scov_ols[x, x] / scov_ols[y, y])

  # ---------- Bootstrap ----------
  if (do_bootstrap) {
    alpha <- 1 - ci_level

    # --- SEM bootstrap via lavaan ---
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

    message('Bootstrapping (SEM + OLS) may take several minutes…')
    boot.med.res <- lavaan::bootstrapLavaan(med.res, R = R,
                                            type = "ordinary",
                                            FUN = lavaan.med.boot.fun)
    boot.sem.mat <- if (is.list(boot.med.res)) do.call(rbind, boot.med.res) else
      as.matrix(boot.med.res)

    ci.upsilon_sem <- stats::quantile(boot.sem.mat[, 1],
                                      probs = c(alpha/2, 1 - alpha/2),
                                      na.rm = TRUE)
    ci.adj.upsilon_sem <- stats::quantile(boot.sem.mat[, 2],
                                          probs = c(alpha/2, 1 - alpha/2),
                                          na.rm = TRUE)

    # --- OLS bootstrap (pairs bootstrap on rows) ---
    n <- nrow(data)
    boot_ols <- replicate(R, {
      idx <- sample.int(n, n, replace = TRUE)
      db  <- data[idx, , drop = FALSE]

      # Fit bootstrap OLS models; swallow failures to NA
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
          a_val <- coefs[x]
          vc_a  <- vcv[x, x]
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
          b_val <- coefs[m]
          vc_b  <- vcv[m, m]
        }
      }

      scov_b <- try(stats::cov(db[, c(x, y, m)], use = "pairwise.complete.obs"),
                    silent = TRUE)
      if (inherits(scov_b, "try-error")) scov_b <- matrix(NA_real_, 3, 3,
                                                          dimnames = list(c(x, y, m), c(x, y, m)))

      # Compute bootstrap statistics (could be NA if any piece failed)
      ups_b  <- (a_val^2 * b_val^2) * scov_b[x, x] / scov_b[y, y]
      adj_b  <- ((a_val^2 - vc_a) * (b_val^2 - vc_b)) *
        scov_b[x, x] / scov_b[y, y]
      c(ups_b, adj_b)
    })

    boot_ols <- t(boot_ols)
    ci.upsilon_ols <- stats::quantile(boot_ols[, 1],
                                      probs = c(alpha/2, 1 - alpha/2),
                                      na.rm = TRUE)
    ci.adj.upsilon_ols <- stats::quantile(boot_ols[, 2],
                                          probs = c(alpha/2, 1 - alpha/2),
                                          na.rm = TRUE)

    return(list(
      # SEM point estimates + CIs
      upsilon_sem = as.numeric(upsilon_sem),
      adj.upsilon_sem = as.numeric(adj.upsilon_sem),
      ci.upsilon_sem = ci.upsilon_sem,
      ci.adj.upsilon_sem = ci.adj.upsilon_sem,
      # OLS point estimates + CIs
      upsilon_ols = as.numeric(upsilon_ols),
      adj.upsilon_ols = as.numeric(adj.upsilon_ols),
      ci.upsilon_ols = ci.upsilon_ols,
      ci.adj.upsilon_ols = ci.adj.upsilon_ols
    ))
  }

  # If no bootstrap: return point estimates only
  return(list(
    upsilon_sem = as.numeric(upsilon_sem),
    adj.upsilon_sem = as.numeric(adj.upsilon_sem),
    upsilon_ols = as.numeric(upsilon_ols),
    adj.upsilon_ols = as.numeric(adj.upsilon_ols)
  ))
}









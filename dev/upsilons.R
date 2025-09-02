library(lavaan)

upsilons <- function(data, x, m, y,
                     R = 5000,
                     ci_level = 0.95,
                     do_bootstrap = TRUE) {
  # Mediation model using the provided variable names
  med.model <- paste0(m, " ~ ", x, "\n",
                      y, " ~ ", x, " + ", m)

  # Fit SEM
  med.res <- lavaan::sem(med.model, data)

  # Path coefficients
  a <- lavaan::coef(med.res)[paste0(m, "~", x)]
  b <- lavaan::coef(med.res)[paste0(y, "~", m)]

  # Variance-covariance of parameters
  vc <- lavaan::vcov(med.res)

  # Sample covariance matrix (pull the $cov matrix from the list)
  scov <- lavaan::inspectSampleCov(med.model, data)$cov

  # Upsilon
  upsilon <- (a^2 * b^2) * scov[x, x] / scov[y, y]

  # Adjusted upsilon
  adj.upsilon <- ((a^2 - vc[paste0(m, "~", x), paste0(m, "~", x)]) *
                    (b^2 - vc[paste0(y, "~", m), paste0(y, "~", m)]) *
                    scov[x, x] / scov[y, y])


  if (do_bootstrap) {
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

    cat('Bootstrapping may take several minutes \n \n')
    boot.med.res <- lavaan::bootstrapLavaan(med.res, R = R,
                                            type = "ordinary",
                                            FUN = lavaan.med.boot.fun)
    boot.mat <- if (is.list(boot.med.res)) do.call(rbind, boot.med.res) else as.matrix(boot.med.res)

    alpha = 0.05

    # na.rm MUST BE TRUE, or bootstrap runs failed or do not converge
    ci.upsilon     <- quantile(boot.mat[, 1], probs = c(alpha/2, 1 - alpha/2),
                               na.rm = TRUE)
    ci.adj.upsilon <- quantile(boot.mat[, 2], probs = c(alpha/2, 1 - alpha/2),
                               na.rm = TRUE)

    return(list(upsilon = as.numeric(upsilon),
                adj.upsilon = as.numeric(adj.upsilon),
                ci.upsilon = ci.upsilon,
                ci.adj.upsilon = ci.adj.upsilon))
  }

  return(list(upsilon = as.numeric(upsilon),
              adj.upsilon = as.numeric(adj.upsilon)))
}







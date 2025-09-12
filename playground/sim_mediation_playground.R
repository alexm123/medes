# sim_mediation examples


pop_med <- sim_mediation(0.39, 0.59, 0.59, empirical = TRUE)
sim_med <- sim_mediation(tau_prime = 0.39, alpha = 0.59, beta = 0.59,
                         n = 500, empirical = FALSE)
tail(sim_med)
S <- var(sim_med)


upsilons(sim_med, x = "x", m = "m", y = "y", engine="ols")
upsilons_cov(S, x = "x", m = "m", y = "y", n = 500)



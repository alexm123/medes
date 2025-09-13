# --- Parallel setup ---
library(doParallel)
library(foreach)
library(doRNG)

n_workers <- max(1L, parallel::detectCores() - 1L)
cl <- makeCluster(n_workers)
on.exit(stopCluster(cl), add = TRUE)
registerDoParallel(cl)

# Reproducible streams across workers
RNGkind("L'Ecuyer-CMRG")
set.seed(42)
registerDoRNG(42)

sample_sizes <- c(50, 100, 250, 500)
coefficients <- c(0, 0.14, 0.39)

S <- 10  # number of simulated datasets per parameter combo (coverage)
R_mc <- 20 # internal reps for upsilons_mc; increase for stab, decrease for speed

# Build param × replication grid (unchanged)
grid <- params %>% tidyr::crossing(repl = seq_len(S)) %>% as.data.frame()

# --- Parallel run across rows of `grid` ---
results_long <- foreach(
  i = 1:nrow(grid),
  .combine = dplyr::bind_rows,
  .inorder = FALSE,
  .packages = c("tibble","dplyr","tidyr","purrr"),
  .export   = c("sim_one","sim_mediation","upsilons_mc","upsilon_pop","R_mc"),
  .options.snow = list(preschedule = FALSE)   # better load balancing
) %dorng% {
  row <- grid[i, ]
  out <- sim_one(alpha = row$alpha, beta = row$beta,
                 tau_prime = row$tau_prime, n = row$n)

  tibble::tibble(
    alpha  = row$alpha,
    beta   = row$beta,
    tau_prime = row$tau_prime,
    n      = row$n,
    repl   = row$repl,
    adj_upsilon_ols        = as.numeric(out$adj_upsilon_ols),
    ci_adj_upsilon_ols_low = as.numeric(out$ci_adj_upsilon_ols_low),
    ci_adj_upsilon_ols_high= as.numeric(out$ci_adj_upsilon_ols_high),
    upsilon_pop            = upsilon_pop(row$alpha, row$beta)
  )
}

# Coverage flag (same as before)
results_long <- results_long %>%
  mutate(covered_ols = (ci_adj_upsilon_ols_low <= upsilon_pop &
                          upsilon_pop <= ci_adj_upsilon_ols_high))

# --- Summaries (unchanged) ---
coverage_summary <- results_long %>%
  group_by(alpha, beta, tau_prime, n, upsilon_pop) %>%
  summarise(
    mean_adj_upsilon_ols = mean(adj_upsilon_ols),
    coverage_ols         = mean(covered_ols),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 5))) %>%
  arrange(alpha, beta, tau_prime, n)

upsilons_mc(d, x="x", y="y", m="m", do_monte_carlo = TRUE, engine="ols")

upsilons(d, x="x", y="y", m="m", do_bootstrap = TRUE, engine="ols", R = 5000)






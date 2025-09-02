library(dplyr)
library(glue)
library(lavaan)
library(MASS)
library(purrr)

####################################################################
# STANDARDISED SIMULATION FUNCTION
# Fixing residual variances to get variances of 1
####################################################################
# Assumes you already have:
# - d_fake (the "real" data to fit for pop_cov)
# - df_params with columns: N, pop_alpha, pop_beta, pop_tau_prime
# - upsilons_sem_ols() available in the environment
# - packages lavaan, MASS, glue, dplyr loaded
run_simulation_upsilon <- function(pop_tau_prime,
                                   pop_alpha,
                                   pop_beta,
                                   sample_size,
                                   num_reps) {
  # Build the lavaan model string with the supplied population parameters.
  model_string <- glue("
    m ~ {pop_alpha} * x
    y ~ {pop_tau_prime} * x + {pop_beta} * m
    x ~~ 1 * x
    m ~~ {1 - pop_alpha^2} * m
    y ~~ {1 - (pop_tau_prime + pop_beta * pop_alpha)^2 - pop_beta^2 * (1 - pop_alpha^2)} * y
  ")
  d_fake <- data.frame(
    x = rnorm(sample_size),
    m = rnorm(sample_size),
    y = rnorm(sample_size)
  )
  # Fit to d_fake to get the implied covariance matrix for (x, m, y)
  fit <- lavaan::lavaan(model = model_string, data = d_fake)
  pop_cov <- lavaan::lavInspect(fit, "cov.all")
  # Ensure variables are in order (x, m, y) for sampling
  pop_cov <- pop_cov[c("x","m","y"), c("x","m","y")]

  # Do num_reps simulations; each time compute upsilons on a new sample
  sims <- replicate(num_reps, {
    sim_mat <- MASS::mvrnorm(
      n = sample_size,
      mu = rep(0, 3),
      Sigma = pop_cov,
      empirical = FALSE
    )
    sim_df <- as.data.frame(sim_mat)
    names(sim_df) <- c("x","m","y")

    # compute upsilons (no bootstrap)
    u <- upsilons_sem_ols(sim_df, x = "x", m = "m", y = "y", do_bootstrap = FALSE)
    # return a named numeric vector
    c(
      upsilon_sem        = as.numeric(u$upsilon_sem),
      adj_upsilon_sem    = as.numeric(u$adj.upsilon_sem),
      upsilon_ols        = as.numeric(u$upsilon_ols),
      adj_upsilon_ols    = as.numeric(u$adj.upsilon_ols)
    )
  })

  # sims is a matrix: rows = metrics, cols = reps
  sim_means <- rowMeans(sims)

  # Return only the parameter settings and simulation means
  dplyr::tibble(
    pop_alpha     = pop_alpha,
    pop_beta      = pop_beta,
    pop_tau_prime = pop_tau_prime,
    n             = sample_size,
    sim_upsilon_sem     = sim_means["upsilon_sem"],
    sim_adj_upsilon_sem = sim_means["adj_upsilon_sem"],
    sim_upsilon_ols     = sim_means["upsilon_ols"],
    sim_adj_upsilon_ols = sim_means["adj_upsilon_ols"]
  )
}

# Apply to each row of df_params (non-parallel)
sim_function <- function(params_row) {
  run_simulation_upsilon(
    sample_size   = params_row[["N"]],
    pop_alpha     = params_row[["pop_alpha"]],
    pop_beta      = params_row[["pop_beta"]],
    pop_tau_prime = params_row[["pop_tau_prime"]],
    num_reps      = 1000
  )
}

# Run
system_time <- system.time({
  sim_res <- lapply(seq_len(nrow(df_params)), function(i) sim_function(df_params[i, ]))
})

# Bind to a single data frame
df_sim_res <- dplyr::bind_rows(sim_res)


#rounded to 3rd decimal space
df_sim_res_rounded <- df_sim_res %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))

# Add a difference columns for OLS - SEM

df_sim_res_rounded <- df_sim_res_rounded %>%
  dplyr::mutate(diff_adj_ols_minus_sem = sim_adj_upsilon_ols - sim_adj_upsilon_sem)


saveRDS(df_sim_res, file = "df_sim_res_compare_OLS_SEM_rounded.rds")


















res_sem <- upsilons(d_sim, x = "x", m = "m", y = "y", do_bootstrap = FALSE, R = 1000)
res_ols <- upsilons_ols(d_sim, x = "x", m = "m", y = "y", do_bootstrap = TRUE, R = 1000)

res_sem
res_ols

# Making sure simulation and population effect sizes are as they should be
model1_sim <- lm(m ~ x, data = out$sim_matrix)
alpha_sim <- coef(model1_sim)[["x"]]

model2_sim <- lm(y ~ x + m, data = out$sim_matrix)
tau_prime_sim <- coef(model2_sim)[["x"]]

beta_sim <- coef(model2_sim)[["m"]]

c(alpha, beta, tau_prime)


model1_pop <-lm(m ~ x, data = out$pop_data)
alpha_pop <- coef(model1_pop)[["x"]]

model2_pop <- lm(y ~ x + m, data = out$pop_data)
tau_prime_pop <- coef(model2_pop)[["x"]]
beta_pop <- coef(model2_pop)[["m"]]

c(alpha_pop, beta_pop, tau_prime_pop)



# Making parameters
sample_sizes <- c(20, 50, 100, 200, 500, 1000)
effect_sizes <- c(0.00, 0.14, 0.39, 0.59)

# Containing all the conditions
df_params <- expand.grid(
  N = sample_sizes,
  pop_alpha = effect_sizes,
  pop_beta = effect_sizes,
  pop_tau_prime = effect_sizes
)

# Exclude the non–positive-definite combo (0.59 for all three)
df_params <- subset(
  df_params,
  !(pop_alpha == 0.59 & pop_beta == 0.59 & pop_tau_prime == 0.59)
)

simulation_function <- function(params) {
  run_simulation_small_standardized(
    sample_size  = params[["N"]],
    pop_alpha    = params[["pop_alpha"]],
    pop_beta     = params[["pop_beta"]],
    pop_tau_prime= params[["pop_tau_prime"]],
    num_reps     = 20
  )
}

# Run sequentially (no parallelization)
system_time <- system.time({
  sim_res <- lapply(1:nrow(df_params), function(i) simulation_function(df_params[i, ]))
})

head(sim_res)[[1]]














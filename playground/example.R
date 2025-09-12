devtools::install_github("alexm123/medes")
library(medes)
library(mice)


#Example data for mediation analysis
set.seed(41)
n <- 500

# Independent variable (x)
x <- rnorm(n, mean = 50, sd = 10)

# Mediator (m), partially dependent on x
m  <- 0.5*x + rnorm(n, mean = 0, sd = 5)

# Dependent variable (y), dependent on x and m
y <- 0.3*x + 0.6*m + rnorm(n, mean = 0, sd = 5)

# Create dataframe
d <- data.frame(x, m, y)

d_scaled <- scale(d) |> as.data.frame()

# Introduce MAR 30%
amp <- ampute(
  data = d,
  prop = 0.30,          # overall missing fraction
  mech = "MAR"          # MAR (not MCAR/MNAR)
)

d_mar <- amp$amp


# 5e6 = 5 million
pop_n <- 5e4
pop_x <- rnorm(pop_n, mean = 50, sd = 10)
pop_m <- 0.5*pop_x + rnorm(pop_n, mean = 0, sd = 5)
pop_y <- 0.3*pop_x + 0.6*pop_m + rnorm(pop_n, mean = 0, sd = 5)
pop_d <- data.frame(pop_x, pop_m, pop_y)

# MAR

# View the first few rows
head(d)

lm_test_pop <- lm(data = pop_d, pop_y ~ pop_x)
summary(lm_test)
rm(lm_test_pop)

lm_test <- lm(data = d, y ~ x)
summary(lm_test)
rm(lm_test)

# delete population when done to save memory
rm(pop_n, pop_x, pop_m, pop_y, pop_d)

# Z-score the variables
Z <- as.data.frame(scale(d))

# a: standardized slope from X -> M (equals cor(x, m))
a  <- coef(lm(m ~ x, data = Z))["x"]

# b and c': standardized partial slopes from Y ~ M + X
coef_y <- coef(lm(y ~ m + x, data = Z))
b      <- unname(coef_y["m"])   # M -> Y | X
cprime <- unname(coef_y["x"])   # X -> Y | M

c(a = a, b = b, cprime = cprime)

###############################################################################
###############################################################################
###############################################################################


rsquare_med(d, x = "x", m = "m", y = "y",
                 do_bootstrap = TRUE, R = 100)



rsq_med_TRUE <- rsquare_med_test(pop_d, x = "pop_x", m = "pop_m", y = "pop_y",
                                do_bootstrap = TRUE, engine='both', R = 100)

# Rsquare mediated with confidence intervals
rsq_med_medes <- rsquare_med(d, x = "x", m = "m", y = "y",
                 do_bootstrap = TRUE, R = 100, engine='both')

# Rsquare mediated with confidence intervals on MAR 30% data
rsq_med_medes_MAR <- rsquare_med(d_mar, x = "x", m = "m", y = "y",
                                 do_bootstrap = TRUE, R = 100, engine='sem')


rsq_med_components_test <- rsquare_med_components(d, x = "x", m = "m", y = "y",
                                                  do_bootstrap = FALSE, engine='both',
                                                  R = 100)
rsq_med_components_test_scaled <- rsquare_med_components(d_scaled, x = "x", m = "m", y = "y",
                                            do_bootstrap = FALSE, engine='both',
                                            R = 100)

rsquare_med_test(d_scaled, x = "x", m = "m", y = "y",
                  do_bootstrap = TRUE, R = 100, engine='both')




ups_medes <- upsilons(d, x = "x", m = "m", y = "y",
                      do_bootstrap = TRUE,
                      R = 100, engine='both', stat='both')


res_mbess <- MBESS::mediation(x = d[[1]], mediator = d[[2]], dv = d[[3]],
                 bootstrap = TRUE, B = 100, which.boot="percentile")

# everything <- medes_all(d, x = "x", m = "m", y = "y",
#                         do_bootstrap = FALSE)



# rsq_indirect(d, x = "x", m = "m", y = "y", do_bootstrap = TRUE)
res_sem <- upsilons(d, x = "x", m = "m", y = "y", do_bootstrap = TRUE, R = 100)
res_ols <- upsilons_ols(d, x = "x", m = "m", y = "y", do_bootstrap = TRUE, R = 1000)
res_sem
res_ols

# Generate data!
######################
######################
######################
sim_dat <- sim_mediation(0.4, 0.4, 0.4, 1000)
var(sim_dat)
MBESS::upsilon(x = sim_dat[[1]], mediator = sim_dat[[2]], dv = sim_dat[[3]], bootstrap = FALSE)

upsilons(sim_dat, x = "x", m = "m", y = "y",
         do_bootstrap = FALSE, stat='both')

model1 <- lm(m ~ x, data = sim_dat)
model2 <- lm(y ~ x + m, data = sim_dat)

QuantPsyc::lm.beta(model1)
QuantPsyc::lm.beta(model2)


######################
######################
######################




rsq_indirect(d, "x", "m", "y", R = 1000)


res_both <- upsilons_sem_ols(d, x = "x", m = "m", y = "y", do_bootstrap = FALSE)





MBESS::upsilon(x = d[[1]], mediator = d[[2]], dv = d[[3]], B = 100)




MBESS::mediation(x = d[[1]], mediator = d[[2]], dv = d[[3]])$Effect.Sizes


set.seed(42)
needed <- 5
max_attempts <- 50

res <- vector("list", needed)
got <- 0
attempt <- 0

while (got < needed && attempt < max_attempts) {
  attempt <- attempt + 1
  out <- tryCatch(
    MBESS::upsilon(x = d[[1]], mediator = d[[2]], dv = d[[3]], B = 1000),
    error = function(e) NULL
  )
  # Coerce to data.frame and validate structure/values
  if (!is.null(out)) {
    tbl <- tryCatch(as.data.frame(out), error = function(e) NULL)
    ok_names <- !is.null(tbl) &&
      all(c("Upsilon","Adj Upsilon") %in% rownames(tbl)) &&
      all(c("Estimate","95% ordinary LCL","95% ordinary UCL") %in% colnames(tbl))
    if (ok_names) {
      vals <- unlist(tbl[c("Upsilon","Adj Upsilon"),
                         c("Estimate","95% ordinary LCL","95% ordinary UCL")])
      if (all(is.finite(vals))) {
        got <- got + 1
        res[[got]] <- tbl
      }
    }
  }
}

cat(sprintf("%d successful runs out of %d attempts\n", got, attempt))

# Use only successful runs
valid_res <- res[seq_len(got)]

# Stack the successful tables and average by row/column
stacked <- lapply(
  valid_res,
  function(tbl) tbl[c("Upsilon","Adj Upsilon"),
                    c("Estimate","95% ordinary LCL","95% ordinary UCL")]
)

# helper: average a row across all runs
row_avg <- function(row_name) {
  mat <- do.call(rbind, lapply(stacked, function(t) as.numeric(t[row_name, ])))
  col_means <- colMeans(mat, na.rm = TRUE)
  names(col_means) <- c("Estimate","95% ordinary LCL","95% ordinary UCL")
  col_means
}

# helper: range (lowest LCL, highest UCL) across runs for a row
row_range <- function(row_name) {
  LCLs <- vapply(valid_res, function(t) as.numeric(t[row_name, "95% ordinary LCL"]), numeric(1))
  UCLs <- vapply(valid_res, function(t) as.numeric(t[row_name, "95% ordinary UCL"]), numeric(1))
  c(`min LCL` = min(LCLs, na.rm = TRUE),
    `max UCL` = max(UCLs, na.rm = TRUE))
}

avg_bootstraps <- list(
  Upsilon = c(row_avg("Upsilon"),     row_range("Upsilon")),
  Adj_Upsilon = c(row_avg("Adj Upsilon"), row_range("Adj Upsilon")),
  runs_used = got,
  attempts  = attempt
)

avg_bootstraps

###############################################################################
###############################################################################
###############################################################################







# Function to list top memory consumers in the environment
who_big <- function(env = .GlobalEnv, n = 10) {
  objs <- ls(envir = env)
  sizes <- sapply(objs, function(x) object.size(get(x, envir = env)))
  sizes <- sort(sizes, decreasing = TRUE)

  # return top n as data.frame
  data.frame(
    Object = names(sizes),
    Size_MB = round(sizes / 1024^2, 2)
  )[1:min(n, length(sizes)), ]
}

# Example usage:
who_big()        # top 10
who_big(n = 20)  # top 20



# Toy dataset from "Friendship" Zhang et al...
# --- Published moments -------------------------------------------------------
library(lavaan)
vars <- c("loneliness","social_comp","friendship_quality","social_pref","prox_prestige")

# Correlations from the table (upper/lower triangle filled symmetrically)
R <- matrix(c(
  1.00,  -0.69, -0.41, -0.44, -0.37,
  -0.69,   1.00,  0.34,  0.29,  0.23,
  -0.41,   0.34,  1.00,  0.31,  0.26,
  -0.44,   0.29,  0.31,  1.00,  0.60,
  -0.37,   0.23,  0.26,  0.60,  1.00
), 5, 5, byrow = TRUE, dimnames = list(vars, vars))

mu  <- c(1.73, 2.85, 2.99, 0.00, 0.42)        # means
sdv <- c(0.72, 0.53, 0.75, 1.72, 0.10)        # SDs
N   <- 509

# Covariance matrix from correlations + SDs
S <- diag(sdv) %*% R %*% diag(sdv)
dimnames(S) <- list(vars, vars)
# --- Model mapping -----------------------------------------------------------
Y  <- "loneliness"
M  <- "social_pref"
X1 <- "social_comp"
X2 <- "friendship_quality"
X3 <- "prox_prestige"   # last predictor; we'll fix X3 -> M to zero

# --- Constrained model: X3 -> M fixed to 0 ----------------------------------
model_constrained <- glue::glue('
  {M} ~ a1*{X1} + a2*{X2}
  {Y} ~ c1p*{X1} + c2p*{X2} + c3p*{X3} + b*{M}

  # allow predictors to correlate
  {X1} ~~ {X2} + {X3}
  {X2} ~~ {X3}

  # effects
  ab1 := a1*b
  ab2 := a2*b
  total1 := c1p + ab1
  total2 := c2p + ab2
')

fitC <- sem(model_constrained,
            sample.cov = S,
            sample.mean = mu,
            sample.nobs = N,
            meanstructure = TRUE,
            fixed.x = FALSE)

summary(fitC, standardized = TRUE, rsquare = TRUE)
#parameterEstimates(fitC, standardized = TRUE)[,c("lhs","op","rhs","label","est","se","z","pvalue","std.all")]
modindices(fitC)
# --- Unconstrained model (lets X3 -> M be estimated; should be ~0 if truly absent)
model_free <- glue::glue('
  {M} ~ a1*{X1} + a2*{X2} + a3*{X3}
  {Y} ~ c1p*{X1} + c2p*{X2} + c3p*{X3} + b*{M}

  {X1} ~~ {X2} + {X3}
  {X2} ~~ {X3}

  ab1 := a1*b
  ab2 := a2*b
  ab3 := a3*b
  total1 := c1p + ab1
  total2 := c2p + ab2
  total3 := c3p + ab3
')

fitF <- sem(model_free,
            sample.cov = S,
            sample.mean = mu,
            sample.nobs = N,
            meanstructure = TRUE,
            fixed.x = FALSE)

summary(fitF, standardized = TRUE, rsquare = TRUE)
anova(fitF, fitC)
set.seed(42)
dat <- as.data.frame(mvrnorm(n = N, mu = mu, Sigma = S, empirical = TRUE))
names(dat) <- vars
var(dat)

# Effect sizes with toy dataset!
rsquare_med_toy <- rsquare_med(dat, x = "social_pref",
                               m = "social_comp",
                               y = "loneliness",
                               do_bootstrap = TRUE, R = 100, engine='both')
ups_medes_toy <- upsilons(dat, x = "friendship_quality",
                          m = "social_comp",
                          y = "loneliness",
                          do_bootstrap = TRUE, R = 100, engine='both')






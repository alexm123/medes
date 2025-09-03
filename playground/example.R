# Example data for mediation analysis
set.seed(42)
n <- 2000

# Independent variable (x)
x <- rnorm(n, mean = 50, sd = 10)

# Mediator (m), partially dependent on x
m  <- 0.5*x + rnorm(n, mean = 0, sd = 5)

# Dependent variable (y), dependent on x and m
y <- 0.3*x + 0.6*m + rnorm(n, mean = 0, sd = 5)

# Create dataframe
d <- data.frame(x, m, y)

# View the first few rows
head(d)


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
#
# rsquare_med(d, x = "x", m = "m", y = "y", do_bootstrap = TRUE)
# rsq_indirect(d, x = "x", m = "m", y = "y", do_bootstrap = TRUE)
res_sem <- upsilons(d, x = "x", m = "m", y = "y", do_bootstrap = TRUE, R = 100)
res_ols <- upsilons_ols(d, x = "x", m = "m", y = "y", do_bootstrap = TRUE, R = 100)

upsilons_sem_ols(d, x = "x", m = "m", y = "y", do_bootstrap = TRUE, R = 100)

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



rsq_med_from_paths <- function(a, b, cprime) {
  a^2 * (b^2 + cprime^2) + 2*a*b*cprime
}
rsq_ind_from_paths <- function(a, b) (1 - a^2) * (a * b)^2


rsq_med_from_paths(a =      0.39,
                   b =      0.39,
                   cprime = 0.59)


rsq_ind_from_paths(a =      0.39,
                   b =      0.39)


(0.39^2) * (0.39^2)








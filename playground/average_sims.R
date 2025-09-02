set.seed(42)
needed <- 5
max_attempts <- 50

## -------------------------
## Helper functions
## -------------------------
row_avg <- function(stacked, row_name) {
  mat <- do.call(rbind, lapply(stacked, function(t) as.numeric(t[row_name, ])))
  col_means <- colMeans(mat, na.rm = TRUE)
  names(col_means) <- c("Estimate","95% ordinary LCL","95% ordinary UCL")
  col_means
}

row_range <- function(valid_res, row_name) {
  LCLs <- vapply(valid_res, function(t) as.numeric(t[row_name, "95% ordinary LCL"]), numeric(1))
  UCLs <- vapply(valid_res, function(t) as.numeric(t[row_name, "95% ordinary UCL"]), numeric(1))
  c(`min LCL` = min(LCLs, na.rm = TRUE),
    `max UCL` = max(UCLs, na.rm = TRUE))
}

## -------------------------
## MBESS::upsilon runs
## -------------------------
res_mbess <- vector("list", needed)
got_mbess <- 0; attempt_mbess <- 0

while (got_mbess < needed && attempt_mbess < max_attempts) {
  attempt_mbess <- attempt_mbess + 1
  out <- tryCatch(
    MBESS::upsilon(x = d[[1]], mediator = d[[2]], dv = d[[3]], B = 1000),
    error = function(e) NULL
  )
  if (!is.null(out)) {
    tbl <- tryCatch(as.data.frame(out), error = function(e) NULL)
    ok_names <- !is.null(tbl) &&
      all(c("Upsilon","Adj Upsilon") %in% rownames(tbl)) &&
      all(c("Estimate","95% ordinary LCL","95% ordinary UCL") %in% colnames(tbl))
    if (ok_names) {
      vals <- unlist(tbl[c("Upsilon","Adj Upsilon"),
                         c("Estimate","95% ordinary LCL","95% ordinary UCL")])
      if (all(is.finite(vals))) {
        got_mbess <- got_mbess + 1
        res_mbess[[got_mbess]] <- tbl
      }
    }
  }
}

valid_mbess <- res_mbess[seq_len(got_mbess)]
stacked_mbess <- lapply(valid_mbess, function(tbl) tbl[c("Upsilon","Adj Upsilon"),
                                                       c("Estimate","95% ordinary LCL","95% ordinary UCL")])

## -------------------------
## upsilons() runs
## -------------------------
res_simple <- vector("list", needed)
got_simple <- 0; attempt_simple <- 0

while (got_simple < needed && attempt_simple < max_attempts) {
  attempt_simple <- attempt_simple + 1
  out <- tryCatch(
    upsilons(d, x = "x", m = "m", y = "y", do_bootstrap = TRUE, R = 5000),
    error = function(e) NULL
  )
  if (!is.null(out) && is.finite(out$upsilon) && is.finite(out$adj.upsilon)) {
    got_simple <- got_simple + 1
    res_simple[[got_simple]] <- out
  }
}

ups_vals     <- vapply(res_simple[seq_len(got_simple)], function(r) r$upsilon, numeric(1))
adj_ups_vals <- vapply(res_simple[seq_len(got_simple)], function(r) r$adj.upsilon, numeric(1))

ci_u_mat  <- do.call(cbind, lapply(res_simple[seq_len(got_simple)], function(r) r$ci.upsilon))
ci_au_mat <- do.call(cbind, lapply(res_simple[seq_len(got_simple)], function(r) r$ci.adj.upsilon))

## -------------------------
## Combined summary
## -------------------------
avg_bootstraps <- list(
  MBESS = list(
    Upsilon     = c(row_avg(stacked_mbess,"Upsilon"),     row_range(valid_mbess,"Upsilon")),
    Adj_Upsilon = c(row_avg(stacked_mbess,"Adj Upsilon"), row_range(valid_mbess,"Adj Upsilon")),
    runs_used   = got_mbess,
    attempts    = attempt_mbess
  ),
  upsilons = list(
    avg_upsilon        = mean(ups_vals,     na.rm = TRUE),
    avg_adj_upsilon    = mean(adj_ups_vals, na.rm = TRUE),
    avg_ci_upsilon     = rowMeans(ci_u_mat,  na.rm = TRUE),
    avg_ci_adj_upsilon = rowMeans(ci_au_mat, na.rm = TRUE),
    range_ci_upsilon   = c(min(ci_u_mat[1, ], na.rm = TRUE), max(ci_u_mat[2, ], na.rm = TRUE)),
    range_ci_adj       = c(min(ci_au_mat[1, ], na.rm = TRUE), max(ci_au_mat[2, ], na.rm = TRUE)),
    runs_used          = got_simple,
    attempts           = attempt_simple
  )
)

avg_bootstraps

# R/edm_simplex.R
# ------------------------------------------------------------------------------
# EDM: Simplex forecasting + skill metrics
#
# Purpose
# - Provide a minimal, dependency-light simplex forecaster:
#   - delay embedding in R^E
#   - E+1 nearest neighbors
#   - exponential distance weighting (local averaging)
#
# Contract
# - make_embed(x, E, tau) -> matrix (n_embed x E)
# - simplex_forecast(x, E, tau, Tp) -> list(obs, pred)
# - skill_metrics(obs, pred) -> list(rho, rmse)
# ------------------------------------------------------------------------------

make_embed <- function(x, E, tau = 1) {
  if (!is.numeric(x)) stop("make_embed: x must be numeric.")
  if (any(!is.finite(x))) stop("make_embed: x must be finite (no NA/Inf).")
  if (!is.numeric(E) || E < 2) stop("make_embed: E must be >= 2.")
  if (!is.numeric(tau) || tau < 1) stop("make_embed: tau must be >= 1.")

  n <- length(x)
  m <- n - (E - 1) * tau
  if (m <= 1) stop("make_embed: time series too short for requested embedding.")

  M <- matrix(NA_real_, nrow = m, ncol = E)
  for (j in 1:E) {
    idx_start <- 1 + (j - 1) * tau
    idx_end <- idx_start + m - 1
    M[, j] <- x[idx_start:idx_end]
  }
  M
}

simplex_forecast <- function(x, E, tau = 1, Tp = 1) {
  if (!is.numeric(Tp) || Tp < 1) stop("simplex_forecast: Tp must be >= 1.")
  M <- make_embed(x, E, tau)
  nM <- nrow(M)

  target_idx <- (1:nM) + (E - 1) * tau + Tp
  valid <- which(target_idx <= length(x))
  if (length(valid) < (E + 2)) stop("simplex_forecast: not enough points for forecast.")

  Mv <- M[valid, , drop = FALSE]
  y <- x[target_idx[valid]]
  preds <- rep(NA_real_, length(y))

  for (i in seq_along(y)) {
    d <- sqrt(rowSums((Mv - Mv[i, ])^2))
    d[i] <- Inf

    nn <- order(d)[1:(E + 1)]
    d_nn <- d[nn]
    d1 <- d_nn[1]

    # Standard simplex weighting: closer neighbors dominate exponentially.
    w <- if (d1 == 0) rep(1, length(d_nn)) else exp(-d_nn / d1)
    preds[i] <- sum(w * y[nn]) / sum(w)
  }

  list(obs = y, pred = preds)
}

skill_metrics <- function(obs, pred) {
  if (!is.numeric(obs) || !is.numeric(pred)) stop("skill_metrics: obs/pred must be numeric.")
  if (length(obs) != length(pred)) stop("skill_metrics: obs and pred must have same length.")

  ok <- is.finite(obs) & is.finite(pred)
  if (!any(ok)) return(list(rho = NA_real_, rmse = NA_real_))

  o <- obs[ok]
  p <- pred[ok]

  rho <- suppressWarnings(cor(o, p))
  rmse <- sqrt(mean((o - p)^2))

  list(rho = rho, rmse = rmse)
}

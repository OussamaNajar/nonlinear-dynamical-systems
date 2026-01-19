# R/windowed_smap.R
# ------------------------------------------------------------------------------
# Windowed EDM: S-map diagnostics over time
#
# Purpose
# - Quantify state-dependence and local linearization changes by measuring:
#   (1) theta_opt(t): theta maximizing one-step out-of-sample forecast skill
#   (2) delta_rho(t): improvement over theta=0 (global linear baseline)
#   (3) b(t) summaries at fixed theta (e.g., theta=8): median slope and slope SD
#
# Contract
# - windowed_smap_diagnostics(x, thetas, E, tau, Tp, window, step, theiler, lib_frac, theta_slope)
#   returns a data.frame with:
#     center_idx, theta_opt, rho0, rho_opt, delta_rho, b_median, b_sd
#
# Assumptions
# - x is finite numeric (no NA/Inf)
# - make_embed(), skill_metrics() are available (from R/edm_simplex.R)
# ------------------------------------------------------------------------------

# internal helper: S-map predictions + slopes for ONE window at a fixed theta
smap_window_fit <- function(seg, E, tau, Tp, theta, theiler, lib_frac, ridge = 1e-8) {
  M <- make_embed(seg, E, tau)
  nM <- nrow(M)
  if (nM < (E + 5)) {
    return(list(obs = numeric(0), pred = numeric(0), slope = numeric(0)))
  }

  target_idx <- (1:nM) + (E - 1) * tau + Tp
  valid <- which(target_idx <= length(seg))
  Mv <- M[valid, , drop = FALSE]
  y  <- seg[target_idx[valid]]
  nV <- nrow(Mv)
  if (nV < (E + 10)) {
    return(list(obs = numeric(0), pred = numeric(0), slope = numeric(0)))
  }

  nLib <- max(10, floor(lib_frac * nV))
  lib_ix <- 1:nLib
  pred_ix <- (nLib + 1):nV
  if (length(pred_ix) < 5) {
    return(list(obs = numeric(0), pred = numeric(0), slope = numeric(0)))
  }

  # Design matrix for local linear model: [1, state...]
  Xlib <- cbind(1, Mv[lib_ix, , drop = FALSE])
  ylib <- y[lib_ix]

  preds  <- rep(NA_real_, length(pred_ix))
  slopes <- rep(NA_real_, length(pred_ix))

  for (j in seq_along(pred_ix)) {
    i <- pred_ix[j]

    # distances from prediction point to library points
    d <- sqrt(rowSums((Mv[lib_ix, , drop = FALSE] - Mv[i, ])^2))

    # Theiler exclusion in window time-index coordinates
    time_dist <- abs(lib_ix - i)
    d[time_dist <= theiler] <- Inf

    ord <- order(d)
    ord <- ord[is.finite(d[ord])]
    if (length(ord) < (E + 1)) next

    # localization scale: mean finite distance (robust)
    dbar <- mean(d[ord], na.rm = TRUE)
    if (!is.finite(dbar) || dbar <= 0) next

    w <- exp(-theta * d / dbar)
    w[!is.finite(w)] <- 0
    if (sum(w) <= 0) next

    # Weighted least squares: beta = (X'WX)^-1 X'W y
    W <- w
    XtW <- t(Xlib) * W
    G <- XtW %*% Xlib
    rhs <- XtW %*% ylib

    # Numerical stabilization: tiny ridge on the diagonal to reduce blow-ups when
    # the local design matrix is near-singular (common in strongly autocorrelated embeddings).
    if (!is.numeric(ridge) || ridge < 0) ridge <- 0
    if (ridge > 0) {
      p <- nrow(G)
      lam <- ridge * (sum(diag(G)) / max(1, p))
      G <- G + diag(lam, p)
    }

    beta <- tryCatch(solve(G, rhs), error = function(e) NULL)
    if (is.null(beta) || any(!is.finite(beta))) next

    x_i <- c(1, as.numeric(Mv[i, ]))
    preds[j] <- sum(x_i * beta)

    # Slope proxy: coefficient on the first coordinate (current Z_t)
    slopes[j] <- beta[2]
  }

  list(obs = y[pred_ix], pred = preds, slope = slopes)
}

windowed_smap_diagnostics <- function(x,
                                     thetas = seq(0, 8, by = 0.5),
                                     E, tau = 1, Tp = 1,
                                     window = 800, step = 50,
                                     theiler = 50, lib_frac = 0.8,
                                     theta_slope = 8,
                                     ridge = 1e-8) {
  if (!is.numeric(x) || any(!is.finite(x))) stop("windowed_smap_diagnostics: x must be finite numeric.")
  if (!is.numeric(thetas) || length(thetas) < 2) stop("windowed_smap_diagnostics: invalid thetas.")
  if (!any(abs(thetas - 0) < 1e-12)) stop("windowed_smap_diagnostics: thetas must include 0.")
  if (!exists("make_embed", mode = "function")) stop("windowed_smap_diagnostics: make_embed() not found.")
  if (!exists("skill_metrics", mode = "function")) stop("windowed_smap_diagnostics: skill_metrics() not found.")
  if (!(lib_frac > 0 && lib_frac < 1)) stop("windowed_smap_diagnostics: lib_frac must be in (0,1).")

  n <- length(x)
  if (window < 50 || window >= n) stop("windowed_smap_diagnostics: invalid window.")
  if (step < 1) stop("windowed_smap_diagnostics: step must be >= 1.")

  starts <- seq(1, n - window + 1, by = step)

  out <- data.frame(
    center_idx = integer(0),
    theta_opt  = numeric(0),
    rho0       = numeric(0),
    rho_opt    = numeric(0),
    delta_rho  = numeric(0),
    b_median   = numeric(0),
    b_sd       = numeric(0)
  )

  for (s in starts) {
    e <- s + window - 1
    seg <- x[s:e]

    # baseline theta=0
    fit0 <- smap_window_fit(seg, E, tau, Tp, theta = 0, theiler = theiler, lib_frac = lib_frac, ridge = ridge)
    met0 <- skill_metrics(fit0$obs, fit0$pred)
    rho0 <- met0$rho

    # sweep theta for best rho
    best_theta <- NA_real_
    best_rho <- -Inf

    for (th in thetas) {
      fit <- smap_window_fit(seg, E, tau, Tp, theta = th, theiler = theiler, lib_frac = lib_frac, ridge = ridge)
      met <- skill_metrics(fit$obs, fit$pred)
      if (is.finite(met$rho) && met$rho > best_rho) {
        best_rho <- met$rho
        best_theta <- th
      }
    }

    # slope diagnostics at fixed theta
    fits <- smap_window_fit(seg, E, tau, Tp, theta = theta_slope, theiler = theiler, lib_frac = lib_frac, ridge = ridge)
    b <- fits$slope
    b <- b[is.finite(b)]
    b_med <- if (length(b) == 0) NA_real_ else median(b)
    b_sd  <- if (length(b) < 2) NA_real_ else sd(b)

    out <- rbind(out, data.frame(
      center_idx = as.integer(round((s + e) / 2)),
      theta_opt  = best_theta,
      rho0       = rho0,
      rho_opt    = best_rho,
      delta_rho  = best_rho - rho0,
      b_median   = b_med,
      b_sd       = b_sd
    ))
  }

  out
}

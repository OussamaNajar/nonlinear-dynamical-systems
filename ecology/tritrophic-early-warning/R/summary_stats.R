# R/summary_stats.R
# ==============================================================================
# Early-warning quantification + compact pipeline summary (evidence-first)
#
# Design goals
# - Metrics are explicit about assumptions and failure modes.
# - No hard-coded column names beyond defaults (Tp selectable).
# - Warning-time thresholds are derived from baseline distributions where possible.
# - Output is machine-readable first (CSV); console printing is optional/minimal.
# ==============================================================================

# ---- helpers ---------------------------------------------------------------

require_cols <- function(df, cols, where = "data.frame") {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(sprintf("%s missing required columns: %s", where, paste(missing, collapse = ", ")))
  }
  invisible(TRUE)
}

kendall_trend <- function(t, x) {
  ok <- is.finite(t) & is.finite(x)
  if (sum(ok) < 5) return(NA_real_)
  suppressWarnings(cor(t[ok], x[ok], method = "kendall"))
}

first_crossing_time <- function(time, x, thr, t_star) {
  if (!is.finite(thr) || !is.finite(t_star)) return(NA_real_)
  ok <- is.finite(time) & is.finite(x) & (time < t_star)
  idx <- which(ok & (x > thr))
  if (length(idx) == 0) return(NA_real_)
  time[min(idx)]
}

baseline_quantile <- function(time, x, t_star, lookback_window, q = 0.99) {
  # Baseline region is "well before" the pre-transition lookback window.
  if (!is.finite(t_star)) return(NA_real_)
  ok <- is.finite(time) & is.finite(x) & (time < (t_star - lookback_window))
  if (sum(ok) < 20) return(NA_real_)
  as.numeric(stats::quantile(x[ok], probs = q, na.rm = TRUE, names = FALSE))
}

# ---- main metrics -----------------------------------------------------------

compute_early_warning_metrics <- function(df_ws,
                                         df_wm,
                                         det,
                                         t_star,
                                         lookback_window = 400,
                                         Tp_for_rho = 1,
                                         theta_warn_method = c("baseline_q99", "fixed"),
                                         theta_fixed = 2,
                                         delta_q = 0.99) {

  theta_warn_method <- match.arg(theta_warn_method)

  # Minimal structural requirements
  require_cols(df_ws, c("time"), "df_ws")
  require_cols(df_wm, c("time", "theta_opt", "delta_rho", "b_median", "b_sd"), "df_wm")

  rho_col <- paste0("rho_Tp", Tp_for_rho)
  if (!(rho_col %in% names(df_ws))) {
    stop(sprintf("compute_early_warning_metrics: df_ws missing %s; choose Tp_for_rho that exists.", rho_col))
  }

  if (!is.finite(t_star)) {
    warning("compute_early_warning_metrics: t_star not finite; returning NA metrics.")
    return(list(
      kendall_rho = NA_real_,
      kendall_theta = NA_real_,
      kendall_delta_rho = NA_real_,
      kendall_b_median = NA_real_,
      kendall_b_sd = NA_real_,
      mean_advance_time = NA_real_,
      thresholds = list(theta = NA_real_, delta_rho = NA_real_)
    ))
  }

  # Pre-transition analysis window (where early warning is expected)
  t0 <- t_star - lookback_window

  # --- Simplex skill trend (expect decline) ---
  idx_rho <- is.finite(df_ws$time) & df_ws$time >= t0 & df_ws$time < t_star
  tau_rho <- kendall_trend(df_ws$time[idx_rho], df_ws[[rho_col]][idx_rho])

  # --- S-map diagnostics trends (expect increase in nonlinearity and Δρ) ---
  idx_s <- is.finite(df_wm$time) & df_wm$time >= t0 & df_wm$time < t_star
  tau_theta <- kendall_trend(df_wm$time[idx_s], df_wm$theta_opt[idx_s])
  tau_delta <- kendall_trend(df_wm$time[idx_s], df_wm$delta_rho[idx_s])
  tau_b_med <- kendall_trend(df_wm$time[idx_s], df_wm$b_median[idx_s])
  tau_b_sd  <- kendall_trend(df_wm$time[idx_s], df_wm$b_sd[idx_s])

  # --- Warning-time thresholds (derived, not vibes) ---
  # theta threshold:
  #   - baseline_q99: theta exceeds extreme baseline quantile (robust to scale)
  #   - fixed: legacy heuristic (keep available, but make it explicit)
  theta_thr <- if (theta_warn_method == "baseline_q99") {
    baseline_quantile(df_wm$time, df_wm$theta_opt, t_star, lookback_window, q = 0.99)
  } else {
    as.numeric(theta_fixed)
  }

  # delta_rho threshold: extreme baseline quantile
  delta_thr <- baseline_quantile(df_wm$time, df_wm$delta_rho, t_star, lookback_window, q = delta_q)

  t_theta_warn <- first_crossing_time(df_wm$time, df_wm$theta_opt, theta_thr, t_star)
  t_delta_warn <- first_crossing_time(df_wm$time, df_wm$delta_rho, delta_thr, t_star)

  warn_times <- c(t_theta_warn, t_delta_warn)
  warn_times <- warn_times[is.finite(warn_times)]
  mean_advance <- if (length(warn_times) > 0) t_star - mean(warn_times) else NA_real_

  list(
    kendall_rho = tau_rho,
    kendall_theta = tau_theta,
    kendall_delta_rho = tau_delta,
    kendall_b_median = tau_b_med,
    kendall_b_sd = tau_b_sd,
    mean_advance_time = mean_advance,
    thresholds = list(theta = theta_thr, delta_rho = delta_thr),
    warning_times = list(theta = t_theta_warn, delta_rho = t_delta_warn),
    window = list(t0 = t0, t_star = t_star, lookback = lookback_window),
    rho_channel = rho_col
  )
}

# ---- summary table ----------------------------------------------------------

create_summary_table <- function(config, det, global_metrics, ew_metrics, path) {

  # Safe fetch helpers
  gm <- global_metrics
  gm_simplex1  <- if (!is.null(gm$simplex_Tp1))  gm$simplex_Tp1$rho  else NA_real_
  gm_simplex10 <- if (!is.null(gm$simplex_Tp10)) gm$simplex_Tp10$rho else NA_real_

  out <- data.frame(
    metric = c(
      "sim_t1",
      "sim_dt_assumed_1",
      "drift_start_time",
      "c2_start",
      "c2_end",
      "transition_t_star",
      "rollsd_window_w",
      "threshold_k",
      "simplex_E",
      "simplex_tau",
      "simplex_rho_h1",
      "simplex_rho_h10",
      "smap_E",
      "smap_theta_opt",
      "smap_rho_theta0",
      "smap_rho_thetaopt",
      "windowed_window",
      "windowed_step",
      "windowed_theiler",
      "ew_kendall_rho",
      "ew_kendall_theta",
      "ew_kendall_delta_rho",
      "ew_kendall_b_median",
      "ew_kendall_b_sd",
      "ew_mean_advance_time",
      "ew_theta_threshold",
      "ew_delta_rho_threshold"
    ),
    value = c(
      config$sim$t1,
      1,
      config$sim$params$t_drift,
      config$sim$params$c2_start,
      config$sim$params$c2_end,
      ifelse(is.finite(det$t_star), det$t_star, NA_real_),
      det$w,
      det$k,
      config$edm$simplex$E,
      config$edm$simplex$tau,
      gm_simplex1,
      gm_simplex10,
      config$edm$smap$E,
      gm$smap_theta_opt,
      gm$smap_rho_at_0,
      gm$smap_rho_at_opt,
      config$edm$windowed$window,
      config$edm$windowed$step,
      config$edm$windowed$theiler,
      ew_metrics$kendall_rho,
      ew_metrics$kendall_theta,
      ew_metrics$kendall_delta_rho,
      ew_metrics$kendall_b_median,
      ew_metrics$kendall_b_sd,
      ew_metrics$mean_advance_time,
      ew_metrics$thresholds$theta,
      ew_metrics$thresholds$delta_rho
    )
  )

  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write.csv(out, path, row.names = FALSE)
  invisible(out)
}

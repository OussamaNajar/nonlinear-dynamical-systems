# scripts/run_pipeline.R
# ==============================================================================
# Tri-Trophic Dynamics + EDM Early-Warning Pipeline
#
# What this pipeline produces (research contract)
# - Mechanistic tri-trophic simulation under slow forcing of c2(t)
# - Empirical transition time t* from variance proxy (rolling SD)
# - EDM baselines (global Simplex + global S-map theta sweep)
# - Windowed EDM diagnostics (Simplex skill decay; S-map theta_opt, Δρ, local slopes)
# - Quantified early-warning strength (Kendall trends + mean lead time)
# - Sensitivity table over key nuisance parameters (window sizes)
#
# Why this exists (portfolio + reviewer logic)
# - Anyone can generate plots. This pipeline generates *defensible evidence*:
#   (i) a mechanistic forcing scenario,
#   (ii) detection of a transition,
#   (iii) time-local predictability and nonlinearity diagnostics,
#   (iv) robustness to analysis hyperparameters.
#
# How to run
#   Rscript scripts/run_pipeline.R
#
# Outputs (disposable artifacts; safe to gitignore)
#   results/ : CSV tables + reproducibility metadata
#   figures/ : PNG/PDF figures
# ==============================================================================

suppressPackageStartupMessages({
  library(deSolve)   # ODE integration only
})

# ------------------------------------------------------------------------------
# Console logging (keeps output consistent and greppable)
# ------------------------------------------------------------------------------
log_info <- function(fmt, ...) message(sprintf(fmt, ...))

# ------------------------------------------------------------------------------
# Robust path resolution (no reliance on working directory)
# ------------------------------------------------------------------------------
get_script_path <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) return(normalizePath(sub("^--file=", "", file_arg[1])))
  stop("Run with: Rscript scripts/run_pipeline.R")
}

script_path  <- get_script_path()
project_root <- normalizePath(file.path(dirname(script_path), ".."))
r_dir        <- file.path(project_root, "R")

results_dir  <- file.path(project_root, "results")
figures_dir  <- file.path(project_root, "figures")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Load library modules (all analysis logic lives under R/)
# NOTE: plotting is canonical at R/plotting.R (single source of truth).
# ------------------------------------------------------------------------------
source(file.path(r_dir, "model_hp.R"))
source(file.path(r_dir, "simulate.R"))
source(file.path(r_dir, "transition_detect.R"))
source(file.path(r_dir, "edm_simplex.R"))
source(file.path(r_dir, "edm_smap.R"))
source(file.path(r_dir, "windowed_simplex.R"))
source(file.path(r_dir, "windowed_smap.R"))
source(file.path(r_dir, "plotting.R"))
source(file.path(r_dir, "summary_stats.R"))

# ------------------------------------------------------------------------------
# Experiment definition (treat as a paper appendix)
# ------------------------------------------------------------------------------
config <- list(
  sim = list(
    t0 = 0,
    t1 = 8000,
    dt = 1,
    # Initial state chosen to avoid trivial equilibria while remaining bounded.
    state0 = c(X = 0.6, Y = 0.2, Z = 0.1),
    params = list(
      # Hastings–Powell style tri-trophic structure
      c1 = 5,  h1 = 3,
      h2 = 2,
      m2 = 0.4, m3 = 0.008,
      # Slow forcing in predator–consumer interaction strength
      c2_start = 0.06,
      c2_end   = 0.18,
      t_drift  = 2000
    )
  ),
  detect = list(
    # Rolling SD is a variance proxy; w controls smoothing vs responsiveness.
    w = 300,
    # k controls false-positive rate (robust threshold via median + k*MAD).
    k = 3,
    # Baseline window chosen early before forcing alters the attractor.
    baseline_max_time = 700
  ),
  edm = list(
    series = "Z",  # top predator is typically the clearest early-warning channel
    simplex = list(E = 3, tau = 1, Tp = c(1, 10)),
    smap    = list(E = 3, tau = 1, Tp = 1,
                   thetas = seq(0, 8, by = 0.5),
                   theta_slope = 8),
    windowed = list(
      
      window = 800,
      step = 50,
      theiler = 50,   # avoid temporally adjacent analogs (leakage)
      lib_frac = 0.8  # explicit library/test split per window (out-of-sample)
    )
  ),
  # Minimal sensitivity analysis: robustness to two key hyperparameters:
  # - detection window w
  # - EDM window length
  sensitivity = list(
    detect_w = c(200, 300),
    window   = c(600, 800, 1000),
    lookback_window = 400
  )
)

# Ensure drift schedule has access to final horizon (hp_rhs_drift expects t1).
config$sim$params$t1 <- config$sim$t1

# ------------------------------------------------------------------------------
# Banner
# ------------------------------------------------------------------------------
cat("\n")
cat(strrep("=", 78), "\n", sep = "")
cat("TRI-TROPHIC + EDM EARLY-WARNING PIPELINE\n")
cat(strrep("=", 78), "\n", sep = "")
cat(sprintf("Project root: %s\n", project_root))
cat(sprintf("Results dir:  %s\n", results_dir))
cat(sprintf("Figures dir:  %s\n", figures_dir))
cat(strrep("=", 78), "\n\n", sep = "")

# ==============================================================================
# 1) Mechanistic simulation (ground truth dynamics under controlled forcing)
# ==============================================================================
log_info("[1/9] Simulating tri-trophic ODE with slow c2(t) forcing...")
times <- seq(config$sim$t0, config$sim$t1, by = config$sim$dt)
df <- simulate_hp(times, config$sim$state0, config$sim$params)
write.csv(df, file.path(results_dir, "sim.csv"), row.names = FALSE)

# ==============================================================================
# 2) Empirical transition detection (variance proxy)
# ==============================================================================
log_info("[2/9] Detecting transition using rolling SD (variance proxy)...")
series_full <- df[[config$edm$series]]

det <- detect_transition(
  time = df$time,
  series = series_full,
  w = config$detect$w,
  k = config$detect$k,
  baseline_max_time = config$detect$baseline_max_time
)

write.csv(
  data.frame(
    t_star = det$t_star,
    threshold = det$threshold,
    w = det$w,
    k = det$k,
    baseline_max_time = det$baseline_max_time
  ),
  file.path(results_dir, "transition.csv"),
  row.names = FALSE
)

log_info("    t* (detected) = %s", ifelse(is.finite(det$t_star), sprintf("%.0f", det$t_star), "NA"))

# ==============================================================================
# 3) Core diagnostic figures (minimal set, paper-ready)
# ==============================================================================
log_info("[3/9] Writing diagnostic figures...")
plot_timeseries(df, file.path(figures_dir, "01_timeseries_xyz.png"))
if ("c2" %in% names(df)) plot_c2(df, file.path(figures_dir, "02_c2_drift.png"))
plot_rollsd(df$time, det$rsd, det$threshold, det$t_star,
            file.path(figures_dir, "03_transition_rollsd.png"))

# Phase space is an interpretability tool: shows attractor geometry pre/post t*
plot_phase_space(df, file.path(figures_dir, "04_phase_space.png"), t_split = det$t_star)

# ==============================================================================
# 4) EDM analysis segment (post-burn-in / post-drift)
# ==============================================================================
t_drift <- config$sim$params$t_drift
idx_post <- df$time >= t_drift
if (!any(idx_post)) stop("No post-drift segment: check t_drift vs simulation horizon.")

df_post <- df[idx_post, , drop = FALSE]
series  <- df_post[[config$edm$series]]

log_info("[4/9] EDM segment: t >= %d (n=%d)", t_drift, length(series))

# ==============================================================================
# 5) Global Simplex baselines (sanity: predictability exists at some horizon)
# ==============================================================================
log_info("[5/9] Global Simplex baselines...")
global_metrics <- list()
E   <- config$edm$simplex$E
tau <- config$edm$simplex$tau

for (Tp in config$edm$simplex$Tp) {
  res <- simplex_forecast(series, E = E, tau = tau, Tp = Tp)
  met <- skill_metrics(res$obs, res$pred)

  global_metrics[[paste0("simplex_Tp", Tp)]] <- met

  write.csv(data.frame(rho = met$rho, rmse = met$rmse),
            file.path(results_dir, sprintf("simplex_Tp%d_metrics.csv", Tp)),
            row.names = FALSE)

  write.csv(data.frame(obs = res$obs, pred = res$pred),
            file.path(results_dir, sprintf("simplex_Tp%d_series.csv", Tp)),
            row.names = FALSE)

  plot_skill(res$obs, res$pred,
             file.path(figures_dir, sprintf("05_simplex_Tp%d.png", Tp)),
             sprintf("Simplex forecast (h=%d, rho=%.3f)", Tp, met$rho))

  log_info("    h=%d: rho=%.3f | rmse=%.4f", Tp, met$rho, met$rmse)
}

# ==============================================================================
# 6) Global S-map theta sweep (state-dependence diagnostic)
# ==============================================================================
log_info("[6/9] Global S-map theta sweep...")
theta_df <- theta_sweep(
  series,
  thetas = config$edm$smap$thetas,
  E = config$edm$smap$E,
  tau = config$edm$smap$tau,
  Tp = config$edm$smap$Tp
)
write.csv(theta_df, file.path(results_dir, "smap_theta_sweep.csv"), row.names = FALSE)

best_theta_idx <- which.max(theta_df$rho)
global_metrics$smap_theta_opt    <- theta_df$theta[best_theta_idx]
global_metrics$smap_rho_at_opt   <- theta_df$rho[best_theta_idx]
global_metrics$smap_rho_at_0     <- theta_df$rho[theta_df$theta == 0]

log_info("    theta_opt = %.1f (rho_opt=%.3f)", global_metrics$smap_theta_opt, global_metrics$smap_rho_at_opt)
log_info("    theta=0   = rho %.3f", global_metrics$smap_rho_at_0)

# ==============================================================================
# 7) Windowed Simplex (predictability loss under forcing)
# ==============================================================================
win     <- config$edm$windowed$window
stp     <- config$edm$windowed$step
theiler <- config$edm$windowed$theiler
lib_frac<- config$edm$windowed$lib_frac

log_info("[7/9] Windowed Simplex (W=%d, step=%d)...", win, stp)

df_ws <- windowed_simplex_skill(
  series,
  E = config$edm$simplex$E,
  tau = config$edm$simplex$tau,
  Tp_vec = config$edm$simplex$Tp,
  window = win,
  step = stp,
  theiler = theiler,
  lib_frac = lib_frac
)

df_ws$time <- df_post$time[df_ws$center_idx]
write.csv(df_ws, file.path(results_dir, "simplex_skill_windowed.csv"), row.names = FALSE)

plot_windowed_simplex(
  df_ws,
  Tp_vec = config$edm$simplex$Tp,
  path_pdf = file.path(figures_dir, "06_simplex_skill_windowed.pdf"),
  t_star = det$t_star
)

# ==============================================================================
# 8) Windowed S-map (time-local nonlinearity + stability proxies)
# ==============================================================================
log_info("[8/9] Windowed S-map diagnostics...")

df_wm <- windowed_smap_diagnostics(
  series,
  thetas = config$edm$smap$thetas,
  E = config$edm$smap$E,
  tau = config$edm$smap$tau,
  Tp = config$edm$smap$Tp,
  window = win,
  step = stp,
  theiler = theiler,
  lib_frac = lib_frac,
  theta_slope = config$edm$smap$theta_slope
)

df_wm$time <- df_post$time[df_wm$center_idx]
write.csv(df_wm, file.path(results_dir, "smap_windowed.csv"), row.names = FALSE)

plot_windowed_smap(
  df_wm,
  path_theta_pdf = file.path(figures_dir, "07_smap_theta_opt_windowed.pdf"),
  path_gain_pdf  = file.path(figures_dir, "08_smap_delta_rho_windowed.pdf"),
  t_star = det$t_star
)

plot_windowed_smap_slopes(
  df_wm,
  path_b_pdf  = file.path(figures_dir, "09_smap_slope_median.pdf"),
  path_sd_pdf = file.path(figures_dir, "10_smap_slope_sd.pdf"),
  t_star = det$t_star
)

# ==============================================================================
# 9) Quantified early-warning strength + summary panel
# ==============================================================================
log_info("[9/9] Early-warning metrics + summary panel...")

ew_metrics <- compute_early_warning_metrics(
  df_ws, df_wm, det, det$t_star,
  lookback_window = config$sensitivity$lookback_window
)

create_summary_table(
  config, det, global_metrics, ew_metrics,
  file.path(results_dir, "summary.csv")
)

plot_summary_panel(
  df, det, df_ws, df_wm,
  path = file.path(figures_dir, "11_summary_panel.png"),
  Tp = 1
)

# ==============================================================================
# Sensitivity analysis (robustness to nuisance hyperparameters)
# ==============================================================================
log_info("[9b/9] Sensitivity analysis (detect_w × EDM window)...")

sens_grid <- expand.grid(
  detect_w = config$sensitivity$detect_w,
  window   = config$sensitivity$window,
  stringsAsFactors = FALSE
)

sens_rows <- vector("list", nrow(sens_grid))

for (i in seq_len(nrow(sens_grid))) {
  w_i   <- sens_grid$detect_w[i]
  win_i <- sens_grid$window[i]

  det_i <- detect_transition(
    time = df$time,
    series = series_full,
    w = w_i,
    k = config$detect$k,
    baseline_max_time = config$detect$baseline_max_time
  )

  df_ws_i <- windowed_simplex_skill(
    series,
    E = config$edm$simplex$E,
    tau = config$edm$simplex$tau,
    Tp_vec = config$edm$simplex$Tp,
    window = win_i,
    step = stp,
    theiler = theiler,
    lib_frac = lib_frac
  )
  df_ws_i$time <- df_post$time[df_ws_i$center_idx]

  df_wm_i <- windowed_smap_diagnostics(
    series,
    thetas = config$edm$smap$thetas,
    E = config$edm$smap$E,
    tau = config$edm$smap$tau,
    Tp = config$edm$smap$Tp,
    window = win_i,
    step = stp,
    theiler = theiler,
    lib_frac = lib_frac,
    theta_slope = config$edm$smap$theta_slope
  )
  df_wm_i$time <- df_post$time[df_wm_i$center_idx]

  ew_i <- compute_early_warning_metrics(
    df_ws_i, df_wm_i, det_i, det_i$t_star,
    lookback_window = config$sensitivity$lookback_window
  )

  sens_rows[[i]] <- data.frame(
    detect_w = w_i,
    window   = win_i,
    t_star   = det_i$t_star,
    kendall_rho       = ew_i$kendall_rho,
    kendall_theta     = ew_i$kendall_theta,
    kendall_delta_rho = ew_i$kendall_delta_rho,
    kendall_b_median  = ew_i$kendall_b_median,
    kendall_b_sd      = ew_i$kendall_b_sd,
    mean_advance_time = ew_i$mean_advance_time
  )
}

sens_df <- do.call(rbind, sens_rows)
write.csv(sens_df, file.path(results_dir, "sensitivity_summary.csv"), row.names = FALSE)
log_info("    Saved: %s", file.path(results_dir, "sensitivity_summary.csv"))

# ==============================================================================
# Reproducibility metadata (minimum viable provenance)
# ==============================================================================
writeLines(capture.output(sessionInfo()), con = file.path(results_dir, "session_info.txt"))

cat("\n")
cat(strrep("=", 78), "\n", sep = "")
cat("PIPELINE COMPLETE\n")
cat(strrep("=", 78), "\n", sep = "")
cat(sprintf("Results: %s\n", results_dir))
cat(sprintf("Figures: %s\n", figures_dir))
cat(strrep("=", 78), "\n\n", sep = "")

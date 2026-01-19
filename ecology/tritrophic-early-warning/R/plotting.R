# R/plotting.R
# ==============================================================================
# Figure generation for tri-trophic + EDM early-warning pipeline (base R)
#
# Design goals
# - Deterministic, script-safe figure writing (always closes devices).
# - Minimal global graphics side effects (par restored).
# - Figures are treated as *evidence*: each plot supports an inference claim.
#
# What this file is NOT
# - A style library.
# - A GUI plotting toolkit.
# - A place to tune colors to "look nice".
#
# Contract
# - All plot_* functions validate inputs and write to disk.
# - save_png/save_pdf evaluate an expression in a device context and close it.
# ==============================================================================

# ------------------------------------------------------------------------------
# Device wrappers (hard requirement: never leave devices open)
# ------------------------------------------------------------------------------
save_png <- function(path, expr, width = 1200, height = 800, res = 150) {
  if (!is.character(path) || length(path) != 1) stop("save_png: path must be a single string.")
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

  png(path, width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)

  force(expr)
}

save_pdf <- function(path, expr, width = 8, height = 5) {
  if (!is.character(path) || length(path) != 1) stop("save_pdf: path must be a single string.")
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

  pdf(path, width = width, height = height, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)

  force(expr)
}

# ------------------------------------------------------------------------------
# Small helper: highlight pre-transition region to support "early warning" reading
# ------------------------------------------------------------------------------
shade_pre_transition <- function(t_star) {
  if (!is.finite(t_star)) return(invisible(NULL))
  usr <- par("usr")
  rect(usr[1], usr[3], t_star, usr[4], col = grDevices::adjustcolor("gray70", alpha.f = 0.15), border = NA)
  invisible(NULL)
}

# ------------------------------------------------------------------------------
# 1) Mechanistic dynamics: X(t), Y(t), Z(t)
# Claim supported: system exhibits structured dynamics under slow forcing.
# ------------------------------------------------------------------------------
plot_timeseries <- function(df, path) {
  req <- c("time", "X", "Y", "Z")
  if (!all(req %in% names(df))) stop("plot_timeseries: df must contain columns: time, X, Y, Z")

  save_png(path, width = 1400, height = 900, {
    oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)
    par(mfrow = c(3, 1), mar = c(4.5, 5, 2.5, 2))

    plot(df$time, df$X, type = "l", xlab = "", ylab = "X(t)", main = "Prey")
    grid()
    plot(df$time, df$Y, type = "l", xlab = "", ylab = "Y(t)", main = "Consumer")
    grid()
    plot(df$time, df$Z, type = "l", xlab = "Time", ylab = "Z(t)", main = "Top predator")
    grid()
  })
}

# ------------------------------------------------------------------------------
# 2) Forcing schedule: c2(t)
# Claim supported: exogenous slow drift is explicitly controlled and documented.
# ------------------------------------------------------------------------------
plot_c2 <- function(df, path) {
  req <- c("time", "c2")
  if (!all(req %in% names(df))) stop("plot_c2: df must contain columns: time, c2")

  save_png(path, {
    oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)
    par(mar = c(5, 5, 3, 2))
    plot(df$time, df$c2, type = "l", lwd = 2,
         xlab = "Time", ylab = expression(c[2](t)),
         main = expression("Slow forcing: " * c[2](t)))
    grid()
  })
}

# ------------------------------------------------------------------------------
# 3) Transition detection diagnostic: rolling SD + robust threshold + t*
# Claim supported: an empirical transition time can be estimated from variance proxy.
# ------------------------------------------------------------------------------
plot_rollsd <- function(time, rsd, thr, t_star, path) {
  if (!is.numeric(time) || !is.numeric(rsd) || length(time) != length(rsd)) {
    stop("plot_rollsd: time and rsd must be numeric vectors of the same length.")
  }
  if (!is.numeric(thr) || length(thr) != 1) stop("plot_rollsd: thr must be a scalar numeric.")
  if (!is.numeric(t_star) || length(t_star) != 1) stop("plot_rollsd: t_star must be a scalar numeric (NA allowed).")

  save_png(path, {
    oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)
    par(mar = c(5, 5, 3, 2))

    plot(time, rsd, type = "l", lwd = 2,
         xlab = "Time", ylab = "Rolling SD",
         main = "Variance proxy for regime shift")
    grid()
    abline(h = thr, lty = 2, lwd = 2)

    if (is.finite(t_star)) {
      shade_pre_transition(t_star)
      # redraw on top after shading
      lines(time, rsd, lwd = 2)
      abline(h = thr, lty = 2, lwd = 2)
      abline(v = t_star, lwd = 3)
    }

    legend("topleft",
           legend = c(paste0("thr=", signif(thr, 3)),
                      paste0("t*=", ifelse(is.finite(t_star), signif(t_star, 4), "NA"))),
           bty = "n")
  })
}

# ------------------------------------------------------------------------------
# 4) Forecast skill plot (observed vs predicted)
# Claim supported: model-free reconstruction has predictive information.
# ------------------------------------------------------------------------------
plot_skill <- function(obs, pred, path, title = "Forecast: observed vs predicted") {
  if (!is.numeric(obs) || !is.numeric(pred) || length(obs) != length(pred)) {
    stop("plot_skill: obs and pred must be numeric vectors of equal length.")
  }

  save_png(path, {
    oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)
    par(mar = c(5, 5, 3, 2))

    plot(obs, type = "l", lwd = 1.5,
         ylab = "Value", xlab = "Index", main = title)
    lines(pred, lwd = 1.5, lty = 2)
    grid()
    legend("topright", legend = c("obs", "pred"), lty = c(1, 2), bty = "n")
  })
}

# ------------------------------------------------------------------------------
# 5) Windowed Simplex forecast skill ρ(t)
# Claim supported: predictability changes systematically as forcing progresses.
# ------------------------------------------------------------------------------
plot_windowed_simplex <- function(df_skill, Tp_vec, path_pdf, t_star = NA_real_) {
  if (!("time" %in% names(df_skill))) stop("plot_windowed_simplex: df_skill must contain column 'time'.")
  for (Tp in Tp_vec) {
    col <- paste0("rho_Tp", Tp)
    if (!(col %in% names(df_skill))) stop(sprintf("plot_windowed_simplex: missing column '%s'.", col))
  }

  save_pdf(path_pdf, width = 10, height = 6, {
    oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)
    par(mar = c(5, 5, 3, 2))

    y0 <- df_skill[[paste0("rho_Tp", Tp_vec[1])]]
    plot(df_skill$time, y0, type = "l", lwd = 2,
         xlab = "Time", ylab = expression(rho),
         main = "Windowed Simplex forecast skill")
    grid()

    if (length(Tp_vec) > 1) {
      for (Tp in Tp_vec[-1]) {
        lines(df_skill$time, df_skill[[paste0("rho_Tp", Tp)]], lwd = 2, lty = 2)
      }
    }

    if (is.finite(t_star)) {
      shade_pre_transition(t_star)
      # redraw on top
      lines(df_skill$time, y0, lwd = 2)
      if (length(Tp_vec) > 1) {
        for (Tp in Tp_vec[-1]) lines(df_skill$time, df_skill[[paste0("rho_Tp", Tp)]], lwd = 2, lty = 2)
      }
      abline(v = t_star, lwd = 3)
    }

    legend("bottomleft", legend = paste0("h=", Tp_vec),
           lty = c(1, rep(2, length(Tp_vec) - 1)),
           bty = "n")
  })
}

# ------------------------------------------------------------------------------
# 6) Windowed S-map diagnostics: θ_opt(t) and Δρ(t)
# Claims supported:
# - θ_opt(t) increases => stronger state-dependence / local nonlinear behavior.
# - Δρ(t) > 0 => localized model improves over global linear baseline.
# ------------------------------------------------------------------------------
plot_windowed_smap <- function(df_smap, path_theta_pdf, path_gain_pdf, t_star = NA_real_) {
  req <- c("time", "theta_opt", "delta_rho")
  if (!all(req %in% names(df_smap))) {
    stop("plot_windowed_smap: df_smap must contain columns: time, theta_opt, delta_rho")
  }

  save_pdf(path_theta_pdf, width = 10, height = 6, {
    oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)
    par(mar = c(5, 5, 3, 2))

    plot(df_smap$time, df_smap$theta_opt, type = "l", lwd = 2,
         xlab = "Time", ylab = expression(theta[opt]),
         main = expression("S-map nonlinearity: " * theta[opt](t)))
    grid()
    if (is.finite(t_star)) {
      shade_pre_transition(t_star)
      lines(df_smap$time, df_smap$theta_opt, lwd = 2)
      abline(v = t_star, lwd = 3)
    }
  })

  save_pdf(path_gain_pdf, width = 10, height = 6, {
    oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)
    par(mar = c(5, 5, 3, 2))

    plot(df_smap$time, df_smap$delta_rho, type = "l", lwd = 2,
         xlab = "Time", ylab = expression(Delta * rho),
         main = expression("S-map predictive gain: " * Delta * rho(t)))
    grid()
    abline(h = 0, lty = 2)
    if (is.finite(t_star)) {
      shade_pre_transition(t_star)
      lines(df_smap$time, df_smap$delta_rho, lwd = 2)
      abline(h = 0, lty = 2)
      abline(v = t_star, lwd = 3)
    }
  })
}

# ------------------------------------------------------------------------------
# 7) Windowed S-map slope diagnostics: median b(t), dispersion sd(b)
# Claim supported:
# - Slope summaries act as stability/interaction proxies in the local linearization.
#   (Interpretation is system-specific; treat as diagnostic, not theorem.)
# ------------------------------------------------------------------------------
plot_windowed_smap_slopes <- function(df_smap, path_b_pdf, path_sd_pdf, t_star = NA_real_) {
  req <- c("time", "b_median", "b_sd")
  if (!all(req %in% names(df_smap))) {
    stop("plot_windowed_smap_slopes: missing required columns: time, b_median, b_sd")
  }

  save_pdf(path_b_pdf, width = 10, height = 6, {
    oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)
    par(mar = c(5, 5, 3, 2))

    plot(df_smap$time, df_smap$b_median, type = "l", lwd = 2,
         xlab = "Time", ylab = expression(b(t)),
         main = "S-map local slope: median")
    grid()
    abline(h = c(-1, 1), lty = 2)
    if (is.finite(t_star)) {
      shade_pre_transition(t_star)
      lines(df_smap$time, df_smap$b_median, lwd = 2)
      abline(h = c(-1, 1), lty = 2)
      abline(v = t_star, lwd = 3)
    }
  })

  save_pdf(path_sd_pdf, width = 10, height = 6, {
    oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)
    par(mar = c(5, 5, 3, 2))

    plot(df_smap$time, df_smap$b_sd, type = "l", lwd = 2,
         xlab = "Time", ylab = "sd(b)",
         main = "S-map local slope dispersion")
    grid()
    if (is.finite(t_star)) {
      shade_pre_transition(t_star)
      lines(df_smap$time, df_smap$b_sd, lwd = 2)
      abline(v = t_star, lwd = 3)
    }
  })
}

# ------------------------------------------------------------------------------
# 8) Phase space reconstruction (2D delay embedding)
# Claim supported: qualitative attractor geometry differs across regimes.
# Notes:
# - This is *illustrative*, not a substitute for embedding diagnostics.
# - Uses lag-1 plot Z(t) vs Z(t+1); intentionally simple and stable.
# ------------------------------------------------------------------------------
plot_phase_space <- function(df, path, t_split = NULL) {
  req <- c("time", "Z")
  if (!all(req %in% names(df))) stop("plot_phase_space: df must contain columns: time, Z")

  n <- nrow(df)
  if (n < 3) stop("plot_phase_space: df must have at least 3 rows.")

  Zt   <- df$Z[1:(n - 1)]
  Ztp1 <- df$Z[2:n]
  tt   <- df$time[2:n]

  save_png(path, width = 1200, height = 1200, {
    oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)
    par(mar = c(5, 5, 3, 2))

    if (!is.null(t_split) && is.finite(t_split)) {
      idx_pre  <- tt < t_split
      idx_post <- tt >= t_split

      plot(Zt[idx_pre], Ztp1[idx_pre], type = "l", lwd = 0.6,
           xlab = expression(Z[t]), ylab = expression(Z[t+1]),
           main = "Phase space (lag-1): pre vs post transition")
      grid()
      if (any(idx_post)) lines(Zt[idx_post], Ztp1[idx_post], lwd = 0.6, lty = 2)

      legend("topleft",
             legend = c("pre", "post"),
             lty = c(1, 2), bty = "n")
    } else {
      plot(Zt, Ztp1, type = "l", lwd = 0.6,
           xlab = expression(Z[t]), ylab = expression(Z[t+1]),
           main = "Phase space (lag-1)")
      grid()
    }
  })
}

# ------------------------------------------------------------------------------
# 9) Summary panel (single figure for "what happened?" narrative)
# Claim supported: early warning indicators are temporally aligned with forcing and t*.
# ------------------------------------------------------------------------------
plot_summary_panel <- function(df, det, df_ws, df_wm, path, Tp = 1) {
  if (!all(c("time", "Z") %in% names(df))) stop("plot_summary_panel: df must contain columns time,Z.")
  if (!all(c("time") %in% names(df_ws))) stop("plot_summary_panel: df_ws must contain column time.")
  if (!all(c("time", "theta_opt", "delta_rho", "b_median") %in% names(df_wm))) {
    stop("plot_summary_panel: df_wm must contain time, theta_opt, delta_rho, b_median.")
  }

  rho_col <- paste0("rho_Tp", Tp)
  if (!(rho_col %in% names(df_ws))) stop(sprintf("plot_summary_panel: df_ws missing %s", rho_col))

  t_star <- det$t_star

  save_png(path, width = 1600, height = 1200, {
    oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)

    layout(matrix(1:6, nrow = 3, byrow = TRUE))
    par(mar = c(4.2, 4.6, 2.6, 1))

    # 1) Z(t)
    plot(df$time, df$Z, type = "l", lwd = 1.2, xlab = "Time", ylab = "Z(t)", main = "Top predator")
    grid()
    if (is.finite(t_star)) abline(v = t_star, lwd = 2)

    # 2) Rolling SD
    plot(df$time, det$rsd, type = "l", lwd = 1.2, xlab = "Time", ylab = "Rolling SD", main = "Variance proxy")
    grid()
    abline(h = det$threshold, lty = 2)
    if (is.finite(t_star)) abline(v = t_star, lwd = 2)

    # 3) Simplex skill
    plot(df_ws$time, df_ws[[rho_col]], type = "l", lwd = 1.2, xlab = "Time", ylab = expression(rho),
         main = paste0("Simplex skill (h=", Tp, ")"))
    grid()
    if (is.finite(t_star)) abline(v = t_star, lwd = 2)

    # 4) theta_opt
    plot(df_wm$time, df_wm$theta_opt, type = "l", lwd = 1.2, xlab = "Time",
         ylab = expression(theta[opt]), main = "S-map nonlinearity")
    grid()
    if (is.finite(t_star)) abline(v = t_star, lwd = 2)

    # 5) delta_rho
    plot(df_wm$time, df_wm$delta_rho, type = "l", lwd = 1.2, xlab = "Time",
         ylab = expression(Delta * rho), main = "S-map gain")
    grid()
    abline(h = 0, lty = 2)
    if (is.finite(t_star)) abline(v = t_star, lwd = 2)

    # 6) b_median
    plot(df_wm$time, df_wm$b_median, type = "l", lwd = 1.2, xlab = "Time",
         ylab = "b(t)", main = "S-map local slope")
    grid()
    abline(h = c(-1, 1), lty = 2)
    if (is.finite(t_star)) abline(v = t_star, lwd = 2)
  })
}

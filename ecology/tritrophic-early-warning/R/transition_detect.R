# R/transition_detect.R
# ------------------------------------------------------------------------------
# Regime-shift indicator via rolling standard deviation (variance proxy)
#
# Purpose
# - Provide a lightweight variance-based indicator that can flag regime changes.
# - Baseline distribution is estimated from an early-time window; threshold is
#   set using robust statistics (median + k*MAD).
#
# Contract
# - detect_transition(time, series, w, k, baseline_max_time) -> list:
#     rsd, threshold, t_star, w, k, baseline_max_time
#
# Assumptions
# - series contains finite numeric values (no NA/Inf). We enforce this strictly.
# - rolling SD is computed with an O(n) sliding-window variance using running sums.
# ------------------------------------------------------------------------------

rolling_sd <- function(x, w) {
  if (!is.numeric(x)) stop("rolling_sd: x must be numeric.")
  n <- length(x)
  if (!is.numeric(w) || length(w) != 1 || w < 2 || w > n) stop("rolling_sd: invalid window w.")
  if (any(!is.finite(x))) stop("rolling_sd: x must be finite (no NA/Inf).")

  rsd <- rep(NA_real_, n)

  # Running sums for sliding window variance
  # Var = (sum(x^2) - sum(x)^2 / w) / (w - 1)
  s  <- 0.0
  s2 <- 0.0

  # initialize with first window
  for (i in 1:w) {
    s  <- s  + x[i]
    s2 <- s2 + x[i] * x[i]
  }

  v <- (s2 - (s * s) / w) / (w - 1)
  rsd[w] <- sqrt(max(v, 0))

  if (w < n) {
    for (i in (w + 1):n) {
      x_out <- x[i - w]
      x_in  <- x[i]

      s  <- s  - x_out + x_in
      s2 <- s2 - x_out * x_out + x_in * x_in

      v <- (s2 - (s * s) / w) / (w - 1)
      rsd[i] <- sqrt(max(v, 0))
    }
  }

  rsd
}

detect_transition <- function(time, series, w = 300, k, baseline_max_time = NULL) {
  if (!is.numeric(time) || !is.numeric(series)) stop("detect_transition: time/series must be numeric.")
  if (length(time) != length(series)) stop("detect_transition: time and series must have same length.")
  if (any(!is.finite(time))) stop("detect_transition: time must be finite.")
  if (any(!is.finite(series))) stop("detect_transition: series must be finite (no NA/Inf).")
  if (k <= 0) stop("detect_transition: k must be positive.")

  rsd <- rolling_sd(series, w)

  # Baseline selection
  if (is.null(baseline_max_time)) {
    first_ok <- which(is.finite(rsd))[1]
    if (is.na(first_ok)) stop("detect_transition: rolling_sd produced no finite values.")
    baseline_max_time <- time[first_ok]
  }

  idx_base <- which(time <= baseline_max_time & is.finite(rsd))
  if (length(idx_base) < max(20, w)) {
    stop("detect_transition: baseline too small; increase baseline_max_time or reduce w.")
  }

  base <- rsd[idx_base]

  # Robust threshold
  thr <- median(base, na.rm = TRUE) + k * mad(base, constant = 1, na.rm = TRUE)

  idx_cross <- which(rsd > thr)
  t_star <- if (length(idx_cross) == 0) NA_real_ else time[min(idx_cross)]

  list(rsd = rsd, threshold = thr, t_star = t_star, w = w, k = k, baseline_max_time = baseline_max_time)
}

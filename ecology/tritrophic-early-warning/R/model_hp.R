# R/model_hp.R
# ------------------------------------------------------------------------------
# Tri-trophic ODE model with time-varying interaction strength
#
# Purpose
# - Provide a mechanistic tri-trophic system with a controlled slow forcing of a
#   key interaction parameter c2(t). This enables stress-testing:
#   - variance-based shift indicators
#   - forecasting skill under gradual regime change
#
# Contract
# - hp_rhs_drift(t, state, params) returns:
#   list(c(dX, dY, dZ), c(c2 = c2_t))
#
# Assumptions
# - state is named numeric vector with X, Y, Z
# - params contains: c1,h1,h2,m2,m3,c2_start,c2_end,t_drift,t1
# ------------------------------------------------------------------------------

hp_rhs_drift <- function(t, state, params) {
  # Fail fast on obvious misuse
  if (!is.numeric(state) || !all(c("X", "Y", "Z") %in% names(state))) {
    stop("hp_rhs_drift: state must be a named numeric vector with X, Y, Z.")
  }
  if (!is.list(params)) stop("hp_rhs_drift: params must be a list.")

  with(as.list(c(state, params)), {

    if (!is.finite(t1) || !is.finite(t_drift) || t1 <= t_drift) {
      stop("hp_rhs_drift: require finite t1 > t_drift.")
    }
    if (c2_start <= 0 || c2_end <= 0) stop("hp_rhs_drift: c2_start/c2_end must be positive.")

    # Drift schedule:
    # - constant until t_drift (baseline regime)
    # - linear drift to c2_end by t1 (slow forcing)
    if (t <= t_drift) {
      c2_t <- c2_start
    } else {
      frac <- (t - t_drift) / (t1 - t_drift)
      frac <- max(0, min(1, frac))
      c2_t <- c2_start + frac * (c2_end - c2_start)
    }

    dX <- X * (1 - X) - c1 * X * Y / (1 + h1 * X)
    dY <- c1 * X * Y / (1 + h1 * X) - c2_t * Y * Z / (1 + h2 * Y) - m2 * Y
    dZ <- c2_t * Y * Z / (1 + h2 * Y) - m3 * Z

    list(c(dX, dY, dZ), c(c2 = c2_t))
  })
}

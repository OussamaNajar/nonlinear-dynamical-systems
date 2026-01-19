# R/simulate.R
# ------------------------------------------------------------------------------
# Simulation wrapper (ODE integration)
#
# Purpose
# - Provide a single, validated entry for integrating the tri-trophic ODE model.
#
# Contract
# - simulate_hp(times, state0, params) -> data.frame:
#     columns: time, X, Y, Z, (and any "out" variables like c2)
#
# Failure modes
# - stops early if hp_rhs_drift is missing or inputs are inconsistent.
# ------------------------------------------------------------------------------

simulate_hp <- function(times, state0, params) {
  if (!is.numeric(times) || length(times) < 3) stop("simulate_hp: times must be numeric with length >= 3.")
  if (!is.numeric(state0) || !all(c("X", "Y", "Z") %in% names(state0))) {
    stop("simulate_hp: state0 must be named numeric vector with X, Y, Z.")
  }
  if (!is.list(params)) stop("simulate_hp: params must be a list.")

  if (!exists("hp_rhs_drift", mode = "function")) {
    stop("simulate_hp: hp_rhs_drift() not found. Source R/model_hp.R first.")
  }

  out <- deSolve::ode(
    y = state0,
    times = times,
    func = hp_rhs_drift,
    parms = params
  )

  as.data.frame(out)
}

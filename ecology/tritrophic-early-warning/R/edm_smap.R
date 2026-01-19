# R/edm_smap.R
# ------------------------------------------------------------------------------
# EDM: S-map (locally weighted linear regression)
#
# What this file does
# - Implements a *multivariate* S-map forecaster consistent with standard EDM:
#     - delay embedding in R^E via make_embed()
#     - weighted least squares on [1, embedded_state]
#     - theta controls locality (theta=0 => global linear; theta>0 => local)
# - Provides a theta sweep utility for diagnosing state dependence.
#
# Why multivariate matters
# - S-map is defined on the reconstructed state-space. If you embed with E>1,
#   the local linear model should use *all* E coordinates. A 1D regression on
#   only the first coordinate is not an S-map in E dimensions.
#
# Contract
# - Requires make_embed() and skill_metrics() from R/edm_simplex.R
# - smap_forecast(x, E, tau, Tp, theta, theiler) -> list(obs, pred, coef)
#     coef is an (n_pred x (E+1)) matrix of local coefficients (intercept + slopes)
# - theta_sweep(x, thetas, E, tau, Tp, theiler)  -> data.frame(theta, rho, rmse)
# ------------------------------------------------------------------------------

# internal: weighted least squares without forming diag(W)
# beta = (X'WX)^{-1} X'W y
wls_beta <- function(X, y, w) {
  if (!is.matrix(X)) stop("wls_beta: X must be a matrix")
  if (!is.numeric(y) || !is.numeric(w)) stop("wls_beta: y and w must be numeric")
  if (nrow(X) != length(y) || length(y) != length(w)) stop("wls_beta: dimension mismatch")

  w[!is.finite(w)] <- 0
  if (sum(w) <= 0) return(NULL)

  XtW <- t(X) * w
  G   <- XtW %*% X
  rhs <- XtW %*% y

  beta <- tryCatch(solve(G, rhs), error = function(e) NULL)
  if (is.null(beta) || any(!is.finite(beta))) return(NULL)
  beta
}

smap_forecast <- function(x, E, tau = 1, Tp = 1, theta = 2, theiler = 0) {
  if (!is.numeric(x)) stop("smap_forecast: x must be numeric.")
  if (any(!is.finite(x))) stop("smap_forecast: x must be finite (no NA/Inf).")
  if (!is.numeric(E) || E < 2) stop("smap_forecast: E must be >= 2.")
  if (!is.numeric(tau) || tau < 1) stop("smap_forecast: tau must be >= 1.")
  if (!is.numeric(Tp) || Tp < 1) stop("smap_forecast: Tp must be >= 1.")
  if (!is.numeric(theta) || theta < 0) stop("smap_forecast: theta must be non-negative.")
  if (!is.numeric(theiler) || theiler < 0) stop("smap_forecast: theiler must be >= 0.")
  if (!exists("make_embed", mode = "function")) stop("smap_forecast: make_embed() not found.")

  M  <- make_embed(x, E, tau)
  nM <- nrow(M)

  target_idx <- (1:nM) + (E - 1) * tau + Tp
  valid <- which(target_idx <= length(x))
  if (length(valid) < (E + 10)) stop("smap_forecast: not enough points for stable fits.")

  Mv <- M[valid, , drop = FALSE]
  y  <- x[target_idx[valid]]

  nV <- nrow(Mv)
  preds <- rep(NA_real_, nV)
  coef  <- matrix(NA_real_, nrow = nV, ncol = E + 1)
  colnames(coef) <- c("a0", paste0("b", 1:E))

  # library uses all points except the prediction point i (LOO style)
  Xall <- cbind(1, Mv)

  for (i in 1:nV) {
    d <- sqrt(rowSums((Mv - Mv[i, ])^2))
    d[i] <- Inf

    # Theiler exclusion in embedded-index coordinates
    if (theiler > 0) {
      time_dist <- abs((1:nV) - i)
      d[time_dist <= theiler] <- Inf
    }

    finite_d <- d[is.finite(d)]
    if (length(finite_d) < (E + 1)) next

    dbar <- mean(finite_d)
    if (!is.finite(dbar) || dbar <= 0) next

    w <- exp(-theta * d / dbar)
    w[!is.finite(w)] <- 0

    beta <- wls_beta(Xall, y, w)
    if (is.null(beta)) next

    preds[i] <- sum(Xall[i, ] * beta)
    coef[i, ] <- as.numeric(beta)
  }

  list(obs = y, pred = preds, coef = coef)
}

theta_sweep <- function(x, thetas, E, tau = 1, Tp = 1, theiler = 0) {
  if (!is.numeric(thetas) || length(thetas) < 1) stop("theta_sweep: provide at least one theta.")
  if (any(thetas < 0)) stop("theta_sweep: theta values must be non-negative.")
  if (!exists("skill_metrics", mode = "function")) stop("theta_sweep: skill_metrics() not found.")

  out <- data.frame(theta = thetas, rho = NA_real_, rmse = NA_real_)
  for (i in seq_along(thetas)) {
    res <- smap_forecast(x, E = E, tau = tau, Tp = Tp, theta = thetas[i], theiler = theiler)
    met <- skill_metrics(res$obs, res$pred)
    out$rho[i]  <- met$rho
    out$rmse[i] <- met$rmse
  }
  out
}

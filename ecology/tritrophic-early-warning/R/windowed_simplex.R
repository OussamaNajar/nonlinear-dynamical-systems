# R/windowed_simplex.R
# ------------------------------------------------------------------------------
# Windowed EDM: Simplex forecast skill over time
#
# Purpose
# - Track predictability changes under slow forcing by computing out-of-sample
#   forecast skill on a sliding window across the time series.
# - Uses a paper-faithful protocol:
#   - delay-coordinate embedding in R^E (make_embed)
#   - E+1 nearest neighbors (simplex) with exponential distance weighting
#   - explicit library/test split within each window (lib_frac)
#   - Theiler exclusion to avoid temporally adjacent analogs
#
# Contract
# - windowed_simplex_skill(x, E, tau, Tp_vec, window, step, theiler, lib_frac) -> data.frame
#   with columns:
#     center_idx,
#     rho_Tp{Tp}, rmse_Tp{Tp} for each Tp in Tp_vec
#
# Assumptions
# - x is finite numeric (no NA/Inf)
# - make_embed() and skill_metrics() are available (from R/edm_simplex.R)
# ------------------------------------------------------------------------------

windowed_simplex_skill <- function(x,
                                   E = 3, tau = 1, Tp_vec = c(1, 10),
                                   window = 800, step = 50,
                                   theiler = 50, lib_frac = 0.8) {
  if (!is.numeric(x) || any(!is.finite(x))) stop("windowed_simplex_skill: x must be finite numeric.")
  if (!exists("make_embed", mode = "function")) stop("windowed_simplex_skill: make_embed() not found.")
  if (!exists("skill_metrics", mode = "function")) stop("windowed_simplex_skill: skill_metrics() not found.")
  if (!(lib_frac > 0 && lib_frac < 1)) stop("windowed_simplex_skill: lib_frac must be in (0,1).")

  n <- length(x)
  if (window < 50 || window >= n) stop("windowed_simplex_skill: invalid window.")
  if (step < 1) stop("windowed_simplex_skill: step must be >= 1.")

  starts <- seq(1, n - window + 1, by = step)

  out <- data.frame(center_idx = integer(0))
  for (Tp in Tp_vec) {
    out[[paste0("rho_Tp", Tp)]] <- numeric(0)
    out[[paste0("rmse_Tp", Tp)]] <- numeric(0)
  }

  for (s in starts) {
    e <- s + window - 1
    seg <- x[s:e]

    # embedding in window coordinates
    M <- make_embed(seg, E, tau)
    nM <- nrow(M)

    row <- list(center_idx = as.integer(round((s + e) / 2)))

    for (Tp in Tp_vec) {
      # horizon target index within window
      target_idx <- (1:nM) + (E - 1) * tau + Tp
      valid <- which(target_idx <= length(seg))
      Mv <- M[valid, , drop = FALSE]
      y  <- seg[target_idx[valid]]

      nV <- length(valid)
      nLib <- max(10, floor(lib_frac * nV))
      lib_ix <- 1:nLib
      pred_ix <- (nLib + 1):nV

      if (length(pred_ix) < 5) {
        row[[paste0("rho_Tp", Tp)]] <- NA_real_
        row[[paste0("rmse_Tp", Tp)]] <- NA_real_
        next
      }

      preds <- rep(NA_real_, length(pred_ix))

      for (j in seq_along(pred_ix)) {
        i <- pred_ix[j]

        # distances from prediction point i to library points
        d <- sqrt(rowSums((Mv[lib_ix, , drop = FALSE] - Mv[i, ])^2))

        # Theiler exclusion in window index coordinates
        time_dist <- abs(lib_ix - i)
        d[time_dist <= theiler] <- Inf

        # Require enough finite neighbors
        ord <- order(d)
        ord <- ord[is.finite(d[ord])]
        if (length(ord) < (E + 1)) next

        nn <- ord[1:(E + 1)]
        d_nn <- d[nn]
        d1 <- d_nn[1]

        w <- if (d1 == 0) rep(1, length(d_nn)) else exp(-d_nn / d1)

        y_lib <- y[lib_ix]
        y_nn <- y_lib[nn]

        preds[j] <- sum(w * y_nn) / sum(w)
      }

      met <- skill_metrics(y[pred_ix], preds)
      row[[paste0("rho_Tp", Tp)]] <- met$rho
      row[[paste0("rmse_Tp", Tp)]] <- met$rmse
    }

    out <- rbind(out, as.data.frame(row))
  }

  out
}

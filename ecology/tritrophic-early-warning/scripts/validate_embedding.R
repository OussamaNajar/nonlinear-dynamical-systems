# scripts/validate_embedding.R
# ==============================================================================
# Embedding dimension validation (False Nearest Neighbors; Kennel et al. 1992)
#
# Purpose (what claim this supports)
# - Provide empirical support for choosing embedding dimension E used by EDM
#   analyses in the main pipeline (Simplex + S-map).
#
# What this is (and is not)
# - A diagnostic: FNN(E) should drop sharply near the minimum adequate E.
# - Not a theorem: results depend on noise, sampling rate, and the chosen series.
#
# Run
#   Rscript scripts/validate_embedding.R
#
# Outputs
#   results/embedding_validation_fnn.csv
#   figures/00_embedding_validation_fnn.png
# ==============================================================================

log_info <- function(fmt, ...) message(sprintf(fmt, ...))

get_script_path <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) return(normalizePath(sub("^--file=", "", file_arg[1])))
  if (interactive()) return(normalizePath(file.path(getwd(), "scripts", "validate_embedding.R")))
  stop("Run with: Rscript scripts/validate_embedding.R")
}

script_path  <- get_script_path()
project_root <- normalizePath(file.path(dirname(script_path), ".."))
r_dir        <- file.path(project_root, "R")

results_dir  <- file.path(project_root, "results")
figures_dir  <- file.path(project_root, "figures")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Core modules (no extra deps; simulate_hp uses deSolve::ode internally)
source(file.path(r_dir, "model_hp.R"))
source(file.path(r_dir, "simulate.R"))
source(file.path(r_dir, "edm_simplex.R"))  # for make_embed()

# ------------------------------------------------------------------------------
# False nearest neighbors (FNN)
#
# Method sketch
# - For each E in {2, ..., E_max-1}:
#   - Embed x into R^E via delay coordinates (lag tau)
#   - For a subset of points, find nearest neighbor in R^E (excluding itself)
#   - Evaluate whether that neighbor remains close when lifted to R^(E+1)
#
# Notes
# - We start at E=2 because the shared embedding helper make_embed() enforces E>=2.
# - We subsample points (stride) to reduce O(n^2) runtime.
# - We optionally use a Theiler exclusion to reduce temporal adjacency artifacts.
# ------------------------------------------------------------------------------
false_nearest_neighbors <- function(x,
                                   E_max = 10,
                                   tau = 1,
                                   rtol = 15,
                                   atol = 2,
                                   stride = 5,
                                   theiler = 0) {

  if (!is.numeric(x) || any(!is.finite(x))) stop("FNN: x must be finite numeric.")
  if (!is.numeric(E_max) || E_max < 3 || E_max > 15) stop("FNN: E_max must be in [3, 15].")
  if (!is.numeric(tau) || tau < 1) stop("FNN: tau must be >= 1.")
  if (!is.numeric(stride) || stride < 1) stop("FNN: stride must be >= 1.")
  if (!is.numeric(theiler) || theiler < 0) stop("FNN: theiler must be >= 0.")

  n <- length(x)
  E_vals <- 2:(E_max - 1)
  out <- data.frame(E = E_vals, fnn_percent = NA_real_, n_tested = 0L)

  # Attractor scale proxy (global): typical radius from mean in 1D.
  # Used only for the "absolute" Kennel criterion.
  R_A <- sqrt(mean((x - mean(x))^2))

  for (E in E_vals) {
    # Build E-dimensional embedding using shared helper (consistent with EDM codebase)
    M_E <- make_embed(x, E = E, tau = tau)
    m <- nrow(M_E)
    if (m < 100) {
      warning(sprintf("FNN: insufficient embedded points for E=%d (m=%d).", E, m))
      out$fnn_percent[out$E == E] <- NA_real_
      out$n_tested[out$E == E] <- 0L
      next
    }

    # Coordinate in the (E+1)th dimension aligned with rows of M_E:
    # next coordinate is x shifted by E*tau.
    idx_future <- (1:m) + E * tau
    valid <- idx_future <= n
    M_E <- M_E[valid, , drop = FALSE]
    x_future <- x[idx_future[valid]]
    m <- nrow(M_E)
    if (m < 100) {
      out$fnn_percent[out$E == E] <- NA_real_
      out$n_tested[out$E == E] <- 0L
      next
    }

    test_idx <- seq(1, m, by = stride)

    n_false  <- 0L
    n_tested <- 0L

    for (ii in test_idx) {
      d_E <- sqrt(rowSums((M_E - M_E[ii, ])^2))

      # Exclude self + Theiler window (reduces temporal adjacency artifacts)
      d_E[ii] <- Inf
      if (theiler > 0) {
        lo <- max(1, ii - theiler)
        hi <- min(m, ii + theiler)
        d_E[lo:hi] <- Inf
      }

      nn <- which.min(d_E)
      d0 <- d_E[nn]
      if (!is.finite(d0) || d0 <= 0) next

      # Lifted distance in R^(E+1)
      d_E1 <- sqrt(d0^2 + (x_future[ii] - x_future[nn])^2)

      # Kennel criteria (relative + absolute)
      rel <- (d_E1 / d0) > rtol
      abs <- (R_A > 0) && ((d_E1 / R_A) > atol)

      if (rel || abs) n_false <- n_false + 1L
      n_tested <- n_tested + 1L
    }

    out$fnn_percent[out$E == E] <- if (n_tested > 0) 100 * n_false / n_tested else NA_real_
    out$n_tested[out$E == E] <- n_tested
  }

  out
}

# ------------------------------------------------------------------------------
# Configuration (mirrors main pipeline defaults; keep in sync if you change model)
# ------------------------------------------------------------------------------
config <- list(
  t0 = 0,
  t1 = 8000,
  dt = 1,
  state0 = c(X = 0.6, Y = 0.2, Z = 0.1),
  params = list(
    c1 = 5, h1 = 3, h2 = 2,
    m2 = 0.4, m3 = 0.008,
    c2_start = 0.06, c2_end = 0.18,
    t_drift = 2000
  )
)
config$params$t1 <- config$t1

cat("\n", strrep("=", 78), "\n", sep = "")
cat("EMBEDDING DIMENSION VALIDATION (FNN)\n")
cat(strrep("=", 78), "\n\n", sep = "")

times <- seq(config$t0, config$t1, by = config$dt)

log_info("Simulating tri-trophic system...")
df <- simulate_hp(times, config$state0, config$params)

# Analyze the same segment EDM uses (post-drift)
idx_post <- df$time >= config$params$t_drift
Z <- df$Z[idx_post]
log_info("Post-drift segment: n=%d", length(Z))

log_info("Computing FNN(E) for E=2..%d (tau=%d).", 9, 1)
fnn <- false_nearest_neighbors(
  Z, E_max = 10, tau = 1,
  rtol = 15, atol = 2,
  stride = 5, theiler = 50
)

write.csv(fnn, file.path(results_dir, "embedding_validation_fnn.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# Plot (evidence-first; no palette theatrics)
# ------------------------------------------------------------------------------
png(file.path(figures_dir, "00_embedding_validation_fnn.png"),
    width = 1100, height = 750, res = 150)
oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)
par(mar = c(5, 5, 3, 2))

plot(fnn$E, fnn$fnn_percent, type = "b", pch = 19, lwd = 2,
     xlab = "Embedding dimension E",
     ylab = "False nearest neighbors (%)",
     main = "FNN embedding validation (post-drift Z(t))")
grid()
abline(h = 5, lty = 2, lwd = 2)
legend("topright",
       legend = c("FNN %", "5% reference"),
       lty = c(1, 2), pch = c(19, NA), bty = "n")
dev.off()

cat("Saved:\n")
cat(sprintf("  Data: %s\n", file.path(results_dir, "embedding_validation_fnn.csv")))
cat(sprintf("  Plot: %s\n", file.path(figures_dir, "00_embedding_validation_fnn.png")))
cat(strrep("=", 78), "\n\n", sep = "")

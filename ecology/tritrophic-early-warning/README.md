# Tri-trophic early warning via EDM (Hastings–Powell)

A **fully reproducible R pipeline** that simulates a forced Hastings–Powell tri-trophic food chain and tests whether **EDM diagnostics** provide **early warning** of an impending stability loss (equilibrium → sustained oscillations) using **only one observed variable** (top predator **Z**).

---

## What this repository does

**Mechanistic experiment**
1. Simulate a tri-trophic ODE with slow drift in interaction strength **c2(t)**.
2. Detect an empirical transition time (t\*) using rolling SD of **Z(t)** (variance proxy).


** EDM diagnostics (post-drift segment) **
3. Global **Simplex** baselines (h = 1 and h = 10).
4. Global **S-map** theta sweep (state dependence check).
5. Sliding-window **Simplex** skill (h = 1 and h = 10).
6. Sliding-window **S-map**: theta\_opt(t), Δρ(t), and slope summaries b(t).

**Robustness**
7. Sensitivity grid over detection window **w** × EDM window length **W**.

---

## Quick start

### 1) Install dependencies
This project intentionally uses minimal dependencies.
You only need:
- R (>= 4.1 recommended)
- `deSolve`

If you don't have `deSolve`:
```r
install.packages("deSolve")
```

### 2) Run the full pipeline
From the project root:
```bash
Rscript scripts/run_pipeline.R
Rscript scripts/validate_embedding.R
```

**Expected runtime:** ~3-4 minutes on modern hardware

---

## Repository layout

- `scripts/` : entry points you run
- `R/` : modular functions (simulation, EDM, plotting, metrics)
- `figures/` : generated paper-ready plots
- `results/` : generated CSV outputs + provenance
- `docs/` : interpretation guides and assumptions
- `paper/` : manuscript PDF

---

## Outputs

The pipeline generates **12 publication-quality figures** and **12 CSV tables**.

### Key Figures

**Diagnostic Overview:**
- `11_summary_panel.png` : **Six-panel dashboard** showing all coordinated early warning indicators

**System Dynamics:**
- `01_timeseries_xyz.png` : Complete tri-trophic time series (X, Y, Z)
- `02_c2_drift.png` : Parameter forcing schedule c₂(t)
- `03_transition_rollsd.png` : Transition detection via rolling standard deviation
- `04_phase_space.png` : Phase-space comparison (pre- vs post-transition)

**EDM Validation:**
- `00_embedding_validation_fnn.png` : False Nearest Neighbors analysis (E=3 optimal)
- `05_simplex_Tp1.png` : Simplex forecast visualization (h=1)
- `05_simplex_Tp10.png` : Simplex forecast visualization (h=10)

**Windowed Diagnostics (PDFs):**
- `06_simplex_skill_windowed.pdf` : Forecast skill degradation over time
- `07_smap_theta_opt_windowed.pdf` : Optimal nonlinearity parameter θ*(t)
- `08_smap_delta_rho_windowed.pdf` : State-dependence gain Δρ(t)
- `09_smap_slope_median.pdf` : Local linearization slope median
- `10_smap_slope_sd.pdf` : Local linearization slope dispersion

### Results (CSV)
- `results/summary.csv` : single-source-of-truth for the run
- `results/sensitivity_summary.csv` : robustness grid
- `results/transition.csv` : detected t^{\ast} and detection parameters
- `results/embedding_validation_fnn.csv` : FNN(E) table
- `results/session_info.txt` : R environment for reproducibility

---

## Key configuration knobs

All main hyperparameters live in `scripts/run_pipeline.R` under the `config <- list(...)` block.

**Embedding parameters:**
- `config$edm$simplex$E` : Embedding dimension for Simplex (validated via FNN)
- `config$edm$smap$E` : Embedding dimension for S-map
- `config$edm$*$tau` : Time delay for reconstruction

**Forecast horizons:**
- `config$edm$simplex$Tp` : Vector of prediction horizons (e.g., c(1, 10))
- `config$edm$smap$Tp` : Single prediction horizon for S-map

**Windowing:**
- `config$edm$windowed$window` : Analysis window size (tested: 600, 800, 1000)
- `config$edm$windowed$step` : Stride between windows
- `config$edm$windowed$theiler` : Temporal exclusion radius (prevents autocorrelation)
- `config$edm$windowed$lib_frac` : Fraction of window used as library

**Transition detection:**
- `config$detect$w` : Rolling window size for variance proxy
- `config$detect$k` : Threshold multiplier (× baseline variability)
- `config$detect$baseline_max_time` : Duration of baseline regime

Adjusting these parameters allows exploration of sensitivity and robustness.

---

## How to interpret the outputs

This project uses a **controlled, noise-free** ODE simulation with smooth parameter drift. As a result:

**Predictability degradation:**
- Longer forecast horizons (h=10) show earlier and stronger degradation than one-step (h=1)
- This is the **clearest destabilization signature** in clean deterministic systems
- After transition, skill partially recovers as system settles onto oscillatory attractor

**Local linearization (S-map slopes):**
- Median slope **b(t)** approaching |b| ≈ 1 indicates weakening contraction
- Increasing dispersion reflects heterogeneity in local dynamics
- These are **supporting diagnostics** rather than primary indicators

**State dependence (θ\* and Δρ):**
- Can be **weak** in this setting because global linear models (θ=0) may already perform well
- Increases during destabilization but less dramatic than forecast skill changes
- Use as **confirmatory evidence** alongside other indicators

**Important caveats:**
- Results are specific to this tri-trophic system and forcing regime
- Threshold values are heuristic but explicitly documented
- Sensitivity analysis demonstrates robustness across parameter ranges
- Real-world applications would require noise robustness testing

See `docs/interpretation_guide.md` for detailed guidance on result interpretation.

---

## Reproducibility notes

Every run generates:
- `results/session_info.txt` with R version, platform, and package versions
- Complete CSV outputs that can be independently analyzed
- Publication-ready figures with consistent styling

All results are **regenerable** from source. The pipeline has no hidden state or manual steps.

To verify reproducibility:
1. Clone repository
2. Install `deSolve`
3. Run `Rscript scripts/run_pipeline.R`
4. Compare your `results/summary.csv` with included version

---

## Project scope and limitations

**What this project demonstrates:**
- Equation-free early warning under partial observability
- Multiple coordinated EDM diagnostics
- Sensitivity analysis across detection parameters
- Reproducible research pipeline

**What this project does not claim:**
- Universal threshold values for all systems
- Robustness to strong observational noise (not tested)
- Real-time operational deployment (batch processing only)
- Performance beyond n ≈ 50k points without optimization

For implementation details, computational complexity, and extension points, see `TECHNICAL_DETAILS.md`.

---

## References

**EDM Methodology:**
- Sugihara & May (1990). "Nonlinear forecasting as a way of distinguishing chaos from measurement error in time series." *Nature*, 344:734–741.
- Sugihara (1994). "Nonlinear forecasting for the classification of natural time series." *Phil. Trans. R. Soc. A*, 348:477–495.
- Sugihara et al. (2012). "Detecting causality in complex ecosystems." *Science*, 338:496–500.
- Ye & Sugihara (2016). "Information leverage in interconnected ecosystems." *Science*, 353:922–925.

**Early Warning Signals:**
- Dakos et al. (2012). "Methods for detecting early warnings of critical transitions in time series illustrated using simulated ecological data." *PLoS ONE*, 7:e41010.

**Hastings-Powell Model:**
- Hastings & Powell (1991). "Chaos in a three-species food chain." *Ecology*, 72(3):896–903.

---

## Citation

If you use this code or methodology, please cite:

```bibtex
@software{najar2025_edm_early_warning,
  author = {Najar, Mohamed Oussama},
  title = {Early Warning Signals in Tri-Trophic Dynamics via Empirical Dynamic Modeling},
  year = {2025},
  url = {https://github.com/OussamaNajar/tritrophic-early-warning}
}
```

---

## License

MIT License - see `LICENSE` file for details.

---

## Contact

**Author:** Mohamed Oussama Najar  
**Institution:** UC Santa Cruz, Department of Applied Mathematics

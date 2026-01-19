# Early Warning Signal Interpretation Guide

This document provides guidance for interpreting the outputs of the
**tri-trophic early warning pipeline** based on Empirical Dynamic Modeling (EDM).
It is intended to help users correctly interpret diagnostic signals,
understand system-specific behavior, and avoid common misinterpretations.

This guide complements the manuscript and repository documentation and
should be read alongside the figures in `figures/` and summary outputs in `results/`.

---

## Scope and Interpretation Philosophy

This project studies early warning of **loss of equilibrium stability**
in a **deterministic, smoothly forced nonlinear system** using only a **single observed variable**.

Key principles:

- Diagnostics are **empirical**, not mechanistic.
- The rolling-variance transition time \(t^\*\) is used **only as a temporal benchmark**, not as the EDM signal itself.
- EDM indicators are evaluated for **relative timing and coordination**, not universal thresholds.
- Weak or null signals are **scientifically meaningful**, not failures.

---

## System Behavior Overview

Under slow parameter drift, the Hastings–Powell tri-trophic system exhibits:

- A long quasi-stationary equilibrium phase
- A destabilization interval with increasing variability and sensitivity
- A transition to sustained oscillations (limit-cycle dynamics)
- Partial recovery of predictability after the new attractor forms

This structure explains why different diagnostics behave differently over time.

---

## Primary EDM Diagnostics

### 1. Forecast Skill Degradation (Simplex)

**What to examine:**
- Changes in forecast skill *prior to* \(t^\*\)
- Stronger effects at longer horizons (e.g. \(h=10\) vs \(h=1\))
- Partial recovery after oscillations emerge

**Interpretation:**
- Long-horizon forecasts are sensitive to weak local expansion rates.
- Earlier degradation at \(h=10\) indicates weakening stability before the transition.
- Recovery after the transition reflects geometric regularity of the limit cycle.

**Key figure:**  
`06_simplex_skill_windowed.pdf`

---

### 2. Local Linearization Behavior (S-map Slopes)

**What to examine:**
- Median slope approaching \(|b(t)| \approx 1\)
- Increasing within-window dispersion of slopes
- Changes concentrated during the destabilization interval

**Interpretation:**
- The S-map coefficient \(b(t)\) is a **local sensitivity proxy** in reconstructed space.
- Values approaching \(|b| \approx 1\) indicate weakened contraction.
- Increased dispersion reflects growing heterogeneity in local linear behavior.
- These slopes are **not** Jacobian eigenvalues and should not be interpreted mechanistically.

**Key figures:**
- `09_smap_slope_median.pdf`
- `10_smap_slope_sd.pdf`

---

### 3. State Dependence (S-map θ and Δρ)

**What to examine:**
- Whether localized models outperform global linear fits
- Temporal changes in θ\* or Δρ near destabilization

**Observed behavior in this system:**
- θ\* remains near zero for much of the record
- Δρ is small or inconsistent

**Interpretation:**
- Global linear models already perform well in this smooth, deterministic system.
- Weak state-dependence signals are expected and **do not invalidate** other EDM diagnostics.
- Early warning does **not require strong nonlinearity**.

**Key figures:**
- `07_smap_theta_opt_windowed.pdf`
- `08_smap_delta_rho_windowed.pdf`

---

## Role of Variance-Based Transition Detection

A rolling standard deviation of \(Z(t)\) is used to define an empirical transition time \(t^\*\).

**Important clarifications:**
- This is **not** the primary EDM result.
- It provides a **benchmark for temporal alignment**, not validation.
- EDM diagnostics are interpreted relative to \(t^\*\), not derived from it.

**Key figure:**  
`03_transition_rollsd.png`

---

## Summary Panel Interpretation

The summary panel (`11_summary_panel.png`) should be interpreted jointly:

- Forecast degradation appears during destabilization
- S-map slope behavior changes in the same interval
- State-dependence remains weak
- Variance rises sharply near the transition

**Key insight:**  
No single indicator dominates; the signal lies in their **coordination in time**.

---

## Common Misinterpretations (Avoid These)

### “Forecast skill is too high, so the method failed”
Incorrect. High baseline predictability reflects good state-space reconstruction.
The **timing and horizon-dependence** of changes carry the signal.

---

### “θ\* doesn’t increase, so the system isn’t nonlinear”
Incorrect. θ\* measures **predictive benefit of localization**, not intrinsic nonlinearity.
This system weakens stability without developing strong localized nonlinear structure.

---

### “Variance detection is the only thing that works”
Incorrect. Variance provides a useful benchmark, but EDM diagnostics
anticipate and contextualize the transition rather than merely detecting it.

---

## Kendall Tau Trend Interpretation

Kendall τ statistics quantify **monotonic trends**, not effect size.

In this system:
- τ values may be weak or mixed
- This reflects smooth forcing and gradual destabilization
- Visual alignment and multi-indicator consistency are more informative

Null or weak τ values should be **reported honestly**, not adjusted.

---

## Recommended Reporting Language

> “Early warning signals emerged during the destabilization interval preceding the empirical transition. Long-horizon forecast skill degraded earlier than short-horizon skill, and S-map slope behavior indicated weakening local contraction and increasing heterogeneity. State-dependence measures were comparatively weak, consistent with smoothly forced dynamics in which global linear models retain high predictive skill.”

---

## When Results Differ

Differences may arise due to:
- Forcing rate
- Window size
- Observation noise
- Partial observability
- Transition definition

Diagnostics must always be interpreted **relative to system context**, not theoretical expectation.

---

## Final Validation Checklist

Before final interpretation:

- [ ] \(t^\*\) detected in a plausible range
- [ ] Forecast skill shows temporal structure
- [ ] S-map slopes change during destabilization
- [ ] Weak state dependence is documented, not hidden
- [ ] Figures and text tell the same story

If all are satisfied → **analysis is internally consistent and scientifically valid**.

---

## Bottom Line

This project demonstrates that:

- Early warning can arise from **predictability changes** even when nonlinearity is weak
- Variance, forecast skill, and local linearization provide **complementary information**
- Honest reporting of weak signals strengthens—not weakens—the analysis

Your pipeline is not failing to match theory.  
It is revealing **which mechanisms matter in this system**.

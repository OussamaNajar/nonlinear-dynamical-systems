# Assumptions and Heuristics

This project is a **controlled methodological demonstration** of early-warning diagnostics using Empirical Dynamic Modeling (EDM).  
It is not intended as a universal, parameter-free early-warning detector.

All assumptions, thresholds, and heuristics are made **explicit** to ensure transparency, reproducibility, and correct interpretation.

---

## Scope of Validity

The results and interpretations in this repository are valid under the following conditions:

- **Deterministic or weakly stochastic dynamics**
- **Smooth, quasi-static parameter forcing**
- **Uniform sampling in time**
- **Partial observability** (single observed variable)
- **Sufficient data density** to populate reconstructed state space

Outside this regime (e.g., strong noise, sparse sampling, abrupt forcing), performance and interpretation may differ.

---

## Transition Definition

### Operational Definition of Transition Time (`t*`)

The transition time `t*` is defined **operationally**, not theoretically.

- `t*` is detected using a **rolling standard deviation (SD)** of the observed variable
- The rolling SD is computed over a sliding window of width `w`
- A transition is flagged when rolling SD exceeds a baseline threshold

**Threshold formula:**

```
threshold = median(baseline) + k × MAD(baseline)
```

where:
- `baseline` is an early, pre-transition segment of the time series (t ≤ `baseline_max_time`)
- `MAD` is the median absolute deviation
- `k` is a user-defined sensitivity multiplier (default: k = 3)
- `w` is the rolling window width (default: w = 300)

**In this project:**
- `baseline_max_time = 700`
- `k = 3`
- `w = 300`
- Result: `t* = 3905`

This definition serves **only as a temporal reference** for aligning diagnostics.

**It does NOT claim:**
- Detection of a mathematical bifurcation
- Identification of a unique "true" transition point
- Universality across systems

---

## Early-Warning Diagnostics

### Advance-Warning Time

An "advance warning time" is defined as:

```
advance time = t* − t_warning
```

where `t_warning` is the first time an indicator crosses a predefined heuristic threshold before `t*`.

Advance times are reported **relative to the same t* definition** to ensure internal consistency.

---

### Heuristic Warning Thresholds

The following warning thresholds are **explicit heuristics**, not theoretically derived constants.

They are used to enable comparison across runs and sensitivity analysis.

#### θ_opt (State-Dependence)

**Warning is triggered when:**

```
θ_opt > 2
```

for the first time prior to `t*`.

**Interpretation:**
- Indicates emergence of locally state-dependent structure
- Measures when localized nonlinear models become beneficial
- **In this system:** θ_opt often remains ≈ 0 (near-linear dynamics)
- **Absence of a θ_opt warning is expected and meaningful**

---

#### Δρ (Predictive Gain)

**Warning is triggered when:**

```
Δρ > 99th percentile of early baseline Δρ
```

**Interpretation:**
- Measures whether localized (nonlinear) models outperform global linear ones
- Δρ = ρ(θ_opt) − ρ(θ=0)
- **In highly predictable or near-linear systems:** Δρ may remain near zero or negative
- This indicator is therefore **system-dependent**

---

## What These Heuristics Are — and Are Not

### These heuristics ARE:

- Explicit and reproducible
- Tunable by the user
- Designed for sensitivity testing
- Appropriate for controlled demonstrations

### These heuristics are NOT:

- Universal thresholds
- Optimal across all systems
- Claims of statistical optimality
- Substitutes for mechanistic stability analysis

**Users are encouraged to:**

1. Adjust `w`, `k`, and warning thresholds
2. Examine robustness using `results/sensitivity_summary.csv`
3. Report how conclusions change under parameter variation

---

## Interpretation Philosophy

This project follows a **data-honest interpretation strategy**:

- Indicators are evaluated based on **observed behavior**, not theoretical expectation
- Weak or absent signals are **reported transparently**
- No attempt is made to "force" indicators to match canonical early-warning narratives

### System-Specific Behavior

In particular, this tri-trophic system exhibits:

- **Near-perfect baseline predictability** (ρ > 0.99)
  - Limits forecast-degradation signals
  - High baseline = less room for observable decline

- **Near-linear dynamics** (θ_opt ≈ 0)
  - Global linear models already perform optimally
  - Limits state-dependence measures
  - Expected behavior for smoothly forced deterministic systems

- **Variance-based indicators dominate**
  - Rolling SD provides 1256 timestep advance warning
  - Primary reliable indicator in this regime

**This behavior is scientifically informative, not a failure of the method.**

---

## Embedding Parameters

### Embedding Dimension (E)

**Choice:** E = 3

**Validation:** False Nearest Neighbors (FNN) analysis
- FNN drops from ~7% (E=2) to <0.1% (E=3)
- Remains near zero for E > 3
- E = 3 is sufficient to unfold the attractor

**Held constant** across all analyses to attribute temporal changes to forcing, not embedding re-optimization.

---

### Time Delay (τ)

**Choice:** τ = 1 (native sampling rate)

**Justification:**
- System is densely sampled (Δt = 1)
- No evidence of undersampling
- Standard choice for regularly sampled data

---

## Windowed Analysis Parameters

### EDM Analysis Window

**Choice:** W = 800 timesteps

**Tested:** W ∈ {600, 800, 1000}

**Result:** Qualitative patterns robust across all window sizes (see `sensitivity_summary.csv`)

**Trade-off:**
- Larger W → smoother trends, more stable estimates, delayed signals
- Smaller W → higher temporal resolution, noisier estimates

---

### Step Size

**Choice:** step = 50 timesteps

**Result:** 120 windows for n = 6000 timesteps

**Purpose:** Balance temporal resolution with computational cost

---

### Theiler Window

**Choice:** theiler = 50 timesteps

**Purpose:** Exclude temporally adjacent neighbors in forecast skill estimation

**Justification:**
- Prevents autocorrelation artifacts
- Must be > max(τ × E, autocorrelation_length)
- Standard practice in EDM applications

---

## Assumptions About the System

### What We Assume:

1. **Stationarity of dynamics** (except for slow parameter drift)
   - No sudden regime changes
   - No structural breaks
   - Smooth forcing only

2. **Observability**
   - Single variable (Z) contains sufficient information
   - E = 3 dimensional reconstruction captures essential dynamics

3. **Determinism**
   - No observational noise
   - No process noise
   - Dynamics fully governed by equations

### What We Do NOT Assume:

1. **Knowledge of governing equations**
   - EDM is equation-free
   - No mechanistic model required

2. **Linearity**
   - Methods work for nonlinear systems
   - But also reveal when systems remain linear

3. **Strong nonlinearity**
   - Early warning does not require chaotic dynamics
   - Linear destabilization can produce signals

---

## Parameter Sensitivity

All critical parameters have been tested for sensitivity:

**Detection parameters:**
- `w` ∈ {200, 300}
- `k` ∈ {2, 3, 4} (implicit in threshold choice)

**EDM parameters:**
- `W` ∈ {600, 800, 1000}

**Results:**
- Qualitative conclusions robust
- Quantitative advance times vary
- See `results/sensitivity_summary.csv` for details

**Recommendation:** Always report sensitivity analysis alongside point estimates.

---

## Known Limitations

### 1. Computational Cost
- O(n²) scaling limits practical applicability to n ≈ 50,000
- Larger datasets require downsampling or optimization

### 2. Heuristic Thresholds
- θ_opt > 2 and 99th percentile cutoffs are not universal
- Require domain-specific tuning

### 3. Single-Variable Observability
- Assumes Z(t) contains sufficient information
- Multivariate extensions may improve performance

### 4. Deterministic Focus
- Noise robustness not extensively tested
- Real-world applications will have measurement error

### 5. Smooth Forcing Assumption
- Abrupt parameter changes may violate quasi-static approximation
- Sensitivity to forcing rate not systematically explored

---

## Reproducibility Checklist

To reproduce or extend this work:

- [ ] Use same random seed (if stochastic elements added)
- [ ] Verify E = 3, τ = 1, theiler = 50
- [ ] Use same detection parameters (w = 300, k = 3, baseline ≤ 700)
- [ ] Apply same warning thresholds (θ_opt > 2, Δρ > 99th percentile)
- [ ] Run sensitivity analysis over parameter grids
- [ ] Report both point estimates and ranges

All parameters are documented in `scripts/run_pipeline.R`.

---

## Summary

**Transition detection is:**
- Heuristic and explicit
- Based on rolling variance
- Relative, not absolute

**Advance-warning metrics are:**
- Defined relative to t*
- Threshold-dependent
- Adjustable and transparent

**Thresholds are:**
- Explicit heuristics
- Not universal constants
- Designed for sensitivity testing

**Interpretation prioritizes:**
- Honesty over conformity
- Observed behavior over theoretical expectation
- Transparency over optimization

**Users should treat this repository as:**
- A reference implementation
- An analysis template
- A demonstration of methodology

**NOT as:**
- A black-box detector
- A universal early-warning system
- A claim of parameter-free detection

---

## For Further Reading

**On threshold selection:**
- Dakos et al. (2012) - Early warning indicator methods
- Scheffer et al. (2009) - Critical transitions framework

**On EDM methodology:**
- Sugihara & May (1990) - Simplex forecasting
- Sugihara et al. (2012) - S-map and causality detection
- Ye et al. (2015) - Equation-free forecasting

**On interpretation:**
- See `docs/interpretation_guide.md` for detailed guidance
- See `TECHNICAL_DETAILS.md` for implementation specifics
- See `results/sensitivity_summary.csv` for robustness assessment

---

## Contact

For questions about assumptions, thresholds, or interpretation:
- Open an issue: github.com/OussamaNajar/tritrophic-early-warning/issues
- Refer to: `docs/interpretation_guide.md`

Contributions and extensions welcome via pull requests.

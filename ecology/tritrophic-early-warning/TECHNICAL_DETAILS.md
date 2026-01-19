# Technical Details

This document provides implementation details, algorithmic complexity, parameter sensitivity, and extension points for the EDM early warning system.

---

## Implementation Details

### Core Algorithms

#### 1. Delay-Coordinate Embedding

**Algorithm:** Takens' embedding theorem reconstruction

```
Input: Time series x(t), embedding dimension E, time delay τ
Output: Reconstructed state vectors X(t) ∈ R^E

X(t) = [x(t), x(t-τ), x(t-2τ), ..., x(t-(E-1)τ)]
```

**Complexity:** O(n) where n = length of time series  
**Memory:** O(n × E) for embedding matrix  
**Implementation:** `R/edm_simplex.R::make_embed()`

**Parameter Selection:**
- **E (embedding dimension):** Validated via False Nearest Neighbors (FNN)
  - E=3 chosen based on FNN < 5% threshold
  - Computational cost scales linearly with E
- **τ (time delay):** Set to 1 (native sampling rate)
  - For undersampled data, use τ = ceiling(1/sampling_rate)
  - Rule of thumb: τ ≈ first zero-crossing of autocorrelation

#### 2. Simplex Forecasting

**Algorithm:** k-nearest neighbor averaging in reconstructed space

```
For each point x(t):
  1. Find k = E+1 nearest neighbors in embedding space
  2. Compute distances d_i using Euclidean metric
  3. Weight by w_i = exp(-d_i / d_min)
  4. Predict: x(t+h) = Σ w_i · x_i(t+h) / Σ w_i
```

**Complexity:** O(n² × E) per forecast horizon  
**Memory:** O(n × E) for distance computation  
**Implementation:** `R/edm_simplex.R::simplex_forecast()`

**Optimizations:**
- **Theiler window:** Excludes temporally adjacent points (50 timesteps)
  - Prevents autocorrelation artifacts
  - Reduces effective n by ~50
- **Exponential weighting:** Emphasizes nearest neighbors
  - More stable than uniform k-NN
  - Less sensitive to k choice

#### 3. S-map Local Linearization

**Algorithm:** Locally weighted linear regression

```
For each point x(t):
  1. Weight all points by w_i = exp(-θ × d_i / d_mean)
  2. Solve: B = argmin_B Σ w_i (x(t+1) - B·X(t))²
  3. Extract coefficient b = ∂x(t+1)/∂x(t)
```

**Complexity:** O(n² × E²) with ridge regression  
**Memory:** O(E²) for design matrix inversion  
**Implementation:** `R/edm_smap.R`

**Parameters:**
- **θ (nonlinearity):** Controls localization strength
  - θ = 0: Global linear model
  - θ > 0: Local models (higher = more local)
  - Swept over [0, 8] in increments of 0.5
- **Ridge regularization:** λ = 10⁻⁸
  - Prevents singular matrices near attractors
  - Adjustable for different system scales

---

## Computational Performance

### Benchmark (macOS Sonoma, M-series)

| Operation | n=6000 | Complexity | Time |
|-----------|--------|------------|------|
| **Embedding** | E=3 | O(n) | <0.1s |
| **Simplex (h=1)** | Full | O(n²×E) | ~15s |
| **Simplex (h=10)** | Full | O(n²×E) | ~15s |
| **S-map θ sweep** | θ∈[0,8] | O(17×n²×E²) | ~45s |
| **Windowed analysis** | W=800, step=50 | O(120×n×E²) | ~120s |
| **Full pipeline** | All | Combined | ~200s |

### Scalability

**Time series length:**
- Linear scaling for embedding
- Quadratic scaling for EDM (distance computation)
- Practical limit: n ≈ 50,000 on desktop
- For n > 50k: Consider downsampling or algorithmic improvements

**Embedding dimension:**
- Linear memory overhead
- Quadratic computation overhead (E² in S-map)
- Practical limit: E ≤ 10 for reasonable performance

**Window size:**
- Affects number of windows: O(n / step)
- Larger windows = more stable estimates, fewer samples
- Recommended: W > 10 × E × (E+1)

---

## Parameter Sensitivity

### Critical Parameters

#### Transition Detection (rolling variance)
```r
detect = list(
  w = 300,              # Window size (timesteps)
  k = 3,                # Threshold multiplier (× baseline SD)
  baseline_max_time = 700  # Baseline duration
)
```

**Sensitivity:**
- **w:** Smaller = earlier detection, more false positives
  - Tested: [200, 300] (see `sensitivity_summary.csv`)
  - Recommendation: w ≈ 5-10% of total series length
- **k:** Lower = more sensitive, earlier warning
  - k=2: Aggressive (more false alarms)
  - k=3: Balanced (used in baseline)
  - k=4: Conservative (later detection)

#### EDM Windowed Analysis
```r
edm = list(
  windowed = list(
    window = 800,       # Analysis window size
    step = 50,          # Stride between windows
    theiler = 50        # Temporal exclusion radius
  )
)
```

**Sensitivity:**
- **window:** Tested [600, 800, 1000]
  - Results robust across all three (see sensitivity analysis)
  - Larger = smoother trends, delayed signal
- **step:** Affects temporal resolution
  - step=50: 120 windows for n=6000
  - Smaller step = higher resolution, more computation
- **theiler:** Prevents autocorrelation leakage
  - Must be > max(τ × E, autocorrelation_length)
  - Too small: Spurious predictability
  - Too large: Insufficient data

---

## Extension Points

### 1. Noise Robustness

**Current:** Deterministic dynamics only  
**Extension:** Add observational or process noise

```r
# In model_hp.R, modify output:
dZ + noise_amplitude * rnorm(1)
```

**Expected impact:**
- Forecast skill will degrade
- Early warning lead time may reduce
- Recommend increasing window size proportionally

### 2. Multivariate Observables

**Current:** Single scalar (Z only)  
**Extension:** Observe multiple variables

```r
# Create multivariate embedding:
X_multi = cbind(
  make_embed(X, E=3, tau=1),
  make_embed(Y, E=3, tau=1)
)
```

**Benefits:**
- Improved reconstruction quality
- Earlier warning signals
- More robust to partial observability

### 3. Real-Time Forecasting

**Current:** Batch processing  
**Extension:** Online/streaming implementation

```r
# Pseudocode for online mode:
maintain_sliding_buffer(size = window + (E-1)*tau)
for new_point in stream:
  update_buffer(new_point)
  if buffer_full:
    compute_metrics()
    check_thresholds()
```

**Challenges:**
- Computational latency (must process faster than data rate)
- Concept drift (system properties changing)
- False alarm rate control

### 4. Adaptive Thresholding

**Current:** Fixed percentile thresholds  
**Extension:** Dynamic threshold adjustment

```r
# Example: CUSUM-like accumulation
accumulator = 0
for t in time:
  signal = compute_metric(t)
  accumulator = max(0, accumulator + (signal - baseline_mean))
  if accumulator > adaptive_threshold:
    trigger_warning()
```

**Benefits:**
- Reduced false alarm rate
- Better performance under nonstationarity

---

## Validation Checklist

Before claiming early warning performance on new systems:

- [ ] **Embedding validation:** FNN analysis shows clear minimum
- [ ] **Baseline stability:** Forecast skill ρ > 0.8 in pre-transition regime
- [ ] **Sensitivity analysis:** Results hold across 2+ window sizes
- [ ] **Visual inspection:** Windowed plots show coordinated trends
- [ ] **Statistical significance:** Kendall τ p-values < 0.05 (if needed)
- [ ] **Transition verification:** Independent indicator confirms t*
- [ ] **Parameter documentation:** All thresholds explicitly reported
- [ ] **Negative controls:** Verify no false warnings in stable regime

---

## Known Limitations

1. **Computational cost:** O(n²) limits scalability to n ≈ 50k
2. **Stationarity assumption:** Requires consistent sampling rate
3. **Partial observability:** Single variable may miss some transitions
4. **Heuristic thresholds:** Require domain-specific tuning
5. **Deterministic focus:** Noise robustness not extensively tested

---

## References & Further Reading

**Core EDM Theory:**
- Takens (1981): "Detecting strange attractors in turbulence"
- Sugihara & May (1990): "Nonlinear forecasting as a way of distinguishing chaos from measurement error"

**Early Warning Applications:**
- Dakos et al. (2008): "Slowing down as an early warning signal"
- Scheffer et al. (2009): "Early-warning signals for critical transitions"
- Sugihara et al. (2012): "Detecting causality in complex ecosystems"

**Implementation:**
- rEDM package (C++ implementation): github.com/SugiharaLab/rEDM
- EDM tutorial: sugihara.bio.ucsc.edu

---

## Contact & Contributions

For questions, bug reports, or extensions:
- Open an issue: github.com/OussamaNajar/tritrophic-early-warning/issues

Contributions welcome via pull requests.

# Nonlinear Dynamical Systems

Research-grade computational studies of **nonlinear dynamical systems** focused on:
- **stability & attractor structure**
- **predictability limits** (sensitivity to initial conditions / noise)
- **regime behavior & transitions**
- **model-based + data-driven analysis**

This repository is **methods-first**: each study is organized around transferable tooling
(e.g., simulation, parameter sweeps, stability diagnostics, regime detection), with
applications spanning ecology and finance.

---

## Why this matters

Many real systems—ecosystems, markets, climate, engineered control loops—share the same hard problems:
- nonlinear feedback creates multi-regime behavior
- prediction skill collapses beyond a horizon (chaos / noise)
- distribution shift/regime change dominates error
- “better fitting” is not the same as “more predictable”

This repo is built to answer: **what is structurally predictable, what is not, and why**.

---

## Structure

- `ecology/`  
  Ecological population dynamics and tri-trophic food-chain models (oscillations, chaos, bifurcation-like behavior).

- `finance/` 
  Regime-aware time series studies and predictability diagnostics for market-like systems (structure-first, not “alpha hype”).

Each study is self-contained with:
- `src/` code
- `experiments/` parameter sweeps/configs
- `figures/` generated plots
- `paper/` report/manuscript-quality writeup
- `README.md` explaining results, assumptions, and how to reproduce

---

## Featured study

### Tri-Trophic Ecological Dynamics (Hastings–Powell family)
Location: `ecology/tri_trophic_dynamics/`

A computational study of a nonlinear tri-trophic food-chain model emphasizing:
- attractor geometry and long-run behavior
- parameter sensitivity and regime transitions
- predictability limits under perturbations
- interpretation for management / decision context

---

## Reproducibility

Each study directory contains its own `README.md` with:
- environment requirements
- how to run simulations
- how to reproduce figures and summary results

---

## Context

Parts of this work were developed during graduate training in scientific computing and ecological modeling.
All artifacts included here are curated to be **standalone, reproducible, and industry-readable**.

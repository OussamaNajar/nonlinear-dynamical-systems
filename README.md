# Nonlinear Dynamical Systems

Research-grade computational studies of **nonlinear dynamical systems** focused on:
- **stability & attractor structure**
- **predictability limits** (sensitivity to initial conditions / noise)
- **regime behavior & transitions**
- **model-based + data-driven analysis**

This repository is **methods-first**: each study is organized around transferable tooling
(e.g., simulation, parameter sweeps, stability diagnostics, regime detection), with extensions to other domains (e.g. finance, control, climate).

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
  Ecological population dynamics and tri-trophic food-chain models (oscillations, chaos, regime shifts).

Each study is self-contained with:
- `R/` or `src/` code
- `scripts/` runnable entrypoints
- `results/` generated outputs
- `figures/` generated plots
- `paper/` manuscript-quality writeup
- `README.md` explaining results, assumptions, and reproduction steps

---

## Featured study

### Tri-Trophic Ecological Dynamics (Hastings–Powell family)
Location: `ecology/tritrophic-early-warning/`

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




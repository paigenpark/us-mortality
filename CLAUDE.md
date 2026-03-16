# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Research project analyzing US life expectancy at birth (e0) by state (1959-2020) using structural time series (STS) / state-space models. The core model is a Random Walk with Drift decomposing mortality shocks into permanent and transitory components.

## Environment Setup

```bash
# Install pixi if needed: curl -fsSL https://pixi.sh/install.sh | sh
pixi install
pixi run setup   # runs renv::restore() to install R packages
```

R >=4.5.2 is managed by pixi. R packages are managed by renv (renv.lock).

## Running Code

```bash
pixi run Rscript code/1-data-prep.R    # Prepare e0 data from life tables
pixi run Rscript code/2-modeling.R      # Fit STS models (main analysis)
pixi run Rscript code/3-spatial_analysis.R  # Moran's I spatial autocorrelation
pixi run Rscript code/eda.R             # Exploratory analysis & clustering
```

## Architecture

**Pipeline**: `1-data-prep.R` → `2-modeling.R` → `3-spatial_analysis.R`, with `eda.R` as standalone exploration.

**Modeling backends** (in `2-modeling.R`):
- `StructTS` — R built-in; tends to collapse variance to zero
- `MARSS/KEM` — EM algorithm via MARSS package
- `MARSS/BFGS` — quasi-Newton optimizer via MARSS; preferred for stable convergence
- `KFAS` — Kalman filter allowing known sampling variance as fixed observation error

**State-space model**:
- State: μ_t = μ_{t-1} + d + η_t, η_t ~ N(0, q) — level with drift
- Observation: y_t = μ_t + ε_t, ε_t ~ N(0, r) — with measurement error

**Key scripts beyond the pipeline**:
- `kfas_for_paige.R` — KFAS with known sampling variance, simulation studies
- `paige_table_1_unpooled_and_pooled_subset_v2.R` — unpooled vs pooled MARSS on 8-state subset, tests variance constraint ratios
- `nowcasting.R` — real-time mortality estimation with incomplete data
- `simulation.R` — variance estimation accuracy tests

**Data**: Life table CSVs in `data/` (50 states, 9 Census divisions, national). `sampling_error.txt` contains known sampling standard deviations by state.

**Results**: PDF plots and tables output to `results/`.

### =========================================================================
### 2-modeling.R
###
### Purpose: Fit structural time series (STS) models to US state-level life
###          expectancy data (e0). The core model is a random walk with drift
###          (RWD), optionally augmented with observation error. We compare
###          four estimation backends:
###            1. StructTS  -- R's built-in state-space optimizer (least flexible)
###            2. MARSS/KEM -- EM algorithm via the MARSS package
###            3. MARSS/BFGS -- quasi-Newton optimizer via the MARSS package
###            4. KFAS      -- Kalman filter with known sampling variance
###
###          Model results are saved to results/ for downstream table/figure
###          generation in separate .qmd files.
### =========================================================================

### ---- Packages & functions ----
library(tidyverse)
library(here)
library(MARSS)
library(KFAS)

source(here("code", "model-functions.R"))


### =========================================================================
### Load data
### =========================================================================

state_e0 = read.csv(here("data", "combined_e0.csv"), header = TRUE, row.names = 1)

# US aggregate (first row of combined_e0)
us_e0 = state_e0[1, ]
state_e0 = state_e0[-1, ]  # remaining rows are states


### =========================================================================
### Fit models to all 50 states
### =========================================================================

# --- MARSS with EM (KEM) ---
# 6 states fail to converge (8 if 2020 is excluded)
state_kem = as.data.frame(apply(state_e0, 1, rwd_with_obs_error_kem))

# --- MARSS with BFGS ---
# 6 states have near-zero observation variance
state_bfgs = as.data.frame(apply(state_e0, 1, rwd_with_obs_error_bfgs))

# --- StructTS ---
# 13 states have exactly zero obs variance (10 if 2020 is excluded)
state_sts = as.data.frame(apply(state_e0, 1, rwd_with_obs_error))

# --- KFAS (with known sampling variance from bootstrap) ---
boot_results = read.csv(here("data", "boot_results_v2.csv"), header = TRUE, row.names = 1)
state_kfas = as.data.frame(t(sapply(1:nrow(state_e0), function(i) {
    run_kfas(state_e0[i, ], known_sampling_var = boot_results[i, "e0.var"])
  })))
rownames(state_kfas) <- rownames(state_e0)

# --- Full MARSS BFGS objects (for fitted values, AIC, CIs) ---
set.seed(10)
marss_sts_fits <- apply(state_e0, 1, rwd_obs_error_bfgs_out)
marss_rwd_fits <- apply(state_e0, 1, rwd_bfgs)

# US aggregate fits
marss_sts_us <- rwd_obs_error_bfgs_out(us_e0)
marss_rwd_us <- rwd_bfgs(us_e0)

# --- AIC/AICc from MARSS BFGS ---
aic_aicc <- data.frame(
  rwd_aic  = sapply(marss_rwd_fits, function(x) x$AIC),
  sts_aic  = sapply(marss_sts_fits, function(x) x$AIC),
  rwd_aicc = sapply(marss_rwd_fits, function(x) x$AICc),
  sts_aicc = sapply(marss_sts_fits, function(x) x$AICc)
)
rownames(aic_aicc) <- rownames(state_e0)

# --- MARSS parameter CIs ---
bfgs_with_cis <- lapply(marss_sts_fits, MARSSparamCIs)


### =========================================================================
### Identify states with zero / missing observation variance
### =========================================================================

# StructTS: exact zeros
state_no_obs_sts = colnames(state_sts)[which(state_sts[3,] == 0)]
cat("StructTS zero obs var:", state_no_obs_sts, "\n")

# KEM: NAs indicate non-convergence
state_no_obs_kem = colnames(state_kem)[which(is.na(state_kem[3,]))]
cat("KEM non-converged:", state_no_obs_kem, "\n")

# BFGS: near-zero (< 1e-10) treated as effectively zero
state_no_obs_bfgs = colnames(state_bfgs)[which(state_bfgs[3,] <= 1e-10)]
cat("BFGS near-zero obs var:", state_no_obs_bfgs, "\n")


### =========================================================================
### Save all results for downstream .qmd files
### =========================================================================

saveRDS(list(
  state_e0       = state_e0,
  us_e0          = us_e0,
  state_kem      = state_kem,
  state_bfgs     = state_bfgs,
  state_sts      = state_sts,
  state_kfas     = state_kfas,
  aic_aicc       = aic_aicc,
  marss_sts_fits = marss_sts_fits,
  marss_rwd_fits = marss_rwd_fits,
  marss_sts_us   = marss_sts_us,
  marss_rwd_us   = marss_rwd_us,
  bfgs_with_cis  = bfgs_with_cis,
  boot_results   = boot_results,
  state_no_obs_sts  = state_no_obs_sts,
  state_no_obs_kem  = state_no_obs_kem,
  state_no_obs_bfgs = state_no_obs_bfgs
), here("results", "model_fits.rds"))

write.csv(state_kfas, here("results", "kfas_results.csv"))

cat("Done. Results saved to results/model_fits.rds\n")

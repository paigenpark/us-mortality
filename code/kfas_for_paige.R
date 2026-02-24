## kfas_for_paige.R

## this program uses KFAS to estimate sts **with sampling error**
## it starts by simulating what we want, and then it fits the model
## I do this once, and then repeat it 20 times for some plotting.

## I use semi-realistic numbers from Wisconsin (squaring SDs in Table 1)
# True Parameters
true_drift <- 0.112            # delta
true_var_level <- 0.266^2        # Q_level (Permanent shocks)
true_var_idiosyncratic <- 0.106^2 # Q_transitory (Transitory shocks)
known_sampling_var <- 0.067^2    # H (Known measurement error)

set.seed(123)
n <- 62 ## about the number of years

# Generate Components
eta <- rnorm(n, mean = 0, sd = sqrt(true_var_level))
nu  <- rnorm(n, mean = 0, sd = sqrt(true_var_idiosyncratic))
eps <- rnorm(n, mean = 0, sd = sqrt(known_sampling_var))

# Generate the State (Random walk with drift)
mu <- numeric(n)
mu[1] <- 10 # Initial level
for(t in 2:n) {
  mu[t] <- mu[t-1] + true_drift + eta[t]
}

# Generate Observed Series
# y = State + Transitory Shock + Sampling Error
y <- mu + nu + eps

pdf("wisconsin_simulated_e0.pdf")
par(oma = c(0, 0, 2, 0))  # top outer margin
plot.ts(y, main="Simulated Series", col="gray")
lines(mu, col="red", lwd=2) # The "True" underlying signal
legend("topleft", legend=c("Observed", "True Signal"), col=c("gray", "red"), lty=1)
dev.off()

library(KFAS)

# Define the Model Structure
# SSMtrend(2) gives us Level + Slope. 
# We fix Slope variance to 0 to make it a Constant Drift.
# H is fixed to our known sampling variance.
model_spec <-
    SSModel(y ~ SSMtrend(degree = 2,
                         Q = list(matrix(NA), matrix(0))) + 
                SSMcustom(Z = 1, T = 0, R = 1, Q = NA), 
            H = known_sampling_var)

# Initial guesses for variances (log-space)
# fitSSM optimizes the NA values in the order they appear in the model
initial_guesses <- log(c(var(diff(y)), var(y)/10))

library(KFAS)


# --- 1. Re-run Fit with Hessian ---
# We must set hessian = TRUE to get the uncertainty of the variance estimates
fit <- fitSSM(model_spec, inits = initial_guesses, hessian = TRUE)

# --- 2. Extract Parameter Estimates ---
# Variances (from Q matrix)
est_var_level <- fit$model$Q[1, 1, 1]
est_var_idiosyncratic <- fit$model$Q[3, 3, 1]

# Drift (from the smoothed state 'slope')
out <- KFS(fit$model)
est_drift <- tail(out$alphahat[, "slope"], 1)
# Standard error for drift is directly available from the KFS output covariance matrix
est_drift_se <- sqrt(tail(out$V[2, 2, ], 1)) 

# --- 3. Calculate Confidence Intervals (Delta Method / Log-space) ---
# The Hessian is for the 'NA' parameters in the order they appear in the model
# Par 1: log(Var_Level), Par 2: log(Var_Idiosyncratic)
se_log <- sqrt(diag(solve(fit$optim.out$hessian)))
params_log <- fit$optim.out$par

# CI for Level Variance
ci_level <- exp(params_log[1] + c(-1.96, 1.96) * se_log[1])
# CI for Idiosyncratic Variance
ci_idio  <- exp(params_log[2] + c(-1.96, 1.96) * se_log[2])
# CI for Drift (normal space)
ci_drift <- est_drift + c(-1.96, 1.96) * est_drift_se

# --- 4. Comparison Table ---
results <- data.frame(
  Parameter = c("Drift", "Var_Level", "Var_Idiosyncratic"),
  True = c(true_drift, true_var_level, true_var_idiosyncratic),
  Estimated = c(est_drift, est_var_level, est_var_idiosyncratic),
  LWR = c(ci_drift[1], ci_level[1], ci_idio[1]),
  UPR = c(ci_drift[2], ci_level[2], ci_idio[2])
)

print(results)

# --- 5. Visualization ---

pdf("kfas_results_for_wisconsin.pdf")
library(tinyplot)
tinyplot(
    Estimated ~ Parameter,
    data = results,
    ymin = LWR,
    ymax = UPR,
    type = "pointrange",
    pch = 19,
    col = 'dodgerblue',
    grid = TRUE,
    ylim = c(0, .2),
    main = "estimates")
tinyplot_add(True ~ Parameter,
        type = 'p',
        pch = 1,
        cex = 2,
        data = results,
        col = "red")
dev.off()
##################### below is AI genereated I think it's ok
## but hard to check

## --- 6. Replicate 20 times with Hessian-based CIs in SD units and plot ---
set.seed(789)
n_rep <- 20
rep_results <- data.frame(
  Rep = integer(),
  Parameter = character(),
  Estimate_SD = numeric(),
  LWR_SD = numeric(),
  UPR_SD = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:n_rep) {
  # --- Simulate ---
  eta <- rnorm(n, 0, sqrt(true_var_level))
  nu  <- rnorm(n, 0, sqrt(true_var_idiosyncratic))
  eps <- rnorm(n, 0, sqrt(known_sampling_var))
  
  mu <- numeric(n); mu[1] <- 10
  for(t in 2:n) mu[t] <- mu[t-1] + true_drift + eta[t]
  y <- mu + nu + eps
  
  # --- Fit KFAS model ---
  model_spec <- SSModel(y ~ SSMtrend(degree = 2,
                                     Q = list(matrix(NA), matrix(0))) + 
                            SSMcustom(Z = 1, T = 0, R = 1, Q = NA), 
                        H = known_sampling_var)
  initial_guesses <- log(c(var(diff(y)), var(y)/10))
  fit <- fitSSM(model_spec, inits = initial_guesses, hessian = TRUE)
  out <- KFS(fit$model)
  
  # --- Hessian-based SEs in log-space ---
  se_log <- sqrt(diag(solve(fit$optim.out$hessian)))
  params_log <- fit$optim.out$par
  
  # --- Drift ---
  est_drift <- tail(out$alphahat[, "slope"], 1)
  est_drift_se <- sqrt(tail(out$V[2,2,], 1))
  ci_drift <- est_drift + c(-1.96, 1.96) * est_drift_se
  
  # --- Level variance ---
  est_var_level <- fit$model$Q[1,1,1]
  ci_level <- exp(params_log[1] + c(-1.96,1.96)*se_log[1])
  
  # --- Idiosyncratic variance ---
  est_var_idio <- fit$model$Q[3,3,1]
  ci_idio <- exp(params_log[2] + c(-1.96,1.96)*se_log[2])
  
  # --- Convert to SD units ---
  rep_results <- rbind(rep_results,
                       data.frame(
                         Rep = i,
                         Parameter = "Drift",
                         Estimate_SD = est_drift,
                         LWR_SD = ci_drift[1],
                         UPR_SD = ci_drift[2]
                       ))
  rep_results <- rbind(rep_results,
                       data.frame(
                         Rep = i,
                         Parameter = "Level_SD",
                         Estimate_SD = sqrt(est_var_level),
                         LWR_SD = sqrt(ci_level[1]),
                         UPR_SD = sqrt(ci_level[2])
                       ))
  rep_results <- rbind(rep_results,
                       data.frame(
                         Rep = i,
                         Parameter = "Idio_SD",
                         Estimate_SD = sqrt(est_var_idio),
                         LWR_SD = sqrt(ci_idio[1]),
                         UPR_SD = sqrt(ci_idio[2])
                       ))
}

# --- Plot all 3 parameters separately ---
library(tinyplot)
param_names <- unique(rep_results$Parameter)
true_vals_sd <- c(Drift=true_drift, Level_SD=sqrt(true_var_level), Idio_SD=sqrt(true_var_idiosyncratic))
colors <- c(Drift="dodgerblue", Level_SD="forestgreen", Idio_SD="orange")

pdf("kfas_results_for_wisconsin_20_replicates.pdf")
par(mfrow = c(3,1))
for (p in param_names) {
  df <- subset(rep_results, Parameter == p)
  tinyplot(
    Estimate_SD ~ Rep,
    data = df,
    ymin = LWR_SD,
    ymax = UPR_SD,
    type = "pointrange",
    ylim = c(0, 0.8),
    pch = 19,
    col = colors[p],
    grid = TRUE,
    main = paste0(p, " Estimates Across Replications (SD units)")
  )
  abline(h = true_vals_sd[p], col = "red", lwd = 2)
}
dev.off()


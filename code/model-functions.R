### =========================================================================
### model-functions.R
###
### Wrapper functions for structural time series (STS) model estimation.
### Each function wraps a different backend around the same model:
###   State equation:   mu_t = mu_{t-1} + d + eta_t,   eta_t ~ N(0, q)
###   Obs   equation:   y_t  = mu_t + eps_t,           eps_t ~ N(0, r)
###
### "RWD" sets r = 0 (no observation error).
### "STS" estimates both q and r.
### =========================================================================

library(MARSS)
library(KFAS)

### ---- StructTS wrappers ----

# RWD only: fixes slope variance and obs variance to 0
rwd <- function(y){
  sts = StructTS(x = y, type = "trend", fixed = c(NA, 0, 0))
  d.hat = sts$model$a[2]
  var_innov.hat = sts$coef["level"]
  return(list(d = d.hat,
              var_innov = var_innov.hat))
}

# STS: RWD + observation error
rwd_with_obs_error <- function(y, ..., maxit = 5000){
  sts = StructTS(x = y, type = "trend", fixed = c(NA, 0, NA), ...,
                 optim.control = list(maxit = maxit))
  d.hat = sts$model$a[2]
  var_innov.hat = sts$coef["level"]
  var_obs.hat = sts$coef["epsilon"]
  names(d.hat) <- names(var_innov.hat) <- names(var_obs.hat) <- ""

  return(c(d = d.hat,
           var_innov = var_innov.hat,
           var_obs = var_obs.hat))
}

# RWD with AIC/BIC
rwd_aic_bic <- function(y){
  sts = StructTS(x = y, type = "trend", fixed = c(NA, 0, 0))
  d.hat = sts$model$a[2]
  var_innov.hat = sts$coef["level"]

  aic_val = -2*sts$loglik + 2
  bic_val = -2*sts$loglik + log(length(y))

  return(c(d = d.hat,
           var_innov = var_innov.hat,
           aic = aic_val,
           bic = bic_val))
}

# STS with AIC/BIC
rwd_with_obs_aic_bic <- function(y, ..., maxit = 4000){
  sts = StructTS(x = y, type = "trend", fixed = c(NA, 0, NA), ...,
                 optim.control = list(maxit = maxit))
  d.hat = sts$model$a[2]
  var_innov.hat = sts$coef["level"]
  var_obs.hat = sts$coef["epsilon"]
  names(d.hat) <- names(var_innov.hat) <- names(var_obs.hat) <- ""

  aic_val = -2*sts$loglik + 2*2
  bic_val = -2*sts$loglik + 2*log(length(y))

  return(c(d = d.hat,
           var_innov = var_innov.hat,
           var_obs = var_obs.hat,
           aic = aic_val,
           bic = bic_val))
}

### ---- MARSS wrappers (BFGS optimizer) ----

# STS via BFGS: returns drift (d), innovation variance (q), obs variance (r)
rwd_with_obs_error_bfgs <- function(y, ...)
{
  mod.list.2 <- list(B = matrix(1),
                     U = matrix("d"),
                     Q = matrix("q"),
                     Z = matrix(1),
                     A = matrix(0),
                     R = matrix("r"),
                     x0 = matrix("mu"),
                     tinitx = 0)
  out2.marss = MARSS(y, model = mod.list.2, control=list(maxit=1000), method= "BFGS", ...)

  if (out2.marss$convergence == 0) {
    d = coef(out2.marss)$U[1,1]
    var_innov = coef(out2.marss)$Q[1,1]
    var_obs = coef(out2.marss)$R[1,1]
  } else {
    d = NA
    var_innov = NA
    var_obs = NA
  }

  names(d) <- names(var_innov) <- names(var_obs) <- ""
  return(c(d= d,
           var_innov= var_innov,
           var_obs = var_obs))
}

# RWD via BFGS: returns the full MARSS object
rwd_bfgs <- function(y, ...)
{
  mod.list.2 <- list(B = matrix(1),
                     U = matrix("d"),
                     Q = matrix("q"),
                     Z = matrix(1),
                     A = matrix(0),
                     R = matrix(0),
                     x0 = matrix("mu"),
                     tinitx = 0)
  out2.marss = MARSS(y, model = mod.list.2, control=list(maxit=3000), method= "BFGS", ...)

  return(out2.marss)
}

# STS via BFGS: returns the full MARSS object
rwd_obs_error_bfgs_out <- function(y, ...)
{
  mod.list.2 <- list(B = matrix(1),
                     U = matrix("d"),
                     Q = matrix("q"),
                     Z = matrix(1),
                     A = matrix(0),
                     R = matrix("r"),
                     x0 = matrix("mu"),
                     tinitx = 0)
  out2.marss = MARSS(y, model = mod.list.2, control=list(maxit=3000), method= "BFGS", ...)

  return(out2.marss)
}

# Compare RWD vs STS: extract AIC and AICc from both model fits
get_aic_aicc = function(y, ...) {
  rwd_result = rwd_bfgs(y)
  obserr_result = rwd_obs_error_bfgs_out(y)

  return(c(rwd_result$AIC, obserr_result$AIC, rwd_result$AICc, obserr_result$AICc))
}

### ---- MARSS wrappers (KEM / EM algorithm) ----

# STS via EM: returns drift (d), innovation variance (q), obs variance (r)
rwd_with_obs_error_kem <- function(y, ...)
{
  mod.list.2 <- list(B = matrix(1),
                     U = matrix("d"),
                     Q = matrix("q"),
                     Z = matrix(1),
                     A = matrix(0),
                     R = matrix("r"),
                     x0 = matrix("mu"),
                     tinitx = 0)
  out2.marss = MARSS(y, model = mod.list.2, control=list(maxit=1000), method= "kem", ...)

  if (out2.marss$convergence == 0) {
    d = coef(out2.marss)$U[1,1]
    var_innov = coef(out2.marss)$Q[1,1]
    var_obs = coef(out2.marss)$R[1,1]
  } else {
    d = NA
    var_innov = NA
    var_obs = NA
  }

  names(d) <- names(var_innov) <- names(var_obs) <- ""
  return(c(d= d,
           var_innov= var_innov,
           var_obs = var_obs))
}


### ---- KFAS wrapper ----

run_kfas <- function(y, known_sampling_var, ...) {

  y <- as.numeric(y)

  model_spec <- SSModel(y ~ SSMtrend(degree = 2,
                                     Q = list(matrix(NA), matrix(0))) +
                            SSMcustom(Z = 1, T = 0, R = 1, Q = NA),
                        H = known_sampling_var)

  initial_guesses <- log(c(var(diff(y)), var(y) / 10))

  fit <- fitSSM(model_spec, inits = initial_guesses, hessian = TRUE, ...)
  out <- KFS(fit$model)

  TT         <- length(y)
  d.hat      <- out$alphahat[TT, "slope"]
  d.var      <- out$V[2, 2, TT]
  d.se       <- sqrt(d.var)
  var_innov  <- fit$model$Q[1, 1, 1]
  var_shock  <- fit$model$Q[3, 3, 1]
  var_samp   <- fit$model$H[1, 1, 1]

  z        <- qnorm(0.975)
  hess     <- fit$optim.out$hessian
  vcov_log <- solve(hess)
  se_log       <- sqrt(diag(vcov_log))
  theta        <- fit$optim.out$par
  ci_var_innov <- exp(theta[1] + c(-z, z) * se_log[1])
  ci_var_shock <- exp(theta[2] + c(-z, z) * se_log[2])

  ci_drift <- d.hat + c(-z, z) * d.se

  names(d.hat) <- names(var_innov) <- names(var_shock) <- names(var_samp) <- ""

  return(c(d            = d.hat,
           var_innov    = var_innov,
           var_samp     = var_samp,
           var_shock    = var_shock,
           d_lo         = ci_drift[1],
           d_hi         = ci_drift[2],
           var_innov_lo = ci_var_innov[1],
           var_innov_hi = ci_var_innov[2],
           var_shock_lo = ci_var_shock[1],
           var_shock_hi = ci_var_shock[2]))
}

marss_fit_for_pooled <-
    function(Y, sigma, sd_ratio = NULL, max_iter = 500, inits = NULL)
{
    ## Y ~ data, 1 or many us states
    ## sigma ~ vector of sampling_sds, 1 or more
    ## sd_ratio of NULL means no constraint. I used 3 for pooled fitting.
    ## not 100% sure of "inits"
    
    # 1. Robust Data Shaping
    # Force to matrix; if it's a vector, make it 1 row (n=1)
    if (is.null(dim(Y))) {
        Y <- matrix(Y, nrow = 1)
    } else {
        Y <- as.matrix(Y)
    }
    n <- nrow(Y)

    # 2. Setup Variance Parameters (Q)
    # Use list-matrix to allow mixing numbers and character names
    Q_mat <- matrix(list(0), 2 * n, 2 * n)
    # Define names/values based on sd_ratio
    if (is.null(sd_ratio)) {
        ## free 
        q_persist_val <- "q_p"
        q_transit_val <- "q_t"
    }
    if (!is.null(sd_ratio)) {
        ## ratio constrained
        var_ratio = sd_ratio^2
        q_persist_val <- paste0(var_ratio, "*q")
        q_transit_val <- "q"
    }
    # Fill the diagonal of the Q list-matrix
    for(i in 1:n) {
        Q_mat[i, i] <- q_persist_val           # Top half: Persistent
        Q_mat[i + n, i + n] <- q_transit_val   # Bottom half: Transient
    }
    
    # 3. Setup Drift Parameters (U)
    # Also a list-matrix to hold character names
    U_mat <- matrix(list(0), 2 * n, 1)
    for(i in 1:n) {
        U_mat[i, 1] <- paste0("u", i)
    }
    ## and stays 0 for i > n?

    # 4. Final Model List
    # Z and B are fixed numbers, so they don't strictly need to be lists, 
    # but R is numeric here, so it's safe.
    model_list <- list(
        Z = cbind(diag(n), diag(n)),       # [I | I]
        B = diag(c(rep(1, n), rep(0, n))), #  persistent (1), transient (0)
        Q = Q_mat,
        U = U_mat,
        R = diag(as.numeric(sigma)^2, n), # Fixed known sampling variance
        A = "zero",
        x0 = matrix(c(Y[,1], rep(0, n)), 2 * n, 1) # Initial states
    )
    
    # 5. Fit
    fit <- MARSS(Y, model = model_list, method = "BFGS", 
                 inits = inits, control = list(maxit = max_iter))
    return(fit)
}
## This is code for Paige.

## It's very rough. The idea is to figure out which specifications of
## StructTS correspond to which models.

## The code below provides some wrapper functions that can be used instead
## of having to figure out what spec to use with StructTS

## There's also examples of simulation to make sure that our wrapper functions are working as we want.

## (1) RWD
N = 400
y <- NULL
y[1] = 0
d = 1; d.hat = NULL;
var_innov <- 2; var_innov.hat = NULL;
eps_innov <- rnorm(N, mean = 0, sd = sqrt(var_innov))
for (i in 1:(N-1))
{
    y[i+1] = y[i] + d + eps_innov[i]
}
plot(y)
sts = StructTS(x = y, type = "trend", fixed = c(NA, 0, 0))
StructTS(x = y, type = "level")
var_innov.hat = sts$coef["level"]
d.hat = sts$model$a[2]
result = cbind("true" = c(var_innov = var_innov,
                 d = d),
               "hat" = round(c(var_innov.hat,
                 d.hat),2))
print(result)

rwd <- function(y)
{
    sts = StructTS(x = y, type = "trend", fixed = c(NA, 0, 0))
    d.hat = sts$model$a[2]
    var_innov.hat = sts$coef["level"]
    return(list(d = d.hat,
                var_innov = var_innov.hat))
}

## (2) RWD with observation error
rwd_with_obs_error <- function(y, ..., maxit = 3000)
{
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

simu_rwd_with_obs_error <- function(N = 100, x0 = 0,
                                    d, var_innov,
                                    var_obs)
{
    eps_innov <- rnorm(N, mean = 0, sd = sqrt(var_innov))
    nt <- eps_obs <- rnorm(N, mean = 0, sd = sqrt(var_obs))
    x <- NULL
    y <- NULL;
    x[1] = x0
    y[1] = x[1] + nt[1]
    for (i in 1:(N-1))
    {
        x[i+1] = x[i] + d + eps_innov[i]
        y[i] = x[i] + nt[i]
    }
    return(y)
}

rwd_with_obs_error_marss <- function(y, ...)
{
    library(MARSS)
    mod.list.2 <- list(B = matrix(1),
                       U = matrix("d"),
                       Q = matrix("q"),
                       Z = matrix(1),
                       A = matrix(0),
                       R = matrix("r"),
                       x0 = matrix("mu"),
                       tinitx = 0)
    out2.marss = MARSS(y, model = mod.list.2, ...)
    d = coef(out2.marss)$U[1,1]
    var_innov = coef(out2.marss)$Q[1,1]
    var_obs = coef(out2.marss)$R[1,1]
    names(d) <- names(var_innov) <- names(var_obs) <- ""
        return(c(d= d,
                 var_innov= var_innov,
                 var_obs = var_obs))
}

rwd_with_obs_error_marss_bfgs <- function(y, ...)
{
  library(MARSS)
  mod.list.2 <- list(B = matrix(1),
                     U = matrix("d"),
                     Q = matrix("q"),
                     Z = matrix(1),
                     A = matrix(0),
                     R = matrix("r"),
                     x0 = matrix("mu"),
                     tinitx = 0)
  out2.marss = MARSS(y, model = mod.list.2, method = "BFGS", ...)
  d = coef(out2.marss)$U[1,1]
  var_innov = coef(out2.marss)$Q[1,1]
  var_obs = coef(out2.marss)$R[1,1]
  names(d) <- names(var_innov) <- names(var_obs) <- ""
  return(c(d= d,
           var_innov= var_innov,
           var_obs = var_obs))
}

d = .1 # the rounded average estimated drift from the state models (original spec was  .4)
var_innov = .1 # the rounded average estimated var_innov from the state models (original spec was 2)
var_obs = 0.001

## standard example
y <- simu_rwd_with_obs_error(N = 1400,
                             d = d, var_innov = var_innov, var_obs = var_obs)
plot(y)

( out <- rwd_with_obs_error(y) )

##         d var_innov   var_obs 
## 0.3794464 2.2167909 1.1101556 


    
out.sts = rwd_with_obs_error(y) 
out.marss = rwd_with_obs_error_marss(y)
out.bfgs = rwd_with_obs_error_marss_bfgs(y)
rbind(out.sts, out.marss, out.bfgs)
##                   d var_innov   var_obs
## out.sts   0.4059950  2.205932 0.9527210
## out.marss 0.4076975  2.198683 0.9361897


## MARSS fit is
## Estimation method: kem 
## Convergence test: conv.test.slope.tol = 0.5, abstol = 0.001
## Estimation converged in 34 iterations. 
## Log-likelihood: -2974.028 
## AIC: 5956.056   AICc: 5956.084   
 
##       Estimate
## R.r      1.108
## U.d      0.381
## Q.q      2.196
## x0.mu   -0.682
## Initial states (x0) defined at t=0

## Standard errors have not been calculated. 
## Use MARSSparamCIs to compute CIs and bias estimates.

### CREATE ILLUSTRATIVE FIGURE

# run multiple simulations 
var_obs_values = seq(0.00001, .5, length.out = 200)
results = matrix(NA, nrow = length(var_obs_values), ncol = 4)
colnames(results) = c("True", "StructTS", "MARSS_kem", "MARSS_BFGS")

for (i in seq_along(var_obs_values)) {
  var_obs <- var_obs_values[i]
  
  y <- simu_rwd_with_obs_error(N = 1400,
                               d = d,var_innov = var_innov, var_obs = var_obs)
  out.sts = rwd_with_obs_error(y) 
  out.marss = rwd_with_obs_error_marss(y)
  out.bfgs = rwd_with_obs_error_marss_bfgs(y)
  
  results[i, "True"] = var_obs
  results[i, "StructTS"] = out.sts["var_obs"]
  results[i, "MARSS_kem"] <- out.marss["var_obs"]
  results[i, "MARSS_BFGS"] <- out.bfgs["var_obs"]
}

# Convert to a data frame
results <- as.data.frame(results)

# figure showing different true observation errors and 3 methods estimates 
library(ggplot2)
ggplot(results, aes(x = True)) +
  geom_line(aes(y=True), color = "black") +
  geom_line(aes(y=StructTS), color = "blue") +
  geom_line(aes(y=MARSS_kem), color = "red") +
  geom_line(aes(y=MARSS_BFGS), color = "green") +
  labs(x = "True Observation Variance", y = "Estimated Observation Variance") +
  scale_color_manual(values = c("black", "blue", "red", "green"), 
                     labels = c("True", "StructTS", "MARSS_kem", "MARSS_BFGS")) +
  theme_minimal()

## rwd only example
y <- simu_rwd_with_obs_error(N = 1000,
                             d = d,var_innov = var_innov, var_obs = 0)
rwd_with_obs_error(y, maxit = 5000)

## test things with other time series data 
dax_matrix <- matrix(EuStockMarkets[, "DAX"], ncol=1)
smi_matrix <- matrix(EuStockMarkets[, "SMI"], ncol=1)
# var_obs is 0 with StructTS 
rwd_with_obs_error(dax_matrix)
rwd_with_obs_error(smi_matrix)
# var_obs is very small, but non-0 with MARSS bfgs
rwd_with_obs_error_marss_bfgs(t(dax_matrix))
rwd_with_obs_error_marss_bfgs(t(smi_matrix))
# neither converged with MARSS kem
rwd_with_obs_error_marss(t(dax_matrix))
rwd_with_obs_error_marss(t(smi_matrix))

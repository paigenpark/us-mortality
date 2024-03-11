## This code creates several wrapper functions including the specifications of StructTS() and MARSS()
## that correspond to models I use in a US mortality paper 

## Since I've been having issues with observation variance estimation in the US data analysis, 
## I also write a simulation to test whether each method is accurately estimating obs variance

## (1) RWD wrapper function
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

## (2) RWD with observation error wrapper functions 

## StructTS
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

## MARSS default (kem)
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

## MARSS BFGS
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


## (3) function to set up simulated y, which follows local linear trend model
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

## setting arbitrary values for drift, var_innov, and var_obs terms 
d = .1 
var_innov = .1
var_obs = .1

## standard example
y <- simu_rwd_with_obs_error(N = 1400,
                             d = d, var_innov = var_innov, var_obs = var_obs)
plot(y)
   
out.sts = rwd_with_obs_error(y) 
out.marss = rwd_with_obs_error_marss(y)
out.bfgs = rwd_with_obs_error_marss_bfgs(y)
rbind(out.sts, out.marss, out.bfgs)

## (4) create illustrative figure

## run multiple simulations with varying var_obs values
var_obs_values = seq(0.01, 0.1, length.out = 100)
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

## convert to a data frame
results <- as.data.frame(results)

## figure plotting estimated vs. true observation errors for 3 methods  
library(ggplot2)
ggplot(results, aes(x = True)) +
  geom_line(aes(y=True, color="True")) +
  geom_line(aes(y=StructTS, color="StructTS")) +
  geom_line(aes(y=MARSS_kem, color="MARSS_kem")) +
  geom_line(aes(y=MARSS_BFGS, color="MARSS_BFGS")) +
  labs(x = "True Observation Variance", y = "Estimated Observation Variance") +
  scale_color_manual(values = c("True" = "black", "StructTS" = "blue", 
                                "MARSS_kem" = "red", "MARSS_BFGS" = "green"),
                     labels = c("MARSS BFGS", "MARSS KEM", "StructTS", "True Observation Variance")) +
  theme_minimal()



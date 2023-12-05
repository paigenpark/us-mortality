# STEP 1: import packages and data
library(matrixStats)
library(pracma)

  # read in life tables for US
path = "/Users/paigepark/Desktop/us-mortality/data"
us_data = read.csv(paste(path, "USA_bltper_1x1.csv", sep = "/"))
us_e0 = us_data$ex[us_data$Age == 0]
#no2020 = us_e0[1:61]

# deal with missing data
  # this step isn't necessary for this project since there won't be any missing data

# STEP 2: set up stationarity check and differencing function
get_differenced_data <- function(data, ...) {
  library(tseries)
  
  # tests the null hypothesis of data being a non-stationary time series
  suppressWarnings(station_test <- adf.test(data, alternative = "stationary"))
  
  while (station_test$p.value > 0.05) {
    print("Data is not stationary, differencing data")
    data <- diff(data)
    suppressWarnings(station_test <- adf.test(data, alternative = "stationary"))
  }
  
  return(data)
}

# STEP 3: set up matrices and state vector to be used in models
makeTrend <- function(x, d=0) {
  n = length(x)
  a = c(x[1], 0) 
  Z = c(1,0)
  T = matrix(c(1, 0, 1, 1), 2L, 2L)
  h = 1
  P <- Pn <- matrix(0, 2L, 2L)
  V <- diag(2L)
  return(list(Z = Z, a = a, P = P, T = T, V = V, h = h, d = d,
              Pn = Pn))
}

# STEP 4: Implement Custom Kalman Filter
# Kalman filter is an iterative process of prediction and correction
# It gives the best estimate of the current state given past observations

kalman_filter <- function(y, model) {
  # get components from model
  Z <- as.matrix(model$Z)
  a <- as.matrix(model$a)
  T <- model$T
  h <- model$h
  P <- model$P
  V <- model$V
  d <- model$d
  
  # initialize
  n = length(y)
  a_updated = matrix(0, 2L, n)
  V_updated = array(0, dim = c(2, 2, n))

  for (t in 1:n) {  
    # prediction step
    if (t == 1) {
    a_pred = a[,t] # at time 1, we use the value for a initialized in the makeTrend() function
    V_pred = V # similarly, we use initialized value for state covariance 
      } else {
    a_pred = t(T %*% a_updated[,t-1] + d) # multiply prior state by transition matrix and add drift 
    V_pred = (T %*% V_updated[,,t-1] %*% t(T)) + (V * h) # multiply prior state covariance by transition
      # matrix squared and add the process covariance
    }
    
    # update step
    innov = y[t] - Z %*% a_pred
    innov_cov = (t(Z) %*% V_pred %*% Z) + h
    k_gain = V_pred %*% Z %*% solve(innov_cov) # determines how much prediction should be updated based on 
      # observed data
    
    a_updated[,t] = a_pred + t(k_gain) %*% innov
    V_updated[,,t] = V_pred - k_gain %*% t(Z) %*% V_pred
  }
  return(list(a = a_updated, V = V_updated))
}

# STEP 5: Implement Custom Kalman Smoother
# improves state estimates by incorporating future observations

kalman_smoother <- function(T, a, P, a_priori, P_priori) {
  # Number of observations
  n <- nrow(a_priori)
  
  # Initialize smoothed estimates with the final filtered estimates
  a_s <- a
  P_s <- P
  
  # Loop over each time step in reverse order from n-1 to 1
  for (t in (n-1):1) {
    # Backward pass
    L <- P[t,,] %*% t(T) %*% solve(P_priori[t+1,,])
    a_s[t,] <- a[t,] + L %*% (a_s[t+1,] - a_priori[t+1,])
    P_s[t,,] <- P[t,,] + L %*% (P_s[t+1,,] - P_priori[t+1,,]) %*% t(L)
  }
  
  return(list(a = a_s, P = P_s))
}

# STEP 6: Parameter estimation (using MLE like StructTS does)

# first we want to define the log-likelihood of the data using the Kalman filter
logLik <- function(params, y, model) {
  # unpack the parameters
  P <- params[1]
  V <- params[2]
  d <- params[3]
  
  # update model with new parameters
  model$P <- P
  model$V <- V
  model$d <- d
  
  # run the Kalman filter
  results <- kalman_filter(y, model)
  
  # calculate the likelihood
  log_like = -0.5 * sum(log(results$V) + (y - results$a)^2 / results$V)
  
  return(-log_like)
}

# running functions on e0 data
differenced_data = get_differenced_data(us_e0)
initial_matrices = makeTrend(differenced_data)

# give initial guess for P, V, and d
initial_guess <- c(1, 1, 0)

# Optimize the negative log-likelihood function
opt <- optim( par = initial_guess, 
              fn = logLik, 
              y = us_e0, 
              model = initial_matrices)

# the MLEs are stored in mle_results$par
mle_P <- mle_results$par[1]
mle_V <- mle_results$par[2]
mle_d <- mle_results$par[3]

# Print the estimated parameters
#print(estimated_parameters)

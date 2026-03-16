### we reproduce table 1 for a subset of states so it will run fast
### unpooled AND pooled
### and compare BIC, AIC etc

### extra credit: we show how it works differently for different ratios
### extra, extra credit: we show unpooled fixed ratio estimators give
## about the same as pooled ratio estimators.

library(data.table)
library(MARSS)

### 0. Read data and select sample

Y_full <- dget("data/state_data.txt")
sigma_full <- dget("data/sampling_error.txt") ## from Paige's old Table 1
## select sample
my_states = rownames(Y_full)[12:19]
## [1] "IN" "IA" "KS" "KY" "LA" "ME" "MD" "MA"
Y <- Y_full[my_states,paste0(1941:2019)]
sigma.vec <- sigma_full[my_states]

### 1. general function to run marss

marss_fit_fun <-
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

### 2. fitting

## objects to save results in
n_regions = nrow(Y)
unpooled_list <- vector("list", n_regions)
names(unpooled_list) <- my_states
ratio_unpooled_list <- vector("list", n_regions)
names(ratio_unpooled_list) <- my_states


### a. unpooled
for (i in 1:nrow(Y))
    unpooled_list[[i]] <- marss_fit_fun(Y[i,],
                                        sigma = sigma.vec[i],
                                        sd_ratio = NULL)

### b. unpooled, constrainted to ratio 3.0
for (i in 1:nrow(Y))
    ratio_unpooled_list[[i]] <- marss_fit_fun(Y[i,],
                                        sigma = sigma.vec[i],
                                        sd_ratio = 3.0)

### c. pooled, ratio 3.0
pooled_3.0 <- marss_fit_fun(Y,
                            sigma = sigma.vec,
                            sd_ratio = 3.0)

### 4. analyze fits

### a. look at unpooled coefs
unpooled_Q <- sapply(unpooled_list,
                            function(x) coef(x)$Q)
dimnames(unpooled_Q) <-
    list(c("q_p", "q_t"), my_states)

print(round(unpooled_Q, 3))
##        IN    IA    KS    KY    LA    ME    MD    MA
## q_p 0.066 0.044 0.046 0.123 0.158 0.094 0.110 0.085
## q_t 0.009 0.006 0.000 0.000 0.001 0.000 0.013 0.027


ratio_unpooled_Q.vec <- sapply(ratio_unpooled_list,
                               function(x) coef(x)$Q)
ratio_unpooled_Q =
    rbind(9*ratio_unpooled_Q.vec, ## persistant
          ratio_unpooled_Q.vec) ## transitory
dimnames(ratio_unpooled_Q) <-
    list(c("q_p", "q_t"), my_states)
print(round(ratio_unpooled_Q, 3))

out_Q_unpooled = t(rbind(unpooled_Q, ratio_unpooled_Q))

unpooled_sd.mat = round(sqrt(out_Q_unpooled), 3)
print(unpooled_sd.mat)

##      q_p   q_t   q_p   q_t
## IN 0.257 0.094 0.261 0.087
## IA 0.210 0.078 0.214 0.071
## KS 0.214 0.000 0.202 0.067
## KY 0.351 0.000 0.324 0.108
## LA 0.398 0.035 0.368 0.123
## ME 0.307 0.000 0.290 0.097
## MD 0.332 0.115 0.334 0.111
## MA 0.291 0.165 0.333 0.111

## looks very good, like we get much better transient q's
## without changing the persistent ones much

## plot results
plot(unpooled_Q[1,],
     unpooled_Q[2,],
     xlab = "process var",
     ylab = "shock var")
text(unpooled_Q[1,],
     unpooled_Q[2,],
     my_states, pos = 4)
points(ratio_unpooled_Q[1,],
       ratio_unpooled_Q[2,],
       col = 'red')
text(ratio_unpooled_Q[1,],
     ratio_unpooled_Q[2,],
     my_states, pos = 1,
     col = 'red')
title("regularization!")
## now the pooled
q_t_pooled = pooled_3.0$par$Q["q",]
q_p_pooled = q_t_pooled * 3^2
points(q_p_pooled, q_t_pooled, col = "darkgreen",
       pch = 19, cex = 2)
text(q_p_pooled, q_t_pooled, "pooled",
       col = "darkgreen", pos = 3)
## could redo for SD

## but this figure already seems very reasonable
## with red (ratio unpooled) regularizign the unpooled
## and the green (pooled) falling right in the middle

### compare AIC, AICc, BIC

## -------- Note: number of parameters needs to be increase
## -------- to account for the "ratio", so
## -------- AIC = 2 k - 2 ln(L)
## -------- means we AIC_adj = AIC + 2 (since k increases by 1)

## easy to do for individual unpooled

ratio_unpooled_AIC.vec <- 2 + sapply(ratio_unpooled_list,
                               AIC) ## we adjust
unpooled_AIC.vec <- sapply(unpooled_list, AIC)
out = cbind(unpooled_AIC.vec,
            ratio_unpooled_AIC.vec)
rownames(out) = my_states
print(round(out,3))

##      unpooled_AIC.vec ratio_unpooled_AIC.vec
## [1,] 39.7631          37.7899               
## [2,] 18.6993          16.7313               
## [3,] 11.0098          9.70219               
## [4,] 71.8305          70.7019               
## [5,] 90.9558          89.5049               
## [6,] 64.9467          64.3503               
## [7,] 77.0873          75.0914               
## [8,] 74.596           73.7092               

## total
colSums(out)
## unpooled_AIC.vec ratio_unpooled_AIC.vec 
##          448.888                453.581 

AIC(pooled_3.0) + 2
## [1] 464.562


## the unpooled estimator wins in terms of AIC
## probably holds for more US states?.
## But probably overfitting and so forecast may
## do worse than pooled model?


## gemini code for grid search

# 1. Define the grid of ratios to test
ratio_grid <- seq(1.5, 4.0, by = 0.1)

# 2. Run the grid search
grid_results <- lapply(ratio_grid, function(r) {
    fit <- marss_fit_fun(Y, sigma = sigma.vec, sd_ratio = r, max_iter = 300)
    
    # Return a clean data frame row
    # Note: Adding +2 to AIC to account for the ratio being a chosen parameter
    data.frame(ratio = r, 
               logLik = fit$logLik, 
               AIC_adj = fit$AIC + 2, 
               q_val = fit$par$Q["q",])
})

# 3. Combine and find the "Best"
grid_df <- do.call(rbind, grid_results)
best_fit <- grid_df[which.min(grid_df$AIC_adj), ]

print(grid_df)
cat("\nOptimal Ratio based on AIC:", best_fit$ratio, "\n")

# 4. Simple plot to see the "valley"
plot(grid_df$ratio, grid_df$AIC_adj, type="b", 
     xlab="SD Ratio", ylab="Adjusted AIC", main="Profile Likelihood of Ratio")


# 4. zoom of Simple plot to see the "valley"
plot(grid_df$ratio, grid_df$AIC_adj, type="b", 
     xlab="SD Ratio", ylab="Adjusted AIC", main="Profile Likelihood of Ratio", 
xlim = c(3,4), ylim = c(464, 466))

## very, very flat AIC_adj suggestion is 3.5 is optimal




### random walk with drift (single-state, no persistent/transitory split)

rw_drift_fit_fun <- function(Y, sigma, max_iter = 500) {
  if (is.null(dim(Y))) {
    Y <- matrix(Y, nrow = 1)
  } else {
    Y <- as.matrix(Y)
  }
  n <- nrow(Y)
  model_list <- list(
    Z = diag(n),
    B = diag(n),
    Q = "diagonal and equal",
    U = matrix(paste0("u", 1:n), n, 1),
    R = diag(as.numeric(sigma)^2, n),
    A = "zero",
    x0 = matrix(Y[, 1], n, 1)
  )
  MARSS(Y, model = model_list, method = "BFGS",
        control = list(maxit = max_iter))
}

######## RANDOM HOLD-OUT IMPUTATION EXPERIMENT ##########

# Set a seed so you get the same random holes every time you run it!
set.seed(42) 

## 1. Create the "Holey" Dataset with Random Drops
Y_matrix <- Y_full[my_states, paste0(1941:2019)]
n_states <- nrow(Y_matrix)
n_years  <- ncol(Y_matrix)
n_total  <- n_states * n_years

drop_fraction <- 0.20 # We will drop 20% of the observations
n_drop <- floor(drop_fraction * n_total)

# Create a logical mask (TRUE = keep, FALSE = drop)
mask <- matrix(TRUE, nrow = n_states, ncol = n_years, dimnames = dimnames(Y_matrix))
drop_indices <- sample(1:n_total, size = n_drop)
mask[drop_indices] <- FALSE

# Apply the mask to create our test dataset
Y_holey <- Y_matrix
Y_holey[!mask] <- NA

## 2. Helper function to apply fixed parameters to holey data
## 2. Helper function to apply fixed parameters to holey data
## 2. Helper function to apply fixed parameters to holey data
nowcast_eval <- function(fit_full, Y_holey_data) {
    fit_impute <- fit_full
    
    # Safely inject the new data (with NAs)
    if (!is.null(fit_impute$model$data)) {
        fit_impute$model$data <- as.matrix(Y_holey_data)
    }
    if (!is.null(fit_impute$marss$data)) {
        fit_impute$marss$data <- as.matrix(Y_holey_data)
    }
    
    # 1. Run the core Kalman filter/smoother to get the hidden STATES (xtT)
    kf_out <- MARSSkf(fit_impute)
    
    # 2. Extract the Z and A matrices from the fitted model
    Z_mat <- coef(fit_impute, type = "matrix")$Z
    A_mat <- coef(fit_impute, type = "matrix")$A
    
    # 3. Calculate expected observations: Y = Z * X + A
    # Z_mat %*% kf_out$xtT gives an N x T matrix. 
    # Adding as.vector(A_mat) adds the state-specific offsets properly to each year.
    fit_impute$ytT <- Z_mat %*% kf_out$xtT + as.vector(A_mat)
    
    return(fit_impute)
}

## 3. Run the imputation using the models ALREADY fit to the full data
# (Assumes unpooled_list, ratio_unpooled_list, pooled_3.0, and rw_full_list exist)

unpooled_nowcast_list <- vector("list", n_states)
ratio_unpooled_nowcast_list <- vector("list", n_states)
rw_nowcast_list <- vector("list", n_states)

for (i in 1:n_states) {
    unpooled_nowcast_list[[i]] <- nowcast_eval(unpooled_list[[i]], Y_holey[i, , drop = FALSE])
    ratio_unpooled_nowcast_list[[i]] <- nowcast_eval(ratio_unpooled_list[[i]], Y_holey[i, , drop = FALSE])
    rw_nowcast_list[[i]] <- nowcast_eval(rw_full_list[[i]], Y_holey[i, , drop = FALSE])
}

pooled_nowcast <- nowcast_eval(pooled_3.0, Y_holey)

## 4. Extract Estimates 
unpooled_nc <- matrix(NA, n_states, n_years, dimnames = dimnames(Y_matrix))
ratio_unpooled_nc <- matrix(NA, n_states, n_years, dimnames = dimnames(Y_matrix))
rw_nc <- matrix(NA, n_states, n_years, dimnames = dimnames(Y_matrix))

# Pooled outputs an N x T matrix natively
pooled_nc <- pooled_nowcast$ytT
rownames(pooled_nc) <- my_states

for (i in 1:n_states) {
    unpooled_nc[i, ]       <- unpooled_nowcast_list[[i]]$ytT[1, ]
    ratio_unpooled_nc[i, ] <- ratio_unpooled_nowcast_list[[i]]$ytT[1, ]
    rw_nc[i, ]             <- rw_nowcast_list[[i]]$ytT[1, ]
}

## 5. Evaluate: RMSE specifically on the held-out values

# We update the RMSE function to handle NAs so we can isolate the dropped points
rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2, na.rm = TRUE))

# We create matrices that ONLY contain values where data was dropped.
# Everything else becomes NA, so the mean() function naturally ignores the data the model was allowed to see.
actual_held_out <- Y_matrix
actual_held_out[mask] <- NA 

pred_unpooled_held_out <- unpooled_nc
pred_unpooled_held_out[mask] <- NA

pred_ratio_held_out <- ratio_unpooled_nc
pred_ratio_held_out[mask] <- NA

pred_pooled_held_out <- pooled_nc
pred_pooled_held_out[mask] <- NA

pred_rw_held_out <- rw_nc
pred_rw_held_out[mask] <- NA

# Calculate state-by-state RMSE
rmse_unpooled <- sapply(1:n_states, function(i) rmse(actual_held_out[i, ], pred_unpooled_held_out[i, ]))
rmse_ratio_unpooled <- sapply(1:n_states, function(i) rmse(actual_held_out[i, ], pred_ratio_held_out[i, ]))
rmse_pooled <- sapply(1:n_states, function(i) rmse(actual_held_out[i, ], pred_pooled_held_out[i, ]))
rmse_rw <- sapply(1:n_states, function(i) rmse(actual_held_out[i, ], pred_rw_held_out[i, ]))

rmse_table <- data.frame(
    state = my_states,
    rw_drift = round(rmse_rw, 3),
    unpooled = round(rmse_unpooled, 3),
    ratio_unpooled_3 = round(rmse_ratio_unpooled, 3),
    pooled_3 = round(rmse_pooled, 3)
)
print(rmse_table)

cat("\nOverall Random Imputation RMSE (20% missing):\n")
cat("  RW with drift:     ", round(rmse(actual_held_out, pred_rw_held_out), 3), "\n")
cat("  Unpooled:          ", round(rmse(actual_held_out, pred_unpooled_held_out), 3), "\n")
cat("  Ratio unpooled 3:  ", round(rmse(actual_held_out, pred_ratio_held_out), 3), "\n")
cat("  Pooled 3:          ", round(rmse(actual_held_out, pred_pooled_held_out), 3), "\n")

######## RANDOM HOLD-OUT VISUALIZATION ##########

# Make sure the "results" directory exists
if (!dir.exists("results")) dir.create("results")

years_all <- as.numeric(colnames(Y_matrix))

for (plot_state in my_states) {
    
    # 1. Extract the full actual time series for the background
    actual_full <- as.numeric(Y_matrix[plot_state, ])
    
    # Identify which years were dropped for this specific state
    dropped_idx <- which(!mask[plot_state, ])
    dropped_years <- years_all[dropped_idx]
    
    # Extract the actual values that were hidden (Ground Truth)
    actual_dropped <- actual_full[dropped_idx]
    
    # Extract predictions just for the dropped years
    pred_unpooled <- as.numeric(pred_unpooled_held_out[plot_state, dropped_idx])
    pred_ratio    <- as.numeric(pred_ratio_held_out[plot_state, dropped_idx])
    pred_pooled   <- as.numeric(pred_pooled_held_out[plot_state, dropped_idx])
    pred_rw       <- as.numeric(pred_rw_held_out[plot_state, dropped_idx])

    # 2. Setup Plot
    pdf(file.path("results", paste0("random_nowcast_", plot_state, ".pdf")), 
        width = 9, height = 5)
    par(mar = c(4, 4, 3, 1))
    
    # Find the min/max of all values so the y-axis scales perfectly
    ylim <- range(c(actual_full, pred_unpooled, pred_ratio, pred_pooled, pred_rw), na.rm = TRUE)
    
    # Plot the continuous actual line in a subtle gray
    plot(years_all, actual_full, type = "l", col = "gray60", lwd = 1.5,
         ylim = ylim, xlab = "Year", ylab = "Life expectancy (e0)",
         main = paste0(plot_state, ": Random Hold-Out Imputation (20% missing)"))
    
    # Add solid small dots for the actual data that the model *was* allowed to see
    kept_idx <- which(mask[plot_state, ])
    points(years_all[kept_idx], actual_full[kept_idx], pch = 16, col = "black", cex = 0.6)
    
    # Add larger open circles for the actual data that was hidden (Ground Truth)
    points(dropped_years, actual_dropped, pch = 1, col = "black", cex = 1.3, lwd = 1.5)

    # 3. Overlay model predictions at the missing spots using colored X's
    # We offset the years slightly on the x-axis so the X's don't perfectly overlap
    points(dropped_years - 0.2, pred_rw,       pch = 4, col = "blue",      cex = 0.9, lwd = 1.5)
    points(dropped_years,       pred_unpooled, pch = 4, col = "red",       cex = 0.9, lwd = 1.5)
    points(dropped_years + 0.2, pred_ratio,    pch = 4, col = "orange",    cex = 0.9, lwd = 1.5)
    points(dropped_years + 0.4, pred_pooled,   pch = 4, col = "darkgreen", cex = 0.9, lwd = 1.5)

    # 4. Legend
    # Pull the exact RMSE values from the table we built in the previous step
    rmse_state <- rmse_table[rmse_table$state == plot_state, ]
    
    legend("topleft", bty = "n", cex = 0.8,
           legend = c("Observed Data", 
                      "Hidden Data (Ground Truth)",
                      paste0("RW drift (RMSE=", rmse_state$rw_drift, ")"),
                      paste0("Unpooled (RMSE=", rmse_state$unpooled, ")"),
                      paste0("Ratio unpooled (RMSE=", rmse_state$ratio_unpooled_3, ")"),
                      paste0("Pooled (RMSE=", rmse_state$pooled_3, ")")),
           col = c("black", "black", "blue", "red", "orange", "darkgreen"),
           pch = c(16, 1, 4, 4, 4, 4),
           pt.lwd = c(1, 1.5, 1.5, 1.5, 1.5, 1.5))
    
    dev.off()
    cat("Saved:", file.path("results", paste0("random_nowcast_", plot_state, ".pdf")), "\n")
}
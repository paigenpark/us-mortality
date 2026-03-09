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

######## forecasting experiment ##########

## split data into training and test sets
train_years <- paste0(1941:2004)
test_years  <- paste0(2005:2019)

Y_train <- Y_full[my_states, train_years]
Y_test  <- Y_full[my_states, test_years]
print(Y_train[, 1:5])
print(Y_test[, 1:5])

## fit three models on training data

## a. unpooled (free ratio)
unpooled_train_list <- vector("list", n_regions)
names(unpooled_train_list) <- my_states
for (i in 1:n_regions)
    unpooled_train_list[[i]] <- marss_fit_fun(Y_train[i, ],
                                              sigma = sigma.vec[i],
                                              sd_ratio = NULL)

## b. unpooled, constrained ratio 3
ratio_unpooled_train_list <- vector("list", n_regions)
names(ratio_unpooled_train_list) <- my_states
for (i in 1:n_regions)
    ratio_unpooled_train_list[[i]] <- marss_fit_fun(Y_train[i, ],
                                                    sigma = sigma.vec[i],
                                                    sd_ratio = 3.0)

## c. pooled, ratio 3
pooled_train <- marss_fit_fun(Y_train,
                              sigma = sigma.vec,
                              sd_ratio = 3.0)

## d. unpooled random walk with drift
rw_train_list <- vector("list", n_regions)
names(rw_train_list) <- my_states
for (i in 1:n_regions)
  rw_train_list[[i]] <- rw_drift_fit_fun(Y_train[i, ],
                                          sigma = sigma.vec[i])

## forecast: use MARSS predict for h-step ahead
n_test <- length(test_years)

## helper: extract forecast from a single-state MARSS fit
forecast_single <- function(fit, h) {
    pred <- predict(fit, type = "ytT", n.ahead = h)
    ## pred$pred is a data.frame with columns .rownames, t, y, estimate, ...
    ## take only the last h rows (the forecast, not fitted values)
    est <- pred$pred$estimate
    tail(est, h)
}

## a. unpooled forecasts
unpooled_fc <- matrix(NA, n_regions, n_test,
                      dimnames = list(my_states, test_years))
for (i in 1:n_regions)
    unpooled_fc[i, ] <- forecast_single(unpooled_train_list[[i]], n_test)

## b. ratio-unpooled forecasts
ratio_unpooled_fc <- matrix(NA, n_regions, n_test,
                            dimnames = list(my_states, test_years))
for (i in 1:n_regions)
    ratio_unpooled_fc[i, ] <- forecast_single(ratio_unpooled_train_list[[i]],
                                              n_test)

## c. pooled forecasts
pooled_pred <- predict(pooled_train, type = "ytT", n.ahead = n_test)
pooled_fc <- matrix(NA, n_regions, n_test,
                    dimnames = list(my_states, test_years))
for (i in 1:n_regions) {
    sub <- pooled_pred$pred[pooled_pred$pred$.rownames == my_states[i], ]
    pooled_fc[i, ] <- tail(sub$estimate, n_test)
}

## d. random walk with drift forecasts
rw_fc <- matrix(NA, n_regions, n_test,
                dimnames = list(my_states, test_years))
for (i in 1:n_regions)
    rw_fc[i, ] <- forecast_single(rw_train_list[[i]], n_test)

## evaluate: RMSE by state and overall
rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2))

rmse_unpooled <- sapply(1:n_regions, function(i)
    rmse(as.numeric(Y_test[i, ]), unpooled_fc[i, ]))
rmse_ratio_unpooled <- sapply(1:n_regions, function(i)
    rmse(as.numeric(Y_test[i, ]), ratio_unpooled_fc[i, ]))
rmse_pooled <- sapply(1:n_regions, function(i)
    rmse(as.numeric(Y_test[i, ]), pooled_fc[i, ]))
rmse_rw <- sapply(1:n_regions, function(i)
    rmse(as.numeric(Y_test[i, ]), rw_fc[i, ]))

rmse_table <- data.frame(
    state = my_states,
    rw_drift = round(rmse_rw, 3),
    unpooled = round(rmse_unpooled, 3),
    ratio_unpooled_3 = round(rmse_ratio_unpooled, 3),
    pooled_3 = round(rmse_pooled, 3)
)
print(rmse_table)
cat("\nOverall RMSE:\n")
cat("  RW with drift:     ", round(rmse(as.numeric(Y_test), rw_fc), 3), "\n")
cat("  Unpooled:          ", round(rmse(as.numeric(Y_test), unpooled_fc), 3), "\n")
cat("  Ratio unpooled 3:  ", round(rmse(as.numeric(Y_test), ratio_unpooled_fc), 3), "\n")
cat("  Pooled 3:          ", round(rmse(as.numeric(Y_test), pooled_fc), 3), "\n")

######## forecast visualization — save one figure per state ##########

years_test <- 2005:2019
context_years <- 1990:2004
years_context <- as.numeric(context_years)

for (plot_state in my_states) {

    actual <- as.numeric(Y_test[plot_state, ])
    actual_context <- as.numeric(Y_full[plot_state, paste0(context_years)])

    fc_rw      <- as.numeric(rw_fc[plot_state, ])
    fc_unpooled <- as.numeric(unpooled_fc[plot_state, ])
    fc_ratio   <- as.numeric(ratio_unpooled_fc[plot_state, ])
    fc_pooled  <- as.numeric(pooled_fc[plot_state, ])

    ylim <- range(c(actual_context, actual, fc_rw, fc_unpooled,
                    fc_ratio, fc_pooled))
    xlim <- range(c(years_context, years_test))

    pdf(file.path("results", paste0("forecast_", plot_state, ".pdf")),
        width = 8, height = 5)
    par(mar = c(4, 4, 3, 1))
    plot(years_context, actual_context, type = "l", col = "black", lwd = 2,
         xlim = xlim, ylim = ylim,
         xlab = "Year", ylab = "Life expectancy (e0)",
         main = paste0(plot_state, ": Actual vs Forecasted e0"))
    lines(years_test, actual, col = "black", lwd = 2)
    points(years_test, actual, col = "black", pch = 16, cex = 0.8)
    abline(v = 2004.5, lty = 2, col = "gray50")

    lines(years_test, fc_rw, col = "blue", lwd = 1.5, lty = 2)
    lines(years_test, fc_unpooled, col = "red", lwd = 1.5, lty = 2)
    lines(years_test, fc_ratio, col = "orange", lwd = 1.5, lty = 2)
    lines(years_test, fc_pooled, col = "darkgreen", lwd = 1.5, lty = 2)

    legend("topleft", bty = "n", cex = 0.8,
           legend = c("Actual",
                      paste0("RW drift (RMSE=",
                             round(rmse(actual, fc_rw), 2), ")"),
                      paste0("Unpooled (RMSE=",
                             round(rmse(actual, fc_unpooled), 2), ")"),
                      paste0("Ratio unpooled (RMSE=",
                             round(rmse(actual, fc_ratio), 2), ")"),
                      paste0("Pooled (RMSE=",
                             round(rmse(actual, fc_pooled), 2), ")")),
           col = c("black", "blue", "red", "orange", "darkgreen"),
           lwd = c(2, 1.5, 1.5, 1.5, 1.5),
           lty = c(1, 2, 2, 2, 2))
    dev.off()
    cat("Saved:", file.path("results", paste0("forecast_", plot_state, ".pdf")), "\n")
}

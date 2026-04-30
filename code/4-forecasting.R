### Rolling-window forecasting analysis
###
### For each origin year, fit four models on data through the origin and
### forecast h-steps ahead. Compare forecasts to actuals across origins
### and horizons.
###
### Models compared:
###   a. Random walk with drift (per region)
###   b. Unpooled MARSS, free sd_ratio (per region)
###   c. Unpooled MARSS, sd_ratio = 3.5 (per region)
###   d. Pooled MARSS, sd_ratio = 3.5 (joint over regions)

library(here)
library(MARSS)

source(here("code", "model-functions.R"))

### 0. Load data -----------------------------------------------------------

Y_full <- read.csv(here("data", "combined_e0.csv"),
                   header = TRUE, row.names = 1, check.names = FALSE)
sigma_full <- read.csv(here("results", "sampling_sds.csv"),
                       header = TRUE, row.names = 1, check.names = FALSE)

my_states <- intersect(rownames(Y_full), rownames(sigma_full))

all_years <- 1959:2019
Y <- as.matrix(Y_full[my_states, paste0(all_years)])
sigma.vec <- setNames(sigma_full[my_states, "x"], my_states)

n_regions <- nrow(Y)

### 1. Helpers -------------------------------------------------------------

## random walk with drift (single-region)
rw_drift_fit_fun <- function(y, sigma, max_iter = 500) {
  if (is.null(dim(y))) {
    y <- matrix(y, nrow = 1)
  } else {
    y <- as.matrix(y)
  }
  n <- nrow(y)
  model_list <- list(
    Z = diag(n),
    B = diag(n),
    Q = "diagonal and equal",
    U = matrix(paste0("u", 1:n), n, 1),
    R = diag(as.numeric(sigma)^2, n),
    A = "zero",
    x0 = matrix(y[, 1], n, 1)
  )
  MARSS(y, model = model_list, method = "BFGS",
        control = list(maxit = max_iter))
}

## h-step forecast from a fit (single- or multi-region)
forecast_fit <- function(fit, h, region_names) {
  pred <- predict(fit, type = "ytT", n.ahead = h)
  out <- matrix(NA, length(region_names), h,
                dimnames = list(region_names, NULL))
  if (length(region_names) == 1) {
    out[1, ] <- tail(pred$pred$estimate, h)
  } else {
    for (r in region_names) {
      sub <- pred$pred[pred$pred$.rownames == r, ]
      out[r, ] <- tail(sub$estimate, h)
    }
  }
  out
}

rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2,
                                              na.rm = TRUE))

### 2. Rolling-window setup ------------------------------------------------

## Use origins so each gets at least max_h years of holdout
max_h        <- 10                       # max forecast horizon
min_train    <- 25                       # min training length
last_year    <- max(all_years)
origin_years <- (min(all_years) + min_train - 1):(last_year - 1)
## predicting at least 1 step beyond origin
origin_years <- origin_years[origin_years <= last_year - 1]
## subsample to 6 evenly spaced windows for testing
origin_years <- origin_years[round(seq(1, length(origin_years), length.out = 6))]

cat("Rolling forecast: origins", min(origin_years), "to", max(origin_years),
    "(", length(origin_years), "windows), horizons 1 to", max_h, "\n")

## storage: list of data frames, one row per (origin, horizon, region, model)
results <- list()

### 3. Loop over origins ---------------------------------------------------

for (origin in origin_years) {

  train_years <- paste0(min(all_years):origin)
  Y_train <- Y[, train_years, drop = FALSE]

  h <- min(max_h, last_year - origin)
  test_years <- paste0((origin + 1):(origin + h))
  Y_test <- Y[, test_years, drop = FALSE]

  cat(sprintf("origin=%d  train_len=%d  h=%d\n",
              origin, length(train_years), h))

  ## fits ----------------------------------------------------------------
  ## a. RW drift, per region
  rw_fc <- matrix(NA, n_regions, h, dimnames = list(my_states, test_years))
  for (i in seq_len(n_regions)) {
    fit <- try(rw_drift_fit_fun(Y_train[i, ], sigma = sigma.vec[i]),
               silent = TRUE)
    if (!inherits(fit, "try-error")) {
      rw_fc[i, ] <- forecast_fit(fit, h, my_states[i])
    }
  }

  ## b. unpooled, free ratio
  unp_fc <- matrix(NA, n_regions, h, dimnames = list(my_states, test_years))
  for (i in seq_len(n_regions)) {
    fit <- try(marss_fit_for_pooled(Y_train[i, ],
                                    sigma = sigma.vec[i],
                                    sd_ratio = NULL),
               silent = TRUE)
    if (!inherits(fit, "try-error")) {
      unp_fc[i, ] <- forecast_fit(fit, h, my_states[i])
    }
  }

  ## c. unpooled, ratio = 3
  ratio_fc <- matrix(NA, n_regions, h, dimnames = list(my_states, test_years))
  for (i in seq_len(n_regions)) {
    fit <- try(marss_fit_for_pooled(Y_train[i, ],
                                    sigma = sigma.vec[i],
                                    sd_ratio = 3.5),
               silent = TRUE)
    if (!inherits(fit, "try-error")) {
      ratio_fc[i, ] <- forecast_fit(fit, h, my_states[i])
    }
  }

  ## d. pooled, ratio = 3
  pooled_fc <- matrix(NA, n_regions, h, dimnames = list(my_states, test_years))
  fit <- try(marss_fit_for_pooled(Y_train,
                                  sigma = sigma.vec,
                                  sd_ratio = 3.5),
             silent = TRUE)
  if (!inherits(fit, "try-error")) {
    pooled_fc[, ] <- forecast_fit(fit, h, my_states)
  }

  ## record ---------------------------------------------------------------
  for (model_name in c("rw_drift", "unpooled", "ratio_unpooled_3",
                       "pooled_3")) {
    fc_mat <- switch(model_name,
                     rw_drift         = rw_fc,
                     unpooled         = unp_fc,
                     ratio_unpooled_3 = ratio_fc,
                     pooled_3         = pooled_fc)
    df <- data.frame(
      origin   = origin,
      horizon  = rep(seq_len(h), each = n_regions),
      region   = rep(my_states, times = h),
      model    = model_name,
      forecast = as.numeric(fc_mat),
      actual   = as.numeric(Y_test),
      stringsAsFactors = FALSE
    )
    results[[length(results) + 1]] <- df
  }
}

forecast_df <- do.call(rbind, results)
forecast_df$error  <- forecast_df$forecast - forecast_df$actual
forecast_df$sq_err <- forecast_df$error^2

### 4. Summaries -----------------------------------------------------------

## RMSE by model and horizon (averaged over origins and regions)
rmse_by_horizon <- aggregate(sq_err ~ model + horizon, data = forecast_df,
                             FUN = function(x) sqrt(mean(x, na.rm = TRUE)))
names(rmse_by_horizon)[3] <- "rmse"
rmse_by_horizon_wide <- reshape(rmse_by_horizon,
                                idvar = "horizon",
                                timevar = "model",
                                direction = "wide")
cat("\n=== RMSE by horizon (avg over origins and regions) ===\n")
print(round(rmse_by_horizon_wide, 3))

## RMSE by model, overall
rmse_overall <- aggregate(sq_err ~ model, data = forecast_df,
                          FUN = function(x) sqrt(mean(x, na.rm = TRUE)))
names(rmse_overall)[2] <- "rmse"
cat("\n=== Overall RMSE by model ===\n")
rmse_overall$rmse <- round(rmse_overall$rmse, 6)
print(rmse_overall)

## RMSE by model and region (avg over origins, horizons)
rmse_by_region <- aggregate(sq_err ~ model + region, data = forecast_df,
                            FUN = function(x) sqrt(mean(x, na.rm = TRUE)))
names(rmse_by_region)[3] <- "rmse"

### 5. Save ----------------------------------------------------------------

saveRDS(list(
  forecast_df       = forecast_df,
  rmse_by_horizon   = rmse_by_horizon_wide,
  rmse_overall      = rmse_overall,
  rmse_by_region    = rmse_by_region,
  origin_years      = origin_years,
  max_h             = max_h,
  my_states         = my_states
), here("results", "rolling_forecasts.rds"))

write.csv(forecast_df,
          here("results", "rolling_forecasts.csv"),
          row.names = FALSE)

cat("\nSaved results to results/rolling_forecasts.rds and .csv\n")

### 6. Quick plot: RMSE vs horizon ----------------------------------------

pdf(here("results", "rolling_rmse_by_horizon.pdf"),
    width = 7, height = 5)
par(mar = c(4, 4, 3, 1))
models <- unique(rmse_by_horizon$model)
cols <- c(rw_drift = "blue", unpooled = "red",
          ratio_unpooled_3 = "orange", pooled_3 = "darkgreen")
ylim <- range(rmse_by_horizon$rmse, na.rm = TRUE)
plot(NA, xlim = c(1, max_h), ylim = ylim,
     xlab = "Forecast horizon (years)", ylab = "RMSE",
     main = "Rolling-window forecast RMSE by horizon")
for (m in models) {
  sub <- rmse_by_horizon[rmse_by_horizon$model == m, ]
  sub <- sub[order(sub$horizon), ]
  lines(sub$horizon, sub$rmse, col = cols[m], lwd = 2)
  points(sub$horizon, sub$rmse, col = cols[m], pch = 16)
}
legend("topleft", legend = models, col = cols[models],
       lwd = 2, pch = 16, bty = "n")
dev.off()
cat("Saved plot to results/rolling_rmse_by_horizon.pdf\n")

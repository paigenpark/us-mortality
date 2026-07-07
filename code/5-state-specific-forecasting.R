### State-specific forecast figures (single train/test split)
###
### Self-contained: fits four models on a single training window
### (1959-2004), forecasts the test period (2005-2019), then saves one
### actual-vs-forecast figure per state. This reconstructs the flat
### per-state forecast matrices that the plotting section needs; the
### rolling-window analysis lives in 4-forecasting.R.

library(here)
library(MARSS)

source(here("code", "model-functions.R"))   # marss_fit_for_pooled()

### 0. Load data -----------------------------------------------------------

Y_full <- read.csv(here("data", "combined_e0.csv"),
                   header = TRUE, row.names = 1, check.names = FALSE)
sigma_full <- read.csv(here("results", "sampling_sds.csv"),
                       header = TRUE, row.names = 1, check.names = FALSE)

my_states <- intersect(rownames(Y_full), rownames(sigma_full))

all_years <- 1959:2019
Y         <- as.matrix(Y_full[my_states, paste0(all_years)])
sigma.vec <- setNames(sigma_full[my_states, "x"], my_states)
n_regions <- nrow(Y)

### 1. Helper functions ----------------------------------------------------
## (mirrors 4-forecasting.R so this script runs standalone)

## random walk with drift (single-region)
rw_drift_fit_fun <- function(y, sigma, max_iter = 500) {
  if (is.null(dim(y))) y <- matrix(y, nrow = 1) else y <- as.matrix(y)
  n <- nrow(y)
  model_list <- list(
    Z = diag(n), B = diag(n),
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
  out  <- matrix(NA, length(region_names), h,
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

rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2, na.rm = TRUE))
}

### 2. Single train/test split ---------------------------------------------

train_years <- 1959:2004
years_test  <- 2005:2019
h           <- length(years_test)

Y_train <- Y[, paste0(train_years), drop = FALSE]
Y_test  <- Y[, paste0(years_test),  drop = FALSE]

### 3. Fit models and forecast ---------------------------------------------

## a. RW drift, per region
rw_fc <- matrix(NA, n_regions, h, dimnames = list(my_states, years_test))
for (i in seq_len(n_regions)) {
  fit <- try(rw_drift_fit_fun(Y_train[i, ], sigma = sigma.vec[i]),
             silent = TRUE)
  if (!inherits(fit, "try-error")) {
    rw_fc[i, ] <- forecast_fit(fit, h, my_states[i])
  }
}

## b. unpooled MARSS, free sd_ratio, per region
unpooled_fc <- matrix(NA, n_regions, h, dimnames = list(my_states, years_test))
for (i in seq_len(n_regions)) {
  fit <- try(marss_fit_for_pooled(Y_train[i, ], sigma = sigma.vec[i],
                                  sd_ratio = NULL),
             silent = TRUE)
  if (!inherits(fit, "try-error")) {
    unpooled_fc[i, ] <- forecast_fit(fit, h, my_states[i])
  }
}

## c. unpooled MARSS, sd_ratio = 3.5, per region
ratio_unpooled_fc <- matrix(NA, n_regions, h,
                            dimnames = list(my_states, years_test))
for (i in seq_len(n_regions)) {
  fit <- try(marss_fit_for_pooled(Y_train[i, ], sigma = sigma.vec[i],
                                  sd_ratio = 3.5),
             silent = TRUE)
  if (!inherits(fit, "try-error")) {
    ratio_unpooled_fc[i, ] <- forecast_fit(fit, h, my_states[i])
  }
}

## d. pooled MARSS, sd_ratio = 3.5, joint over regions
pooled_fc <- matrix(NA, n_regions, h, dimnames = list(my_states, years_test))
fit <- try(marss_fit_for_pooled(Y_train, sigma = sigma.vec, sd_ratio = 3.5),
           silent = TRUE)
if (!inherits(fit, "try-error")) {
  pooled_fc[, ] <- forecast_fit(fit, h, my_states)
}

######## forecast visualization — save one figure per state ##########

context_years <- 1980:2004
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

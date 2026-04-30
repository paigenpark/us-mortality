######## forecast visualization — save one figure per state ##########

years_test <- 2005:2019
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

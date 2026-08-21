results <- readRDS("results_5d.rds")
true_vals <- c(rep(0, 5), 1, 0, 0, 0)

mses <- matrix(NA, nrow=9, ncol=4)
for(p in 1:9) {
  for(m in 1:4) {
    mses[p, m] <- mean((results[, p, m] - true_vals[p])^2, na.rm=TRUE)
  }
}

pdf("figures/mse_comparison_5d.pdf", width=12, height=12)
par(mfrow=c(3, 3), mar=c(5, 5, 4, 1))

titles <- c(expression(mu[1]), expression(mu[2]), expression(mu[3]), 
            expression(mu[4]), expression(mu[5]), expression(sigma^2),
            expression(mu[1] - 2*mu[2]), expression(mu[4] - mu[5]), expression(mu[1] - mu[2] - mu[4]))

methods <- c("IPW", "FL(Unres)", "FL(CorRes)", "FL(IncRes)")
colors <- c("gray", "skyblue", "blue", "red")

for(i in 1:9) {
  bp <- barplot(mses[i, ], names.arg=methods, col=colors,
                main=titles[i], ylab="MSE", las=2, cex.names=1.2, cex.main=2)
  # add text
  text(bp, 0, formatC(mses[i, ], format="e", digits=2), pos=3, cex=1.2, xpd=TRUE)
}

dev.off()
cat("Plot saved to figures/mse_comparison_5d.pdf\n")

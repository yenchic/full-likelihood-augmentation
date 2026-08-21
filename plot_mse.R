results <- read.table("simulation_results.txt", header=TRUE)

# True parameters
n0 <- 50
n1 <- 5000
mu_true <- 0
sigma2_true <- 1
beta_true <- 2
alpha_true <- log(n1 / n0) - mu_true * beta_true - 0.5 * beta_true^2 * sigma2_true

mse <- function(est, truth) {
  mean((est - truth)^2, na.rm = TRUE)
}

mses <- data.frame(
  Parameter = c("mu", "sigma2", "alpha", "beta"),
  IPW = c(
    mse(results[, "mu_IPW"], mu_true),
    mse(results[, "sigma2_IPW"], sigma2_true),
    mse(results[, "alpha_IPW"], alpha_true),
    mse(results[, "beta_IPW"], beta_true)
  ),
  FL = c(
    mse(results[, "mu_FL"], mu_true),
    mse(results[, "sigma2_FL"], sigma2_true),
    mse(results[, "alpha_FL"], alpha_true),
    mse(results[, "beta_FL"], beta_true)
  )
)

# Function to generate plot
generate_plot <- function() {
  par(mfrow=c(2,2), mar=c(3, 4, 3, 1), oma=c(0,0,2,0))
  
  params <- c("mu", "sigma2", "alpha", "beta")
  titles <- c(expression(mu), expression(sigma^2), expression(alpha), expression(beta))
  
  for(i in 1:4) {
    p <- params[i]
    vals <- c(mses$IPW[i], mses$FL[i])
    bp <- barplot(vals, names.arg=c("IPW", "FL"), 
                  col=c("#1f77b4", "#e377c2"), 
                  main=titles[i], ylab="Mean Squared Error",
                  ylim=c(0, max(vals) * 1.25),
                  cex.names=1.2, cex.axis=1.1, cex.lab=1.1, cex.main=1.5)
    
    # Add text labels on top of bars
    text(bp, vals, formatC(vals, format="e", digits=2), pos=3, cex=1.1)
  }
  mtext("MSE Comparison: IPW vs FL Estimators", outer=TRUE, cex=1.5, font=2)
}

# Save as PDF for paper
pdf("figures/mse_comparison.pdf", width=8, height=6)
generate_plot()
dev.off()

# Save as PNG for quick viewing
png("mse_comparison.png", width=800, height=600, res=100)
generate_plot()
dev.off()

cat("Plots generated: figures/mse_comparison.pdf and mse_comparison.png\n")

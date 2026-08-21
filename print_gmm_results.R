res <- readRDS("results_gmm_scaling.rds")
n1_vals <- names(res)

cat("IPW MSE:\n")
ipw_mat <- do.call(rbind, lapply(res, function(x) x$IPW_mean))
colnames(ipw_mat) <- c("mu1", "mu2", "var1", "var2", "w1", "alpha", "beta")
print(round(ipw_mat, 5))

cat("\nFL MSE:\n")
fl_mat <- do.call(rbind, lapply(res, function(x) x$FL_mean))
colnames(fl_mat) <- c("mu1", "mu2", "var1", "var2", "w1", "alpha", "beta")
print(round(fl_mat, 5))

cat("\nImprovement (IPW/FL):\n")
print(round(ipw_mat / fl_mat, 2))

# Save a plot for var1 and w1
pdf("figures/mse_gmm_scaling.pdf", width=10, height=5)
par(mfrow=c(1,2))

n1_num <- as.numeric(n1_vals)

plot(n1_num, ipw_mat[, "var1"], type="b", col="gray", pch=16, lwd=2, 
     ylim=c(0, max(ipw_mat[, "var1"], fl_mat[, "var1"])), 
     xlab="Auxiliary Sample Size (n1)", ylab="MSE of var1", main="MSE vs n1 (GMM var1)")
lines(n1_num, fl_mat[, "var1"], type="b", col="blue", pch=16, lwd=2)
legend("topright", legend=c("IPW", "FL"), col=c("gray", "blue"), lwd=2, pch=16)

plot(n1_num, ipw_mat[, "w1"], type="b", col="gray", pch=16, lwd=2, 
     ylim=c(0, max(ipw_mat[, "w1"], fl_mat[, "w1"])), 
     xlab="Auxiliary Sample Size (n1)", ylab="MSE of w1", main="MSE vs n1 (GMM w1)")
lines(n1_num, fl_mat[, "w1"], type="b", col="blue", pch=16, lwd=2)
legend("topright", legend=c("IPW", "FL"), col=c("gray", "blue"), lwd=2, pch=16)

dev.off()


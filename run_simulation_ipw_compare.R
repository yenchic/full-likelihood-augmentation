set.seed(123)

n_iter <- 1000
n0 <- 50
n1 <- 5000

# True parameters
mu_true <- 0
sigma2_true <- 1
beta_true <- 2
alpha_true <- log(n1 / n0) - mu_true * beta_true - 0.5 * beta_true^2 * sigma2_true

results <- matrix(NA, nrow = n_iter, ncol = 12)
colnames(results) <- c(
  "mu_oldIPW", "sigma2_oldIPW", "alpha_oldIPW", "beta_oldIPW",
  "mu_newIPW", "sigma2_newIPW", "alpha_newIPW", "beta_newIPW",
  "mu_FL", "sigma2_FL", "alpha_FL", "beta_FL"
)

# Negative log-likelihood function for Full-Likelihood (FL)
nll_fl <- function(params, x, a) {
  mu <- params[1]
  log_sigma2 <- params[2]
  sigma2 <- exp(log_sigma2)
  alpha <- params[3]
  beta <- params[4]

  if (sigma2 < 1e-6) {
    return(Inf)
  }

  log_phi <- -0.5 * log(2 * pi * sigma2) - (x - mu)^2 / (2 * sigma2)
  tilt <- a * (alpha + beta * x)

  K <- alpha + mu * beta + 0.5 * beta^2 * sigma2
  log_part <- ifelse(K > 20, K, log1p(exp(K)))

  -sum(log_phi + tilt - log_part)
}

cat("Running simulation with", n_iter, "iterations...\n")

for (i in 1:n_iter) {
  # 1. Data Generation
  x_0 <- rnorm(n0, mean = mu_true, sd = sqrt(sigma2_true))
  a_0 <- rep(0, n0)

  x_1 <- rnorm(n1, mean = mu_true + beta_true * sigma2_true, sd = sqrt(sigma2_true))
  a_1 <- rep(1, n1)

  x <- c(x_0, x_1)
  a <- c(a_0, a_1)

  # 2. Odds model (Logistic Regression)
  glm_fit <- suppressWarnings(glm(a ~ x, family = binomial))
  alpha_ipw <- coef(glm_fit)[1]
  beta_ipw <- coef(glm_fit)[2]
  odds_hat <- exp(alpha_ipw + beta_ipw * x)

  # 3. Old IPW (Target + Reweighted Auxiliary)
  w_old <- (1 - a) + a / odds_hat
  mu_oldIPW <- sum(w_old * x) / sum(w_old)
  sigma2_oldIPW <- sum(w_old * (x - mu_oldIPW)^2) / sum(w_old)

  # 4. New IPW (Reweighted Auxiliary Only)
  w_new <- a / odds_hat
  mu_newIPW <- sum(w_new * x) / sum(w_new)
  sigma2_newIPW <- sum(w_new * (x - mu_newIPW)^2) / sum(w_new)

  # 5. Full-Likelihood (FL) Method
  start_vals <- c(mu_oldIPW, log(sigma2_oldIPW), alpha_ipw, beta_ipw)
  opt <- optim(par = start_vals, fn = nll_fl, x = x, a = a, method = "BFGS")
  mu_fl <- opt$par[1]
  sigma2_fl <- exp(opt$par[2])
  alpha_fl <- opt$par[3]
  beta_fl <- opt$par[4]

  results[i, ] <- c(
    mu_oldIPW, sigma2_oldIPW, alpha_ipw, beta_ipw,
    mu_newIPW, sigma2_newIPW, alpha_ipw, beta_ipw,
    mu_fl, sigma2_fl, alpha_fl, beta_fl
  )

  if (i %% 100 == 0) cat("Completed iteration", i, "\n")
}

write.table(results, "simulation_results_ipw_compare.txt", row.names = FALSE, sep = "\t")
cat("\nResults saved to simulation_results_ipw_compare.txt\n")

# Calculate MSE for each method
mse_oldIPW <- colMeans(sweep(results[, 1:4], 2, c(mu_true, sigma2_true, alpha_true, beta_true))^2, na.rm = TRUE)
mse_newIPW <- colMeans(sweep(results[, 5:8], 2, c(mu_true, sigma2_true, alpha_true, beta_true))^2, na.rm = TRUE)
mse_FL <- colMeans(sweep(results[, 9:12], 2, c(mu_true, sigma2_true, alpha_true, beta_true))^2, na.rm = TRUE)

cat("\n--- Mean Squared Error (MSE) ---\n")
cat(sprintf("%-15s | %-12s | %-12s | %-12s | %-12s\n", "Method", "mu", "sigma2", "alpha", "beta"))
cat("----------------------------------------------------------------------\n")
cat(sprintf("%-15s | %12.5f | %12.5f | %12.5f | %12.5f\n", "old IPW", mse_oldIPW[1], mse_oldIPW[2], mse_oldIPW[3], mse_oldIPW[4]))
cat(sprintf("%-15s | %12.5f | %12.5f | %12.5f | %12.5f\n", "new IPW", mse_newIPW[1], mse_newIPW[2], mse_newIPW[3], mse_newIPW[4]))
cat(sprintf("%-15s | %12.5f | %12.5f | %12.5f | %12.5f\n", "FL", mse_FL[1], mse_FL[2], mse_FL[3], mse_FL[4]))

set.seed(123)

n_iter <- 1000
n0 <- 50
n1 <- 5000

# True parameters
mu_true <- 0
sigma2_true <- 1
beta_true <- 2
alpha_true <- log(n1 / n0) - mu_true * beta_true - 0.5 * beta_true^2 * sigma2_true

results <- matrix(NA, nrow = n_iter, ncol = 8)
colnames(results) <- c("mu_IPW", "sigma2_IPW", "alpha_IPW", "beta_IPW",
                       "mu_FL", "sigma2_FL", "alpha_FL", "beta_FL")

# Negative log-likelihood function for Full-Likelihood (FL)
nll_fl <- function(params, x, a) {
  mu <- params[1]
  log_sigma2 <- params[2]
  sigma2 <- exp(log_sigma2)
  alpha <- params[3]
  beta <- params[4]
  
  # Avoid non-positive variance issues in optimization
  if(sigma2 < 1e-6) return(Inf)
  
  # log phi(x; mu, sigma^2)
  log_phi <- -0.5 * log(2 * pi * sigma2) - (x - mu)^2 / (2 * sigma2)
  
  # tilt applied only when observed
  tilt <- a * (alpha + beta * x)
  
  # log partition function
  K <- alpha + mu * beta + 0.5 * beta^2 * sigma2
  # robust log(1 + exp(K)) computation
  log_part <- ifelse(K > 20, K, log1p(exp(K)))
  
  -sum(log_phi + tilt - log_part)
}

cat("Running simulation with", n_iter, "iterations...\n")

for (i in 1:n_iter) {
  # 1. Data Generation
  # We fix sample sizes n0 and n1 exactly as requested.
  # Target sample (A = 0)
  x_0 <- rnorm(n0, mean = mu_true, sd = sqrt(sigma2_true))
  a_0 <- rep(0, n0)
  
  # Auxiliary sample (A = 1)
  # Derived mathematically in the plan: X|A=1 ~ N(mu + beta*sigma^2, sigma^2)
  x_1 <- rnorm(n1, mean = mu_true + beta_true * sigma2_true, sd = sqrt(sigma2_true))
  a_1 <- rep(1, n1)
  
  x <- c(x_0, x_1)
  a <- c(a_0, a_1)
  
  # 2. IPW Method
  # Fit logistic regression for the odds model
  glm_fit <- glm(a ~ x, family = binomial)
  alpha_ipw <- coef(glm_fit)[1]
  beta_ipw <- coef(glm_fit)[2]
  
  # Compute odds and weights
  # Note: A = 1 for auxiliary, A = 0 for target.
  # The weight formulation is w_i = (1 - A_i) + A_i / O_hat(X_i)
  odds_hat <- exp(alpha_ipw + beta_ipw * x)
  w <- (1 - a) + a / odds_hat
  
  # IPW Weighted MLE for Gaussian
  mu_ipw <- sum(w * x) / sum(w)
  sigma2_ipw <- sum(w * (x - mu_ipw)^2) / sum(w)
  
  # 3. Full-Likelihood (FL) Method
  # Provide IPW estimates as starting values for stability
  # Use log(sigma2) for unconstrained optimization of variance
  start_vals <- c(mu_ipw, log(sigma2_ipw), alpha_ipw, beta_ipw)
  
  opt <- optim(par = start_vals, fn = nll_fl, x = x, a = a, method = "BFGS")
  
  mu_fl <- opt$par[1]
  sigma2_fl <- exp(opt$par[2])
  alpha_fl <- opt$par[3]
  beta_fl <- opt$par[4]
  
  # Save results
  results[i, ] <- c(mu_ipw, sigma2_ipw, alpha_ipw, beta_ipw, mu_fl, sigma2_fl, alpha_fl, beta_fl)
  
  if (i %% 100 == 0) {
    cat("Completed iteration", i, "\n")
  }
}

# Write raw results to txt file
write.table(results, "simulation_results.txt", row.names = FALSE, sep = "\t")

cat("\n--- Simulation Complete ---\n")
cat("Results saved to simulation_results.txt\n\n")

# Calculate MSEs
mse <- function(est, truth) {
  mean((est - truth)^2)
}

mses <- data.frame(
  Parameter = c("mu", "sigma2", "alpha", "beta"),
  IPW_MSE = c(
    mse(results[, "mu_IPW"], mu_true),
    mse(results[, "sigma2_IPW"], sigma2_true),
    mse(results[, "alpha_IPW"], alpha_true),
    mse(results[, "beta_IPW"], beta_true)
  ),
  FL_MSE = c(
    mse(results[, "mu_FL"], mu_true),
    mse(results[, "sigma2_FL"], sigma2_true),
    mse(results[, "alpha_FL"], alpha_true),
    mse(results[, "beta_FL"], beta_true)
  )
)

print(mses)

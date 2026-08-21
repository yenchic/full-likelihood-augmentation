set.seed(456)

n_iter <- 1000
n0 <- 50
n1 <- 5000

# True parameters
beta_true <- 2

# Marginal probability inversion for Unif[-1, 1]
# \int_{-1}^1 0.5 * exp(alpha + beta*x) dx = (n1 / n0)
# exp(alpha) * (exp(beta) - exp(-beta)) / (2 * beta) = n1 / n0
# exp(alpha) = (n1 / n0) * 2 * beta / (exp(beta) - exp(-beta))
alpha_true <- log((n1 / n0) * 2 * beta_true / (exp(beta_true) - exp(-beta_true)))

results <- matrix(NA, nrow = n_iter, ncol = 4)
colnames(results) <- c("alpha_IPW", "beta_IPW", "alpha_FL", "beta_FL")

# Negative log-likelihood function for Full-Likelihood (FL) 
# Note: FL *incorrectly* assumes the target distribution is Gaussian!
nll_fl <- function(params, x, a) {
  mu <- params[1]
  log_sigma2 <- params[2]
  sigma2 <- exp(log_sigma2)
  alpha <- params[3]
  beta <- params[4]
  
  if(sigma2 < 1e-6) return(Inf)
  
  # Gaussian log-density (misspecified model)
  log_phi <- -0.5 * log(2 * pi * sigma2) - (x - mu)^2 / (2 * sigma2)
  tilt <- a * (alpha + beta * x)
  K <- alpha + mu * beta + 0.5 * beta^2 * sigma2
  log_part <- ifelse(K > 20, K, log1p(exp(K)))
  
  -sum(log_phi + tilt - log_part)
}

cat("Running misspecified simulation with", n_iter, "iterations...\n")

for (i in 1:n_iter) {
  # 1. Data Generation
  # Target sample (A = 0): Uniform[-1, 1]
  x_0 <- runif(n0, min = -1, max = 1)
  a_0 <- rep(0, n0)
  
  # Auxiliary sample (A = 1): p(x|A=1) \propto I(x \in [-1, 1]) * exp(beta * x)
  # Inverse transform sampling for this truncated exponential:
  u <- runif(n1)
  x_1 <- (1 / beta_true) * log(u * (exp(beta_true) - exp(-beta_true)) + exp(-beta_true))
  a_1 <- rep(1, n1)
  
  x <- c(x_0, x_1)
  a <- c(a_0, a_1)
  
  # 2. IPW Method
  # IPW only fits the odds model, which is correctly specified!
  glm_fit <- glm(a ~ x, family = binomial)
  alpha_ipw <- coef(glm_fit)[1]
  beta_ipw <- coef(glm_fit)[2]
  
  # 3. Full-Likelihood (FL) Method
  # FL tries to fit the Gaussian models to the uniform data
  start_vals <- c(0, log(1), alpha_ipw, beta_ipw)
  opt <- optim(par = start_vals, fn = nll_fl, x = x, a = a, method = "BFGS")
  
  alpha_fl <- opt$par[3]
  beta_fl <- opt$par[4]
  
  # Save results
  results[i, ] <- c(alpha_ipw, beta_ipw, alpha_fl, beta_fl)
  
  if (i %% 100 == 0) {
    cat("Completed iteration", i, "\n")
  }
}

write.table(results, "simulation_results_mis.txt", row.names = FALSE, sep = "\t")

cat("\n--- Misspecified Simulation Complete ---\n")

mse <- function(est, truth) {
  mean((est - truth)^2)
}

mses <- data.frame(
  Parameter = c("alpha", "beta"),
  IPW_MSE = c(
    mse(results[, "alpha_IPW"], alpha_true),
    mse(results[, "beta_IPW"], beta_true)
  ),
  FL_MSE = c(
    mse(results[, "alpha_FL"], alpha_true),
    mse(results[, "beta_FL"], beta_true)
  )
)

print(mses)

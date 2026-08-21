library(mvtnorm)
set.seed(123)

n_iter <- 500
n0 <- 100
n1 <- 10000

mu0 <- rep(0, 5)
sigma2_0 <- 1
Sigma0 <- sigma2_0 * diag(5)
beta_true <- c(2, 1, 0, 1, 1)

S <- exp(4)
alpha_true <- log(n1 / n0) - log(S)

mu_star <- mu0 + as.vector(Sigma0 %*% beta_true)
Sigma_star <- Sigma0

results <- array(NA, dim=c(n_iter, 9, 4))
dimnames(results)[[2]] <- c("mu1", "mu2", "mu3", "mu4", "mu5", "sigma2", "c1", "c2", "c3")
dimnames(results)[[3]] <- c("IPW", "FL_unres", "FL_c_res", "FL_i_res")

nll_fl_unres <- function(params, x_mat, a_vec) {
  mu <- params[1:5]
  sigma2 <- params[6]
  if(sigma2 <= 1e-5) return(1e100)
  alpha <- params[7]
  beta <- params[8:12]
  
  Sigma <- sigma2 * diag(5)
  log_p_x <- dmvnorm(x_mat, mean=mu, sigma=Sigma, log=TRUE)
  tilt <- a_vec * (alpha + as.vector(x_mat %*% beta))
  
  beta_var_beta <- as.numeric(t(beta) %*% Sigma %*% beta)
  K <- alpha + sum(mu * beta) + 0.5 * beta_var_beta
  log_part <- ifelse(K > 20, K, log1p(exp(K)))
  
  val <- -sum(log_p_x + tilt - log_part)
  if (is.na(val) || is.infinite(val)) return(1e100)
  val
}

nll_fl_c_res <- function(params, x_mat, a_vec) {
  mu <- params[1:5]
  sigma2 <- params[6]
  if(sigma2 <= 1e-5) return(1e100)
  alpha <- params[7]
  b_scale <- params[8]
  beta <- b_scale * c(2, 1, 0, 1, 1)
  
  Sigma <- sigma2 * diag(5)
  log_p_x <- dmvnorm(x_mat, mean=mu, sigma=Sigma, log=TRUE)
  tilt <- a_vec * (alpha + as.vector(x_mat %*% beta))
  
  beta_var_beta <- as.numeric(t(beta) %*% Sigma %*% beta)
  K <- alpha + sum(mu * beta) + 0.5 * beta_var_beta
  log_part <- ifelse(K > 20, K, log1p(exp(K)))
  
  val <- -sum(log_p_x + tilt - log_part)
  if (is.na(val) || is.infinite(val)) return(1e100)
  val
}

nll_fl_i_res <- function(params, x_mat, a_vec) {
  mu <- params[1:5]
  sigma2 <- params[6]
  if(sigma2 <= 1e-5) return(1e100)
  alpha <- params[7]
  b_scale <- params[8]
  beta <- b_scale * c(1, 1, 1, 1, 1)
  
  Sigma <- sigma2 * diag(5)
  log_p_x <- dmvnorm(x_mat, mean=mu, sigma=Sigma, log=TRUE)
  tilt <- a_vec * (alpha + as.vector(x_mat %*% beta))
  
  beta_var_beta <- as.numeric(t(beta) %*% Sigma %*% beta)
  K <- alpha + sum(mu * beta) + 0.5 * beta_var_beta
  log_part <- ifelse(K > 20, K, log1p(exp(K)))
  
  val <- -sum(log_p_x + tilt - log_part)
  if (is.na(val) || is.infinite(val)) return(1e100)
  val
}

start_unres <- c(mu0, sigma2_0, alpha_true, beta_true)
start_c_res <- c(mu0, sigma2_0, alpha_true, 1)
start_i_res <- c(mu0, sigma2_0, alpha_true, 1)

cat("Running 5D simulation...\n")
for(i in 1:n_iter) {
  if (i %% 50 == 0) cat("Iter", i, "\n")
  x0 <- rmvnorm(n0, mean=mu0, sigma=Sigma0)
  a0 <- rep(0, n0)
  x1 <- rmvnorm(n1, mean=mu_star, sigma=Sigma_star)
  a1 <- rep(1, n1)
  x <- rbind(x0, x1)
  a <- c(a0, a1)
  
  # IPW
  glm_fit <- suppressWarnings(glm(a ~ x, family=binomial))
  alpha_ipw <- coef(glm_fit)[1]
  beta_ipw <- coef(glm_fit)[-1]
  O_hat <- exp(alpha_ipw + as.vector(x %*% beta_ipw))
  O_hat <- pmax(O_hat, 1e-10)
  w_ipw <- (1 - a) + a / O_hat
  
  w_sum <- sum(w_ipw)
  mu_ipw <- colSums(w_ipw * x) / w_sum
  
  centered_x <- sweep(x, 2, mu_ipw)
  var_est <- sum(w_ipw * rowSums(centered_x^2)) / w_sum
  sigma2_ipw <- var_est / 5
  
  c1_ipw <- mu_ipw[1] - 2*mu_ipw[2]
  c2_ipw <- mu_ipw[4] - mu_ipw[5]
  c3_ipw <- mu_ipw[1] - mu_ipw[2] - mu_ipw[4]
  
  results[i, , "IPW"] <- c(mu_ipw, sigma2_ipw, c1_ipw, c2_ipw, c3_ipw)
  
  # FL Unrestricted
  opt_unres <- optim(par=start_unres, fn=nll_fl_unres, x_mat=x, a_vec=a, method="BFGS", control=list(maxit=1000))
  mu_unres <- opt_unres$par[1:5]
  s2_unres <- opt_unres$par[6]
  results[i, , "FL_unres"] <- c(mu_unres, s2_unres, mu_unres[1]-2*mu_unres[2], mu_unres[4]-mu_unres[5], mu_unres[1]-mu_unres[2]-mu_unres[4])
  
  # FL Correct Restricted
  opt_c <- optim(par=start_c_res, fn=nll_fl_c_res, x_mat=x, a_vec=a, method="BFGS", control=list(maxit=1000))
  mu_c <- opt_c$par[1:5]
  s2_c <- opt_c$par[6]
  results[i, , "FL_c_res"] <- c(mu_c, s2_c, mu_c[1]-2*mu_c[2], mu_c[4]-mu_c[5], mu_c[1]-mu_c[2]-mu_c[4])
  
  # FL Incorrect Restricted
  opt_i <- optim(par=start_i_res, fn=nll_fl_i_res, x_mat=x, a_vec=a, method="BFGS", control=list(maxit=1000))
  mu_i <- opt_i$par[1:5]
  s2_i <- opt_i$par[6]
  results[i, , "FL_i_res"] <- c(mu_i, s2_i, mu_i[1]-2*mu_i[2], mu_i[4]-mu_i[5], mu_i[1]-mu_i[2]-mu_i[4])
}

saveRDS(results, "results_5d.rds")
cat("Simulation complete. Saved to results_5d.rds\n")

set.seed(123)

n_iter <- 500
n0 <- 50
n1_vals <- c(5000, 10000, 20000, 50000, 100000)

mu_true <- c(0, 3)
sigma2_true <- c(1, 2)
w_true <- c(0.4, 0.6)
beta_true <- 0.2

get_w <- function(g) {
  p <- exp(g) / (1 + exp(g))
  c(p, 1 - p)
}
g_true <- log(w_true[1] / w_true[2])

nll_ipw <- function(params, x, w_ipw) {
  mu <- params[1:2]
  var <- exp(params[3:4])
  w <- get_w(params[5])
  
  if (any(var < 1e-6) || any(var > 100)) return(1e100)
  
  p_mat <- sapply(1:2, function(k) dnorm(x, mu[k], sqrt(var[k])))
  p_x <- as.vector(p_mat %*% w)
  p_x <- pmax(p_x, 1e-300)
  
  val <- -sum(w_ipw * log(p_x))
  if (is.na(val) || is.infinite(val)) return(1e100)
  val
}

nll_fl <- function(params, x, a) {
  mu <- params[1:2]
  var <- exp(params[3:4])
  w <- get_w(params[5])
  alpha <- params[6]
  beta <- params[7]
  
  if (any(var < 1e-6) || any(var > 100)) return(1e100)
  if (abs(beta) > 20) return(1e100)
  if (abs(alpha) > 50) return(1e100)
  
  p_mat <- sapply(1:2, function(k) dnorm(x, mu[k], sqrt(var[k])))
  p_x <- as.vector(p_mat %*% w)
  p_x <- pmax(p_x, 1e-300)
  
  tilt <- a * (alpha + beta * x)
  
  S_temp <- sum(w * exp(beta * mu + 0.5 * beta^2 * var))
  if (is.na(S_temp) || is.infinite(S_temp) || S_temp <= 0) return(1e100)
  
  K <- alpha + log(S_temp)
  log_part <- ifelse(K > 20, K, log1p(exp(K)))
  
  val <- -sum(log(p_x) + tilt - log_part)
  if (is.na(val) || is.infinite(val)) return(1e100)
  val
}

res_list <- list()

for (n1 in n1_vals) {
  cat("Running n1 =", n1, "\n")
  
  S_parts <- w_true * exp(beta_true * mu_true + 0.5 * beta_true^2 * sigma2_true)
  S <- sum(S_parts)
  alpha_true <- log(n1 / n0) - log(S)
  
  mu_star <- mu_true + beta_true * sigma2_true
  sigma2_star <- sigma2_true
  w_star <- S_parts / S
  
  start_ipw <- c(mu_true, log(sigma2_true), g_true)
  
  mses_ipw <- matrix(NA, nrow=n_iter, ncol=7)
  mses_fl <- matrix(NA, nrow=n_iter, ncol=7)
  
  for (iter in 1:n_iter) {
    comp_0 <- sample(1:2, n0, replace=TRUE, prob=w_true)
    x_0 <- rnorm(n0, mu_true[comp_0], sqrt(sigma2_true[comp_0]))
    a_0 <- rep(0, n0)
    
    comp_1 <- sample(1:2, n1, replace=TRUE, prob=w_star)
    x_1 <- rnorm(n1, mu_star[comp_1], sqrt(sigma2_star[comp_1]))
    a_1 <- rep(1, n1)
    
    x <- c(x_0, x_1)
    a <- c(a_0, a_1)
    
    glm_fit <- suppressWarnings(glm(a ~ x, family=binomial))
    alpha_hat <- unname(coef(glm_fit)[1])
    beta_hat <- unname(coef(glm_fit)[2])
    
    O_hat <- exp(alpha_hat + beta_hat * x)
    O_hat <- pmax(O_hat, 1e-10)
    w_ipw <- (1 - a) + a / O_hat
    
    opt_ipw <- optim(par=start_ipw, fn=nll_ipw, x=x, w_ipw=w_ipw, method="BFGS", control=list(maxit=1000))
    est_ipw <- c(opt_ipw$par[1:2], exp(opt_ipw$par[3:4]), get_w(opt_ipw$par[5])[1], alpha_hat, beta_hat)
    
    start_fl <- c(opt_ipw$par, alpha_hat, beta_hat)
    opt_fl <- optim(par=start_fl, fn=nll_fl, x=x, a=a, method="BFGS", control=list(maxit=1000))
    est_fl <- c(opt_fl$par[1:2], exp(opt_fl$par[3:4]), get_w(opt_fl$par[5])[1], opt_fl$par[6:7])
    
    # Label switching fix for IPW
    if (est_ipw[1] > est_ipw[2]) {
      est_ipw[1:2] <- est_ipw[2:1]
      est_ipw[3:4] <- est_ipw[4:3]
      est_ipw[5] <- 1 - est_ipw[5]
    }
    
    # Label switching fix for FL
    if (est_fl[1] > est_fl[2]) {
      est_fl[1:2] <- est_fl[2:1]
      est_fl[3:4] <- est_fl[4:3]
      est_fl[5] <- 1 - est_fl[5]
    }
    
    true_params <- c(mu_true, sigma2_true, w_true[1], alpha_true, beta_true)
    
    mses_ipw[iter, ] <- (est_ipw - true_params)^2
    mses_fl[iter, ] <- (est_fl - true_params)^2
  }
  
  res_list[[as.character(n1)]] <- list(
    IPW_mean = colMeans(mses_ipw, na.rm=TRUE),
    IPW_sd = apply(mses_ipw, 2, sd, na.rm=TRUE),
    FL_mean = colMeans(mses_fl, na.rm=TRUE),
    FL_sd = apply(mses_fl, 2, sd, na.rm=TRUE)
  )
  
  cat("  n1 =", n1, "completed.\n")
}
saveRDS(res_list, "results_gmm_scaling.rds")
print(res_list)

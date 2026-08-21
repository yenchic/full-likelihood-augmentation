set.seed(123)

n_iter <- 500
n0 <- 100
n1_vals <- c(5000, 10000, 20000, 50000, 100000)

run_ising_scaling <- function(N) {
  cat("=== Starting N =", N, "(n0 =", n0, ") ===\n")
  h_true <- rep(0, N)
  num_inter <- N * (N - 1) / 2
  J_true <- rep(0, num_inter)
  if(N==2) J_true <- 1
  if(N==3) J_true <- rep(1, 3)
  if(N==5) J_true[c(1,5,8,10)] <- 1
  
  beta_main <- rep(1, N)
  
  states <- as.matrix(expand.grid(rep(list(c(-1, 1)), N)))
  colnames(states) <- paste0("X", 1:N)
  
  inter_matrix <- matrix(NA, nrow=2^N, ncol=num_inter)
  idx <- 1
  for(i in 1:(N-1)) {
    for(j in (i+1):N) {
      inter_matrix[, idx] <- states[, i] * states[, j]
      idx <- idx + 1
    }
  }
  
  E_0 <- as.vector(states %*% h_true + inter_matrix %*% J_true)
  P_0 <- exp(E_0) / sum(exp(E_0))
  
  res_list <- list()
  
  for (n1 in n1_vals) {
    cat("Running n1 =", n1, "\n")
    alpha_true <- log(n1/n0) - 1.5 
    odds <- exp(alpha_true + as.vector(states %*% beta_main))
    P_1 <- P_0 * odds
    P_1 <- P_1 / sum(P_1)
    
    mses_ipw <- numeric(n_iter)
    mses_fl <- numeric(n_iter)
    
    nll_fl_res <- function(params, counts_0, counts_1) {
      h <- params[1:N]
      J <- params[(N+1):(N+num_inter)]
      alpha <- params[N+num_inter+1]
      b_main <- params[(N+num_inter+2):(2*N+num_inter+1)]
      
      E <- as.vector(states %*% h + inter_matrix %*% J)
      tilt <- alpha + as.vector(states %*% b_main)
      
      max_E <- max(E)
      log_p0 <- E - (max_E + log(sum(exp(E - max_E))))
      
      E_tilt <- E + tilt
      max_E_tilt <- max(E_tilt)
      log_p1 <- E_tilt - (max_E_tilt + log(sum(exp(E_tilt - max_E_tilt))))
      
      -sum(counts_0 * log_p0) - sum(counts_1 * log_p1)
    }
    
    for(iter in 1:n_iter) {
      samp0 <- sample(1:(2^N), n0, replace=TRUE, prob=P_0)
      samp1 <- sample(1:(2^N), n1, replace=TRUE, prob=P_1)
      
      counts_0 <- numeric(2^N)
      counts_1 <- numeric(2^N)
      for(i in 1:(2^N)) {
        counts_0[i] <- sum(samp0 == i)
        counts_1[i] <- sum(samp1 == i)
      }
      
      df0 <- data.frame(state_idx = samp0, a = 0)
      df1 <- data.frame(state_idx = samp1, a = 1)
      df <- rbind(df0, df1)
      for(v in 1:N) {
        df[[paste0("X", v)]] <- states[df$state_idx, v]
      }
      
      form <- as.formula(paste("a ~", paste(paste0("X", 1:N), collapse="+")))
      glm_fit <- suppressWarnings(glm(form, data=df, family=binomial))
      preds <- predict(glm_fit, type="response")
      preds <- pmax(pmin(preds, 1 - 1e-10), 1e-10)
      w <- ifelse(df$a == 0, 1, (1-preds)/preds)
      
      w_counts <- as.numeric(tapply(w, factor(df$state_idx, levels=1:(2^N)), sum))
      w_counts[is.na(w_counts)] <- 0
      
      nll_ipw_fast <- function(params) {
        h <- params[1:N]
        J <- params[(N+1):(N+num_inter)]
        E <- as.vector(states %*% h + inter_matrix %*% J)
        max_E <- max(E)
        log_p <- E - (max_E + log(sum(exp(E - max_E))))
        -sum(w_counts * log_p)
      }
      
      start_ipw <- c(h_true, J_true)
      opt_ipw <- try(optim(start_ipw, nll_ipw_fast, method="BFGS"), silent=TRUE)
      if(inherits(opt_ipw, "try-error")) {
        opt_ipw <- optim(start_ipw, nll_ipw_fast, method="Nelder-Mead")
      }
      mses_ipw[iter] <- mean((opt_ipw$par[(N+1):(N+num_inter)] - J_true)^2)
      
      start_res <- c(h_true, J_true, alpha_true, beta_main)
      opt_res <- try(optim(start_res, nll_fl_res, counts_0=counts_0, counts_1=counts_1, method="BFGS"), silent=TRUE)
      if(inherits(opt_res, "try-error")) {
        opt_res <- optim(start_res, nll_fl_res, counts_0=counts_0, counts_1=counts_1, method="Nelder-Mead")
      }
      mses_fl[iter] <- mean((opt_res$par[(N+1):(N+num_inter)] - J_true)^2)
    }
    
    res_list[[as.character(n1)]] <- c(IPW_mean=mean(mses_ipw), IPW_sd=sd(mses_ipw), FL_mean=mean(mses_fl), FL_sd=sd(mses_fl))
  }
  
  res_df <- do.call(rbind, res_list)
  saveRDS(res_df, paste0("results_ising_scaling_N", N, "_n0_", n0, ".rds"))
  print(res_df)
  
}


set.seed(123)
run_ising_scaling(2)
  # Note: to get the standard error, use sd/sqrt(n_iter)

run_ising_scaling(3)

run_ising_scaling(5)

library(ggplot2)
library(dplyr)
library(gridExtra)

set.seed(123)

n_iter <- 500

# True parameters for the Gaussian Mixture Model
mu_true <- c(0, 3)
sigma2_true <- c(1, 2)
w_true <- c(0.4, 0.6)
beta_true <- 0.2

# True target density
p_target <- function(x) {
  w_true[1] * dnorm(x, mu_true[1], sqrt(sigma2_true[1])) +
  w_true[2] * dnorm(x, mu_true[2], sqrt(sigma2_true[2]))
}

# Auxiliary distribution is also a GMM
mu_aux <- mu_true + beta_true * sigma2_true
w_aux_unnorm <- w_true * exp(beta_true * mu_true + 0.5 * beta_true^2 * sigma2_true)
w_aux <- w_aux_unnorm / sum(w_aux_unnorm)

# Grid for MISE evaluation
x_grid <- seq(-3, 10, length.out = 1000)
dx <- x_grid[2] - x_grid[1]
true_pdf_grid <- p_target(x_grid)

# --- Function to run a single iteration for given n0, n1 ---
run_iter <- function(n0, n1) {
  # 1. Sample target
  comp_0 <- sample(1:2, n0, replace = TRUE, prob = w_true)
  X0 <- rnorm(n0, mu_true[comp_0], sqrt(sigma2_true[comp_0]))
  
  # 2. Sample auxiliary
  comp_1 <- sample(1:2, n1, replace = TRUE, prob = w_aux)
  X1 <- rnorm(n1, mu_aux[comp_1], sqrt(sigma2_true[comp_1]))
  
  # 3. Fit odds model (Logistic Regression)
  X_all <- c(X0, X1)
  A_all <- c(rep(0, n0), rep(1, n1))
  fit <- suppressWarnings(glm(A_all ~ X_all, family = binomial()))
  
  alpha_hat <- coef(fit)[1]
  beta_hat <- coef(fit)[2]
  
  # 4. Target KDE
  kde_target <- density(X0, bw = "nrd0")
  pdf_target_hat <- approx(kde_target$x, kde_target$y, xout = x_grid, rule = 2)$y
  
  # 5. ASR KDE
  kde_aux <- density(X1, bw = "nrd0")
  pdf_aux_hat <- approx(kde_aux$x, kde_aux$y, xout = x_grid, rule = 2)$y
  
  O_hat_grid <- exp(alpha_hat + beta_hat * x_grid)
  
  # Normalizing constant (sample ratio estimator)
  Omega_hat <- n0 / n1
  
  # ASR-KDE carefully dividing by O(x)
  pdf_asr_hat <- (pdf_aux_hat / O_hat_grid) / Omega_hat
  
  # 6. Calculate MSE integral (MISE)
  err_target <- sum((pdf_target_hat - true_pdf_grid)^2 * dx)
  err_asr <- sum((pdf_asr_hat - true_pdf_grid)^2 * dx)
  
  # IPW-KDE
  W_ipw <- 1 / (n0 * exp(alpha_hat + beta_hat * X1))
  kde_ipw <- suppressWarnings(density(X1, weights = W_ipw, bw = kde_aux$bw))
  pdf_ipw_hat <- approx(kde_ipw$x, kde_ipw$y, xout = x_grid, rule = 2)$y
  
  # AIPW-KDE
  pdf_aipw_hat <- 0.5 * pdf_target_hat + 0.5 * pdf_ipw_hat
  
  err_ipw <- sum((pdf_ipw_hat - true_pdf_grid)^2 * dx, na.rm=TRUE)
  err_aipw <- sum((pdf_aipw_hat - true_pdf_grid)^2 * dx, na.rm=TRUE)
  
  c(err_target, err_asr, err_ipw, err_aipw)
}

# --- Batch 1: Scale n0, Fixed n1 = 5000 ---
n0_vals <- c(50, 100, 200, 400, 800, 1600)
n1_fixed <- 50000
mise_target_b1 <- numeric(length(n0_vals))
mise_asr_b1 <- numeric(length(n0_vals))
mise_ipw_b1 <- numeric(length(n0_vals))
mise_aipw_b1 <- numeric(length(n0_vals))
se_target_b1 <- numeric(length(n0_vals))
se_asr_b1 <- numeric(length(n0_vals))
se_ipw_b1 <- numeric(length(n0_vals))
se_aipw_b1 <- numeric(length(n0_vals))

cat("Starting Batch 1...\n")
for (k in 1:length(n0_vals)) {
  n0 <- n0_vals[k]
  errs <- simplify2array(parallel::mclapply(1:n_iter, function(i) run_iter(n0, n1_fixed), mc.cores = 8))
  mise_target_b1[k] <- mean(errs[1, ])
  mise_asr_b1[k] <- mean(errs[2, ])
  mise_ipw_b1[k] <- mean(errs[3, ])
  mise_aipw_b1[k] <- mean(errs[4, ])
  se_target_b1[k] <- sd(errs[1, ]) / sqrt(n_iter)
  se_asr_b1[k] <- sd(errs[2, ]) / sqrt(n_iter)
  se_ipw_b1[k] <- sd(errs[3, ]) / sqrt(n_iter)
  se_aipw_b1[k] <- sd(errs[4, ]) / sqrt(n_iter)
  cat(sprintf("n0 = %4d | Target MISE: %.5f | ASR MISE: %.5f\n", n0, mise_target_b1[k], mise_asr_b1[k]))
}

df_b1 <- data.frame(
  n0 = rep(n0_vals, 4),
  MISE = c(mise_target_b1, mise_asr_b1, mise_ipw_b1, mise_aipw_b1),
  ymin = c(mise_target_b1 - 1.96 * se_target_b1, mise_asr_b1 - 1.96 * se_asr_b1, mise_ipw_b1 - 1.96 * se_ipw_b1, mise_aipw_b1 - 1.96 * se_aipw_b1),
  ymax = c(mise_target_b1 + 1.96 * se_target_b1, mise_asr_b1 + 1.96 * se_asr_b1, mise_ipw_b1 + 1.96 * se_ipw_b1, mise_aipw_b1 + 1.96 * se_aipw_b1),
  Method = rep(c("Target-KDE", "ASR-KDE", "IPW-KDE", "TAIPW-KDE"), each=length(n0_vals))
)

# --- Batch 2: Scale n1, Fixed n0 = 50 ---
n0_fixed <- 50
n1_vals <- c(5000, 10000, 20000, 50000, 100000)
mise_target_b2 <- numeric(length(n1_vals))
mise_asr_b2 <- numeric(length(n1_vals))
mise_ipw_b2 <- numeric(length(n1_vals))
mise_aipw_b2 <- numeric(length(n1_vals))
se_target_b2 <- numeric(length(n1_vals))
se_asr_b2 <- numeric(length(n1_vals))
se_ipw_b2 <- numeric(length(n1_vals))
se_aipw_b2 <- numeric(length(n1_vals))

cat("\nStarting Batch 2...\n")
for (k in 1:length(n1_vals)) {
  n1 <- n1_vals[k]
  errs <- simplify2array(parallel::mclapply(1:n_iter, function(i) run_iter(n0_fixed, n1), mc.cores = 8))
  mise_target_b2[k] <- mean(errs[1, ])
  mise_asr_b2[k] <- mean(errs[2, ])
  mise_ipw_b2[k] <- mean(errs[3, ])
  mise_aipw_b2[k] <- mean(errs[4, ])
  se_target_b2[k] <- sd(errs[1, ]) / sqrt(n_iter)
  se_asr_b2[k] <- sd(errs[2, ]) / sqrt(n_iter)
  se_ipw_b2[k] <- sd(errs[3, ]) / sqrt(n_iter)
  se_aipw_b2[k] <- sd(errs[4, ]) / sqrt(n_iter)
  cat(sprintf("n1 = %6d | Target MISE: %.5f | ASR MISE: %.5f\n", n1, mise_target_b2[k], mise_asr_b2[k]))
}

df_b2 <- data.frame(
  n1 = rep(n1_vals, 4),
  MISE = c(mise_target_b2, mise_asr_b2, mise_ipw_b2, mise_aipw_b2),
  ymin = c(mise_target_b2 - 1.96 * se_target_b2, mise_asr_b2 - 1.96 * se_asr_b2, mise_ipw_b2 - 1.96 * se_ipw_b2, mise_aipw_b2 - 1.96 * se_aipw_b2),
  ymax = c(mise_target_b2 + 1.96 * se_target_b2, mise_asr_b2 + 1.96 * se_asr_b2, mise_ipw_b2 + 1.96 * se_ipw_b2, mise_aipw_b2 + 1.96 * se_aipw_b2),
  Method = rep(c("Target-KDE", "ASR-KDE", "IPW-KDE", "TAIPW-KDE"), each=length(n1_vals))
)


# --- Plotting ---

p1 <- ggplot(df_b1, aes(x = n0, y = MISE, color = Method, shape = Method, group = Method)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.05, linewidth = 0.8) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_x_log10(breaks = n0_vals) +
  scale_y_log10() +
  labs(title = expression("Batch 1: Increasing " * n[0] * " (Fixed " * n[1] * " = 50000)"),
       x = expression(n[0]), y = "MISE") +
  theme_bw() +
  theme(legend.position = "bottom",
        title = element_text(size=12),
        legend.text = element_text(size=10),
        legend.title = element_blank())

p2 <- ggplot(df_b2, aes(x = n1, y = MISE, color = Method, shape = Method, group = Method)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.05, linewidth = 0.8) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_x_log10(breaks = n1_vals) +
  scale_y_log10() +
  labs(title = expression("Batch 2: Increasing " * n[1] * " (Fixed " * n[0] * " = 50)"),
       x = expression(n[1]), y = "MISE") +
  theme_bw() +
  theme(legend.position = "bottom",
        title = element_text(size=12),
        legend.text = element_text(size=10),
        legend.title = element_blank())

pdf("figures/mse_asr_scaling.pdf", width = 10, height = 5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

cat("\nDone! Plot saved to figures/mse_asr_scaling.pdf\n")

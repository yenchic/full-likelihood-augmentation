# Reproducibility Scripts

This directory contains the R scripts necessary to reproduce the simulations, tables, and figures presented in the paper 
> *On efficiency gains via augmenting a tiny sample with a massive auxiliary data*. 

## 1. Gaussian-Logistic Model (Motivating Example)
These scripts correspond to the baseline simulation in Section 2.1 (Figure 1), demonstrating the full efficiency gain in a simple univariate Gaussian target model with a logistic odds model.
- `run_simulation.R`: Runs the Monte Carlo simulation comparing the Inverse Probability Weighting (IPW) and Full-Likelihood (FL) estimators.
- `plot_mse.R`: Processes the output of the simulation and generates the MSE comparison plot.

## 2. Ising Model (Appendix)
These scripts correspond to the complex structural invariance scaling simulations across dimensions ($N=2, 3, 5$). 
- `run_simulation_ising.R`: This unified script encapsulates the full logic for generating the different network topologies (complete graphs for $N=2,3$ and a sparse chain for $N=5$). You can manually adjust the target sample size `n0` at the top of the file (e.g., set `n0 <- 100` or `n0 <- 200`) to replicate the different bottleneck constraints shown in the paper's scaling plots.

## 3. Gaussian Mixture Model (Appendix)
These scripts correspond to the 2-component GMM simulations demonstrating full efficiency gains for structurally invariant parameters (variances), and bounded efficiency for confounded parameters.
- `run_simulation_gmm_scaling.R`: Runs the GMM Monte Carlo simulation across increasing auxiliary sample sizes ($n_1$).
- `print_gmm_results.R`: Processes the GMM simulation results and formats them into the exact tables presented in the appendix.

## 4. 5D Gaussian Model
These scripts correspond to the multivariate extension of the baseline model, demonstrating how efficiency gains persist across multiple dimensions for structurally invariant parameters.
- `run_simulation_5d.R`: Runs the 5-dimensional Gaussian Monte Carlo simulation.
- `plot_mse_5d.R`: Generates the MSE comparison plots across the multiple dimensions.

## 5. Model Mis-specification
These scripts correspond to the robustness analysis, demonstrating how the IPW and FL estimators behave when the target distribution is deliberately mis-specified.
- `run_simulation_mis.R`: Runs the Monte Carlo simulation under model mis-specification.
- `plot_mse_mis.R`: Generates the MSE comparison plots highlighting the robustness trade-offs between IPW and FL.

## 6. Nonparametric Density Estimation (KDE)
These scripts correspond to the nonparametric density estimation experiments, illustrating the efficiency gains of the ASR-KDE and comparing it against IPW-KDE and TAIPW-KDE estimators under varying sample size scaling regimes.
- `run_simulation_asr_kde.R`: Runs the Monte Carlo simulation to evaluate the Mean Integrated Squared Error (MISE) of the KDE estimators, and generates the corresponding MSE scaling plots.
- `run_simulation_ipw_compare.R`: Runs a targeted simulation to compare the standard old IPW estimator, the new IPW estimator, and the FL estimator across multiple parameters. 
- `plot_mse_compare.R`: Processes the output of the IPW comparison simulation and generates the MSE comparison plots.

library(quantreg)
library(zoo)
library(stringr)
library(stats)

options(warn = -1)

# Find last significant lag in [1, dq_max]
QADF_lag_selection <- function(variable, tau, dq_max) {
  for (dq in dq_max:1) {
    lag_formula_str <- paste0(
      "L(variable, ", 1:dq_max, ")", 
      collapse = " + ") 
    model <- quantreg::dynrq(
      formula(str_glue("variable ~ trend(variable) +
                           {lag_formula_str}")),
      tau = tau)
    coeff <- coefficients(summary(model, se = "boot"))
    # If last coefficient significant, choose it
    if (coeff[dq + 2, 4] < 0.05) {
      return(dq)
    }
  }
  return(1)
}


## ADF test statistic, (Koenker, Xiao, 2004)
QADF <- function(variable, taus, dq_max) {
  t_stats <- c()
  alphas_hat <- c()
  t_stats <- c()
  Y_length <- length(variable)
  DF_stats <- c(-2.86, -2.81, -2.75, -2.72, -2.64, -2.58, -2.51, -2.40, -2.28, -2.12)
  for (tau_i in seq_along(taus)) {
    tau = taus[tau_i]
    print(paste("Working with tau =", tau))
    # Select lags
    lag_max <- max(2, QADF_lag_selection(variable, tau, dq_max))
    
    ## Estimate quantile density function with kernel
    # Get limit bandwidth, h -> 0 as n -> infty
    # We use bandwidth.rq from quantreg package
    # It estimates Hall-Sheather kernel bandwidth by default
    h_bandwidth <- quantreg::bandwidth.rq(tau, n = length(variable))
    tau_lower <- tau - h_bandwidth
    tau_upper <- tau + h_bandwidth
    
    # Generate Q(tau + h), Q(tau - h)
    model_lower <- quantreg::dynrq(
      variable ~ trend(variable) + L(variable, 1) + L(d(variable), 1:lag_max),
      tau = tau_lower
    )
    
    model_upper <- quantreg::dynrq(
      variable ~ trend(variable) + L(variable, 1) + L(d(variable), 1:lag_max),
      tau = tau_upper
    )
    
    # This is a central model Q(tau)
    # Use this later
    model <- quantreg::dynrq(
      variable ~ 
        trend(variable) + 
        L(variable, 1) + 
        L(d(variable), 1:lag_max),
      tau = tau)
    
    # Extract coefficients, generate Q(tau + h) - Q(tau - h)
    coeff_upper <- coefficients(model_upper)
    coeff_lower <- coefficients(model_lower)
    coeff_diff <- matrix(coeff_upper - coeff_lower)
    
    # Zero padded variable series
    Y_lag <- na.fill(stats::lag(variable, -1, na.pad = TRUE), 0)
    Y_diff_lags <- na.fill(stats::lag(diff(zoo(variable)), -(1:lag_max)), 0)

    # Matrix with mean of variables
    X <- cbind(zoo(1),
               mean(seq(1:Y_length)),
               mean(Y_lag),
               t(colMeans(Y_diff_lags)))

    # Q(tau + h) - Q(tau - h) = mean(x) @ (alpha(tau + h) - alpha(tau - h))
    X_diff <- X %*% (coeff_diff)
    
    # Quantile density
    # NB!: fofinvQ > 1, by def.
    fofinvF <- (2 * h_bandwidth  / X_diff)
    
    # Generate projection matrix to space X = (1, L(d(y), 1), L(d(y), 2), ...)
    X_space <- cbind(zoo(rep(1, Y_length)), Y_diff_lags)
    X_space <- na.fill(X_space, 0)
    X_space <- data.matrix(as.data.frame(X_space))
    
    # Projection matrix
    PX <- diag(Y_length) - (X_space %*%
                              solve(t(X_space) %*% X_space) %*%
                              t(X_space))
    
    # Extract persistence coefficient
    alpha_hat <- coefficients(model)[3]
    
    # Test alpha = 1
    t_stat <- ((fofinvF / sqrt(tau * (1 - tau))) *
                 sqrt(t(Y_lag) %*% t(PX) %*% Y_lag) *
                 (alpha_hat - 1))
    alphas_hat <- c(alphas_hat, alpha_hat)
    t_stats <- c(t_stats, t_stat)
  }
  t_stats_significant <- which(t_stats < DF_stats[tau_i])
  alphas_hat[t_stats_significant] <- paste0(alphas_hat[t_stats_significant], "*")
  return(data.frame(alpha = alphas_hat))
}
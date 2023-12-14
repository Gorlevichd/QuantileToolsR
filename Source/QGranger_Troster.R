### Quantile Granger Test module
# Done through F-test
library(quantreg)


QGranger_stat <- function(X, Y, q, p, taus) {
  ## Weight matrix
  # It can be calculated without iterating over tau,
  # so we do it outside tau loop 
  # Get information as lags of Y and X
  start_id <- start(Y)[1][1] + q
  T <- length(Y) - q
  Info_t <- window(cbind(stats::lag(zoo(Y), -(1:q)), stats::lag(zoo(X), -(1:p))), 
                    start = start_id)
  W_non_exp <- matrix(0, nrow = nrow(Info_t), ncol = nrow(Info_t))
  for (i in seq_along(colnames(Info_t))) {
    # Take separate column from Info_t
    # Get W for column, update W with sum
    Info_t_col <- Info_t[, i]
    flattened_Info <- matrix(Info_t_col, nrow = length(Info_t_col))
    ones_col <- matrix(1, ncol = length(Info_t_col))
    repeat_matr <- flattened_Info %*% ones_col
    time_diff <- repeat_matr - t(repeat_matr)
    W_non_exp = W_non_exp + (-0.5 * (time_diff)^2)
  }
  W <- exp(W_non_exp)

  psi_matrix <- matrix(0, length(Y) - q, 1)
  test_stat <- 0
  # Iterate over tau, denoted as j in the paper
  for (j in seq_along(taus)) {
    tau <- taus[j]

    ## Psi-matrix
    ## Def:
    ## psi(i, j) = psi(tau_j, residual_i)
    ## psi(tau_j, eps) = 1[eps <= 0] - tau_j
    ## eps = [Y_i - hat(Q_i)]
    # Fit regression with only lags of Y
    Y_only_fit <- dynlm(Y ~ L(Y, 1:q), tau = tau)
    # Get predicted quantile for all Y
    Q_hat <- predict(Y_only_fit)
    # Get residual
    Q_residual <- window(Y, start = start_id) - Q_hat 
    # Get psi(eps)
    psi_Q_residual_j <- 1 * (Q_residual <= 0) - tau
    psi_matrix[, 1] = psi_Q_residual_j
    tau_stat <- t(psi_matrix) %*% (W / T) %*% psi_matrix
    test_stat = test_stat + abs(tau_stat)
  }
  test_stat <- test_stat / ((length(Y) - q) * length(taus))
  return(test_stat)
}


QGranger_CritVal <- function(X, Y, q, p, taus, k = 3) {
  ## We get the critical values with bootstart 
  # Set k = 3 by default, as they indicated it is a good choice
  # (Troster, 2018); (Chernozhukov, Ferndandez-Val, 2005)
  subsample_size <- floor(k * length(Y)^(2/5))
  n_subsamples <- length(Y) - subsample_size - 1
  stats_distribution <- c()
  for (subsample_i in 2:n_subsamples) {
    subsample_Y <- window(Y, 
                          start = subsample_i, 
                          end = subsample_i + subsample_size - 1)
    subsample_X <- window(X,
                          start = subsample_i,
                          end = subsample_i + subsample_size - 1)
    stats_distribution <- c(stats_distribution, 
                            QGranger_stat(subsample_X, subsample_Y, q, p, taus))
  }
  return(stats_distribution)
}

QGranger <- function(X, Y, q, p, taus, k = 3) {
  Qstat <- QGranger_stat(X, Y, q, p, taus)
  QGranger_distr <- QGranger_CritVal(X, Y, q, p, taus)
  mean(QGranger_distr < Qstat)
}
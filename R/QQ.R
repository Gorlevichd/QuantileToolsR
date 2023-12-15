#' Quantile-in-Quantile Regression
#' 
#' Fit Quantile-in-Quantile Regression with Gaussian kernel
#' @param Y Dependent variable
#' @param X Independent variable
#' @param tau_Y List of qunatiles of Y
#' @param tau_X List of quantiles of X
#' @param h Kernel width
#' @return QQ Result 3d matrix with beta(i, j) s. t. Q(Y, tau(i)) = beta(i, j) Q(X, T(j)) + eps(i, j)$
#' @export
QQRegression <- function(Y, X, tau_Y, tau_X, h = 0.05) {
  X_quantiles <- quantile(X, tau_X)
  X_coef_matrix <- matrix(0, 
                          nrow = length(tau_Y), 
                          ncol = length(tau_X))
  for (tau_Y_i in seq_along(tau_Y)) {
    tau_Y_current <- tau_Y[tau_Y_i]
    for (tau_X_i in seq_along(tau_X)) {
      # This is a quantile value, rather than
      q_X_current <- X_quantiles[tau_X_i]
      distance_X_to_q <- X - q_X_current
      kernel_weight <- dnorm(distance_X_to_q / h)
      qqmodel <- quantreg::dynrq(Y ~ 
                                   L(Y, 1) + 
                                   distance_X_to_q,
                                 tau = tau_Y_current,
                                 weight = window(kernel_weight, start = 2))
      X_coef_matrix[tau_Y_i, tau_X_i] <- coefficients(qqmodel)[3]
    }
  }
  rownames(X_coef_matrix) <- tau_Y
  colnames(X_coef_matrix) <- tau_X 
  return(X_coef_matrix)
}
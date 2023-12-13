### Quantile Granger Test module
# Done through F-test
library(quantreg)

QGranger <- function(X, Y, qar_order_max, taus) {
  # Init empty df 
  F_pval_df <- data.frame(tau = taus)
  # Iterate over [AR(1), ... AR(qar_order_max)]
  for (qar_order in 1:qar_order_max) {
    F_pval_q <- c()
    print(paste0("Working with AR(", qar_order, ")"))
    # Iterate over tau in [taus]
    for (tau_i in seq_along(taus)) {
      tau = taus[tau_i]
      # Get model with X variable and Y lags
      Y_with_X_rq <- quantreg::dynrq(
        d(Y) ~ L(d(Y), 1:qar_order) + L(d(X), 1),
        tau = tau)
      # Model without X, just Y lags
      just_Y_rq <- quantreg::dynrq(
        d(Y) ~ L(d(Y), 1:qar_order),
        tau = tau)
      # H0: X lags do not affect Y
      # Check if just_Y_rq is statistically the same as Y_with_X_rq
      F_pval_tau <- quantreg::anova.rq(
        Y_with_X_rq, 
        just_Y_rq,
        se = "boot"
      )[["table"]]$pvalue
      F_pval_q <- c(F_pval_q, F_pval_tau)
    }
    F_pval_q_df <- data.frame(F_pval_q)
    colnames(F_pval_q_df) <- str_glue("AR({qar_order})")
    F_pval_q_df$tau <- taus
    F_pval_df <- merge(F_pval_df, 
                       F_pval_q_df, 
                       by = "tau", 
                       all = TRUE)
  }
  return(F_pval_df)
}
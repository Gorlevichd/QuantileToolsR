#' Quantile Error Correction model with step selection
#' 
#' Fits QECM and selects lags
#' @param p Dependent variable
#' @param q1 Independent variable 1
#' @param q2 Inpendendt variable 2
#' @param p_lag Maximum lag of p
#' @param q1_lag Maximum lag of q1
#' @param q2_lag Maximum lag of q2
#' @param tau_list List of quantiles for fit
#' @return Dataframe with fitted coefficients
#' @export
QECM <- function(p, q1, q2, 
                 p_lag, q1_lag, q2_lag,
                 tau_list) {
  p_lag_var <- paste0(
    "L(d(p), ", 1:p_lag, ")")
  # q1, q2 - explanatory
  q1_lag_var <- paste0(
    "L(d(q1), ", 1:q1_lag, ")")
  q2_lag_var <- paste0(
    "L(d(q2), ", 1:q2_lag, ")")
  ECM_var <- c("L(coint_eps, 1)", 
               stringr::str_glue("{p_lag_var}"),  
               stringr::str_glue("{q1_lag_var}"),
               stringr::str_glue("{q2_lag_var}"),
               "Intercept_coint", "q1_coint", "q2_coint")
  result_df <- data.frame(var = ECM_var)

  # Iterate over tau
  for (tau in tau_list) {
    print(paste0("Working with tau = ", tau))
    ## This is a two-step ECM approach
    # Step 1: Estimate Long-run cointegrating model, OLS
    coint_model <- quantreg::dynrq(
      p ~ q1 + q2,
      tau = tau)
    # Extract coefficients + bootstrap confidence interval
    coint_result_matr <- coefficients(summary(coint_model, se = "boot"))
    # Get residuals to plug them into Step 2
    coint_eps <- residuals(coint_model)

    # Step 2: Construct ECM model with residuals
    # For stepwise to work, we need to explicitly add all the lags
    # Otw it deletes a whole variable
    # p is target
    p_lag_formula <- paste0(
      "L(d(p),", 1:p_lag, ")", 
      collapse = " + ")
    # q1, q2 - explanatory
    q1_lag_formula <- paste0(
      "L(d(q1),", 1:q1_lag, ")", 
      collapse = " + ")
    q2_lag_formula <- paste0(
      "L(d(q2),", 1:q2_lag, ")",
      collapse = " + ")
    ECM_str_formula <- paste0("d(p) ~ L(coint_eps, 1) + ", 
                              stringr::str_glue(
                              "{p_lag_formula} + 
                              {q1_lag_formula} + 
                              {q2_lag_formula}")
                              )
    
    # Select best model with stepwise AIC
    # start = 10, otw different length because of lags
    ECM_model <- step(quantreg::dynrq(
      formula(ECM_str_formula),
      tau = tau, 
      start = max(p_lag, q1_lag, q2_lag) + 2), trace = 0)
    # Extract coefficients, join with cointeg rating coefficients
    ECM_result_matr <- coefficients(summary(ECM_model, se = "boot"))
    ECM_coeffs <- matrix(ECM_result_matr[, 1])
    coint_coeffs <- matrix(-coint_result_matr[, 1])
    ECM_coeffs <- rbind(ECM_coeffs, coint_coeffs)
    colnames(ECM_coeffs) <- tau

    # Add signifiance symbols to values
    ## 1% = ***
    ## 5% = **
    ## 10% = *
    ECM_pval <- matrix(ECM_result_matr[, 2])
    coint_pval <- matrix(coint_result_matr[, 2])
    ECM_pval <- rbind(ECM_pval, coint_pval)
    ECM_pval_1_perc <- which(ECM_pval < 0.01)
    ECM_pval_5_perc <- which((ECM_pval < 0.05) & ECM_pval >= 0.01)
    ECM_pval_10_perc <- which((ECM_pval < 0.1) & (ECM_pval >= 0.05))
    ECM_coeffs = round(ECM_coeffs, digits = 4)
    ECM_coeffs[ECM_pval_1_perc] <- paste0(ECM_coeffs[ECM_pval_1_perc], "***")
    ECM_coeffs[ECM_pval_5_perc] <- paste0(ECM_coeffs[ECM_pval_5_perc], "**")
    ECM_coeffs[ECM_pval_10_perc] <- paste0(ECM_coeffs[ECM_pval_10_perc], "*")
    coint_terms <- paste0(c("Intercept",
                            attr(coint_model$terms, "term.labels")), 
                          "_coint")
    ECM_df <- data.frame(ECM_coeffs)
    colnames(ECM_df) <- tau
    ECM_df$var <- c("Intercept", 
                        attr(ECM_model$terms , "term.labels"),
                        coint_terms)
    result_df <- merge(result_df, ECM_df, all = TRUE, on = "var")
  }
  return(result_df)
}
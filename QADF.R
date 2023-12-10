library(quantreg)
library(readxl)
library(zoo)
library(stats)
library(xts)

data <- read_excel("EXDATA.xlsx")
min(data$Date)
max(data$Date)
variable <- ts(log(data$ER))
quantile(log(data$ER), c(0.06, 0.1, 0.13))
### QAR ADF Test


## Lag selection
# I select by median AIC, as it seems to better catch the 
# generalized performance of the model
qardf.select_lags <- function(ts, pmax, tau) {
  AIC_best <- Inf
  diff_p_best <- NULL
  for (diff_p in 1:pmax) {
    dfm <- dynrq(ts ~ trend(ts) + L(ts) + 
                      L(diff(ts), 1:diff_p), 
                 tau = tau)
    median_AIC <- median(AIC(dfm))
    if (median_AIC < AIC_best) {
      AIC_best <- median_AIC
      diff_p_best <- diff_p
    }
  }
  return (diff_p_best)
}


## Test Statistic
ts <- variable
pmax <- 7
tau_seq <- seq(0.1, 0.9, 0.1)
tr_v <- 1:length(ts)
hs <- FALSE
tau <- 0.5
diff_AR_order <- qardf.select_lags(ts, pmax, tau_seq)

## Estimate Density Function
# Estimate kernel bandwidth
h_width = bandwidth.rq(tau, length(ts), hs=hs)

# Estimate upper and lower model
tau_lower = tau - h_width
tau_upper = tau + h_width
model <- dynrq(ts ~ trend(ts) + 
                    L(ts) + 
                    L(diff(ts), 1:diff_AR_order), 
               tau = c(tau_lower, tau, tau_upper))
at_means <- data.frame(t(colMeans(model.frame(model))))
predict(model, newdata=data.frame(at_means))
# X_mean <- matrix(colMeans(model$x), nrow = 1)
  # predicted_quant <- X_mean %*% alpha
  # predicted_quant_upper <- X_mean %*% alpha_upper
  # predicted_quant_lower <- X_mean %*% alpha_lower
  # alpha_diff <- (alpha_upper - alpha_lower)
  # denominator <- X_mean %*% (alpha_diff)
  # 
  # f_ofinvF <- (2 * h_width)
  # # Calculate t-stat
  # Reference (Koenker, Xiao, 2004)
  # X <- matrix(X, nrow = length(ts))
  # Px <- X %*% solve(t(X) %*% X) %*% t(X)
  # t_n <- ((f_ofinvF / (sqrt(tau[tau_i] * (1 - tau[tau_i])))) * 
  #           (t(Y_lag) %*% Px %*% Y_lag)^(1/2) * (alpha - 1))
  # print(t_n)
  break
}
library(stringr)
library(plotly)
library(reshape2)
# Module with plotter functions

# Plotly plotter for timeseries
# Plotly kinda looks nicer than ggplot2
# and standard plotting library, 
# so why not

#' Time-series plots with linear trend
#' 
#' Plots plotly graph for y in time with linear trend
#' @param y series
#' @param dates date-index
#' @return plotly line plot
#' @export
ts_plot_ly <- function(y, dates) {
    y_index <- dates
    varnm <- deparse(substitute(y))
    linear_trend <- lm(y ~ dates)
    plotly::plot_ly(type = "scatter", mode = "lines") %>%
      add_trace(y = y, 
                x = dates, 
                name = varnm) %>%
      add_trace(
        y = predict(linear_trend),
        x = dates,
        name = stringr::str_glue("{varnm}, linear trend")) %>%
      layout(
        xaxis = (list(title = "Date")),
        yaxis = (list(title = varnm)))
}


#' Quantile-in-Quantile Plot
#' 
#' Plots Quantile-in-Quantile parameters againts quantiles of $X$ and $Y$
#' @param QQresult QQresult matrix
#' @return Surface Plot with tau_X and tau_Y on x and y axis and corresponding beta
#' @export
QQplot <- function(QQresult) {
    tau_X <- colnames(QQresult)
    tau_Y <- rownames(QQresult)
    plotly::plot_ly() %>%
      add_surface(x = tau_X, 
        y = tau_Y, 
        z = QQresult)
}
library(stringr)
library(plotly)
library(reshape2)
# Module with plotter functions

# Plotly plotter for timeseries
# Plotly kinda looks nicer than ggplot2
# and standard plotting library, 
# so why not
ts_plot_ly <- function(y, dates) {
    y_index <- dates
    varnm <- deparse(substitute(y))
    linear_trend <- lm(y ~ dates)
    plot_ly(type = "scatter", mode = "lines") %>%
    add_trace(y = y, 
            x = dates, 
            name = varnm) %>%
    add_trace(
        y = predict(linear_trend),
        x = dates,
        name = str_glue("{varnm}, linear trend")) %>%
    layout(
        xaxis = (
            list(
                title = "Date"
            )
        ),
        yaxis = (
            list(
                title = varnm
            )
        )
    )
}


QQplot <- function(QQresult) {
    tau_X <- colnames(QQresult)
    tau_Y <- rownames(QQresult)
    plot_ly() %>%
        add_surface(x = tau_X, 
        y = tau_Y, 
        z = QQresult)
}
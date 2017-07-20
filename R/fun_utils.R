#' Generic Plot
#' @keywords internal
#' @export
plot.pclm <- function(x, ...) {
  with(x, {
    logic1  <- is.null(input$offset)
    breaksR <- input$breaks[-1]
    breaksL <- rev(rev(input$breaks)[-1]) + 1
    n       <- length(breaksR)
    leng    <- breaksR - breaksL + 1
    hist    <- input$dta / leng
    Age     <- min(breaksL):max(breaksR) - 1
    # estimated distribution
    log_scale <- if (!logic1) 'y' else ''
    y_lim     <- c(min(fitted.values), 
                   min(max(fitted.values$value*2), max(fitted.values)))
    
    plot(Age, fitted.values$value, type = "n", axes = FALSE,
         ylab = '', xlab = '', log = log_scale, ylim = y_lim,
         main = 'Observed vs. Estimated data', ...)
    axis(1); axis(2)
    # histogram
    if (logic1) {
      for (i in 1:n) {
        polygon(x = c(breaksR[i], breaksL[i], 
                      breaksL[i], breaksR[i]),
                y = c(0, 0, hist[i], hist[i]),
                angle = seq(45, 180, length = n)[i],
                col = '#f7f4f9')
      }
    }
    # confidence intervals
    polygon(x = c(rev(Age), Age), 
            y = c(rev(fitted.values$lower.bound), 
                  fitted.values$upper.bound),
            col = rgb(1, 0, 0, alpha = 0.1), border = "red")
    # Fitted line
    lines(Age, fitted.values$value, lwd = 2, col = "#0c2c84")
    legend('topleft', legend = c('Fitted Values', 'Conf.Intervals'), 
           col = c("#0c2c84", "#fcc5c0"), lty = 1, lwd = c(2, 4), bty = 'n')
  })
}



#' @keywords internal
#' @export
print.pclm <- function(x, ...){
  with(x$input,
       {
         cat('\nPenalized Composite Link Model\n------------------------------')
         cat('\nNumber of input groups:', length(breaks) - 1)
         cat('\nNumber of fitted values:', nrow(x$fitted.values))
         cat('\nEstimated smoothing parameter:', lambda)
         cat('\nAIC:', x$goodness.of.fit$AIC)
         cat('\nBIC:', x$goodness.of.fit$BIC, '\n')
       })
}


#' @keywords internal
#' @export
summary.pclm <- function(object, ...){
  with(object$input,
       {
         cat('\nPenalized Composite Link Model\n------------------------------')
         cat('\nNumber of input groups:', length(breaks) - 1)
         cat('\nNumber of fitted values:', nrow(object$fitted.values))
         cat('\nEstimated smoothing parameter:', lambda)
         cat('\nAIC:', object$goodness.of.fit$AIC)
         cat('\nBIC:', object$goodness.of.fit$BIC, '\n')
       })
}
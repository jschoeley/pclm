#' Generic Plot
#' @keywords internal
#' @export
plot.pclm <- function(x, ...) {
  with(x, {
    breaksR = input$breaks[-1]
    breaksL = rev(rev(input$breaks)[-1]) + 1
    n = length(breaksR)
    leng = breaksR - breaksL + 1
    hist = input$y / leng
    Age  = min(breaksL):max(breaksR) - 1
    # estimated distribution
    plot(Age, fitted.values$value, type = "l", 
         ylab = "Number of deaths", xlab = "Age", 
         main = "Age at death distribution", 
         col = "red", lwd = 3)
    # histogram
    for (i in 1:n) {
     polygon(x = c(breaksR[i], breaksL[i], 
                   breaksL[i], breaksR[i]),
             y = c(0, 0, hist[i], hist[i]),
             angle = seq(45, 180, length = n)[i])
    }
    # confidence intervals
    lines(Age, fitted.values$lower.bound, col = "pink") 
    lines(Age, fitted.values$upper.bound, col = "pink")
  })
}


#' @keywords internal
#' @export
print.pclm <- function(x, ...){
  with(x$input,
       {
         cat('\nPenalized Composite Link Model\n------------------------------')
         cat('\nNumber of input groups:', length(breaks))
         cat('\nNumber of fitted values:', nrow(x$fitted.values))
         cat('\nEstimated smoothing parameter:', lambda)
         cat('\nAIC:', x$goodness.of.fit$AIC)
         cat('\nBIC:', x$goodness.of.fit$BIC)
       })
}


#' @keywords internal
#' @export
summary.pclm <- function(object, ...){
  with(object$input,
       {
         cat('\nPenalized Composite Link Model\n------------------------------')
         cat('\nNumber of input groups:', length(breaks))
         cat('\nNumber of fitted values:', nrow(object$fitted.values))
         cat('\nEstimated smoothing parameter:', lambda)
         cat('\nAIC:', object$goodness.of.fit$AIC)
         cat('\nBIC:', object$goodness.of.fit$BIC)
       })
}
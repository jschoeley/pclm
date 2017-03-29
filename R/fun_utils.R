
#' Generic Plot
#' @keywords internal
#' @export
plot.pclm <- function(x, ...) {
  with(x, {
       n = length(breaksR)
       leng = breaksR - breaksL + 1
       hist = x$input$y / leng
       Age  = min(breaksL):max(breaksR) - 1
       # estimated distribution
       plot(Age, gamma, type = "l", ylab = "Number of deaths", xlab = "Age", 
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
       lines(Age, eta_lower, col = "pink") 
       lines(Age, eta_upper, col = "pink")
  })
}




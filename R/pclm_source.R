
#' Fit Univariate Penalized Composite Link Model
#' 
#' @param y Vector of observed counts
#' @param breaks Breaks
#' @param lambda Smoothing parameter. If \code{lambda = NULL} (default) an 
#' algorithm will find the optimal value.
#' @param deg Order of differences of the components of b
#' @param show a logical value. Indicates whether iteration details should be shown
#' @return A \code{pclm} object
#' @examples  
#' library(pclm)
#' 
#' # Select data
#' country = "SWE"
#' year    = 2014
#' sex     = "Total"
#' 
#' deaths <- hmddx[hmddx$country == country & hmddx$Year == year, ]
#' 
#' # Aggregate them artificially in 5-years age classes with 85+ 
#' # add last bin with 0 counts to complete the histogram
#' breaks = c(seq(5,85,5),115,130)
#' deaths$Groups_Counts <- c(rep(1:17, each = 5), rep(18, 26))
#' aggredate_dx <- aggregate(deaths$Total, FUN = "sum",
#'                           by = list(deaths$Groups_Counts))$x
#' aggredate_dx <-  c(aggredate_dx, 0)
#' 
#' 
#' # Fit model
#' mod = pclm(y = aggredate_dx, breaks)
#' 
#' # Generic plot
#' plot(mod)
#' # Add real data on top
#' lines(0:110, deaths$Total, type = "o", col = "blue") 
#' 
#' @export
pclm <- function(y, breaks, lambda = NULL, deg = 2, show = FALSE){
  
  lambda.hat = lambda
  if (is.null(lambda)) {
  # Find lambda if NULL
    opt <- optimise(f = objective_fun, y = y, breaks = breaks,
                    what = 'aic', interval = c(0, 10^7))
    lambda.hat <- opt$minimum
  }
  # solve the PCLM 
  mod = fit_pclm(y, breaks, lambda = lambda.hat, deg = 2, show)
  cat('lambda.hat, ED & AIC:', lambda.hat, mod$trace, mod$aic, '\n')
  
  out <- structure(class = "pclm", mod)
  return(out)  
}


#' Objective function (to be optimised)
#' @keywords internal
#' 
objective_fun <- function(x, y, breaks, what = 'aic') {
  fit <- fit_pclm(y, breaks, lambda = x, deg = 2)
  if (what == 'aic') out = fit$aic else out = fit$bic
  return(out)
}


#' Fit univariate pclm 
#' 
#' @keywords internal
fit_pclm <- function(y, breaks, lambda = 1, deg = 2, show = FALSE){
  input <- c(as.list(environment())) # save all the input for later use
  breaksR = breaks
  breaksL <- c(1, rev(rev((breaksR + 1))[-1]))
  n = length(y)
  m = rev(breaks)[1]
  C = matrix(0, nrow = n, ncol = m)
  for (i in 1:n) C[i, breaksL[i]:breaksR[i]] =  1
  X = diag(m)
  
  # Some preparations
  nx <- dim(X)[2]
  D  <- diff(diag(nx), diff = deg)
  bstart <- log(sum(y) / nx);
  b <- rep(bstart, nx);
  
  # Perform the iterations
  for (it in 1:50) {
    b0  <- b
    eta <- X %*% b
    gam <- exp(eta)
    mu  <- C %*% gam
    w   <- c(1 / mu, rep(lambda, nx - deg)) 
    Gam <- gam %*% rep(1, nx)
    Q   <- C %*% (Gam * X)
    z   <- c(y - mu + Q %*% b, rep(0, nx - deg))
    Fit <- lsfit(rbind(Q, D), z, wt = w, intercept = F) 
    b   <- Fit$coef
    
    db <- max(abs(b - b0))
    if (show)  cat(it, " ", db, "\n")
    if (db < 1e-6) break
  }
  
  # Regression diagnostics
  R  <- t(Q) %*% diag(c(1 / mu)) %*% Q
  H  <- solve(R + lambda * t(D) %*% D, R)
  H0 <- solve(R + lambda * t(D) %*% D) # variance-covariance matrix Bayesian approach
  H1 <- H0 %*% R %*% H0 # variance-covaraince matrix sandwitch estimator
  
  trace <- sum(diag(H))
  ok <- y > 0
  dev <- 2 * sum(y[ok] * log(y[ok] / mu[ok]))
  aic <- dev + 2 * trace
  bic <- dev + log(length(y)) * trace
  
  # standard errors using sandwitch estimator 
  # (use alternatively mod$H0 for Bayesian approach)
  s.e. = sqrt(diag(H1))
  # confidence intervals
  lower = exp(eta - 2*s.e.)
  upper = exp(eta + 2*s.e.)
  
  out <- list(input = input,
              gamma = gam, mu = mu, H0 = H0, H1 = H1, 
              trace = trace, dev = dev, aic = aic, bic = bic,
              eta = eta, eta_upper = upper, eta_lower = lower, s.e. = s.e.,
              breaksL = breaksL, breaksR = breaksR)
  return(out)
}


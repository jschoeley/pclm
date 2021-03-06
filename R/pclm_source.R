
#' Fit Univariate Penalized Composite Link Model
#' 
#' @param dta Vector of observed counts
#' @param breaks Breaks in the data set. For example: if we have 3 bins 
#' \code{[0-5), [5-10) and [10-15]} the breaks will be defined as a vector 
#' \code{c(0, 5, 10, 16)}
#' @param offset Offset term.
#' @param lambda Smoothing parameter. If \code{lambda = NULL} (default) an 
#' algorithm will find the optimal value.
#' @param show Logical value. Indicates whether iteration details should be shown or not.
#' Default: \code{FALSE}.
#' @param ci.level Level of significance for computing confidence intervals. 
#' Default: \code{0.05}.
#' @param objective.fun Objective function used in optimisation process.
#' Choose between \code{AIC} or \code{BIC}. Default: \code{AIC}
#' @param opt.interval Interval to be used in one-dimensional optimisation of lambda. 
#' Default: \code{c(0, 10^5)}.
#' @param ... Other parameters to pe passed in \code{fit_pclm} internal function.
#' @return A \code{pclm} object.
#' @seealso \code{\link{pclm_mx}}
#' @examples  
#' library(pclm)
#' 
#' # Select data
#' year = 2014
#' sex  = "Total"
#' 
#' deaths <- hmdDx[hmdDx$Year == year, ]
#' 
#' # Aggregate them artificially in 5-years age classes with 85+ 
#' # add last bin with 0 counts to complete the histogram
#' breaks = c(seq(0, 85, by = 5), 115, 130)
#' deaths$Groups_Counts <- c(rep(1:17, each = 5), rep(18, 26))
#' aggredate_dx <- aggregate(deaths$Total, FUN = "sum",
#'                           by = list(deaths$Groups_Counts))$x
#' aggredate_dx <-  c(aggredate_dx, 0)
#' 
#' # Fit model
#' mod = pclm(aggredate_dx, breaks)
#' mod
#' fitted.values(mod)
#' 
#' # Generic plot
#' plot(mod)
#' # Add real data on top
#' lines(0:110, deaths$Total, type = "o", col = "blue") 
#' 
#' @importFrom stats lsfit optimise qnorm
#' @export
#' 
pclm <- function(dta, breaks, offset = NULL, lambda = NULL, 
                 show = FALSE, ci.level = 0.05,
                 objective.fun = 'AIC', opt.interval = c(0, 10^5), ...){
  
  if (is.null(lambda)) {
  # Find lambda if NULL
    opt <- optimise(f = objective_fun, 
                    dta = dta, 
                    breaks = breaks,
                    offset = offset,
                    what = objective.fun, 
                    interval = opt.interval)
    lambda <- opt$minimum
  }
  # solve the PCLM 
  mdl <- fit_pclm(dta, breaks, offset, lambda, show, ci.level)
  out <- structure(class = "pclm", mdl)
  return(out)  
}


#' Objective function (to be optimised)
#' @keywords internal
#' 
objective_fun <- function(x, dta, breaks, offset, what) {
  fit <- fit_pclm(dta, breaks, offset, lambda = x)
  if (what == 'AIC') out = fit$goodness.of.fit$AIC else out = fit$goodness.of.fit$BIC
  return(out)
}

#' Fit univariate pclm 
#' 
#' @keywords internal
fit_pclm <- function(dta, breaks, offset = NULL, lambda, 
                     show = FALSE, ci.level = 0.05,
                     deg = 2, iter = 50, tol = 1e-6){
  input <- c(as.list(environment())) # save all the input for later use
  
  breaksR = breaks[-1]
  breaksL = rev(rev(breaks)[-1]) + 1
  m = max(breaks) - min(breaks)
  X = diag(m)
  
  # Build C matrix
  C = matrix(0, nrow = length(dta), ncol = m)
  for (j in 1:length(dta)) C[j, breaksL[j]:breaksR[j] - min(breaks)] = 1
  if (!is.null(offset)) C <- C %*% diag(c(offset))
  
  # Some preparations
  nx <- dim(X)[2]
  D  <- diff(diag(nx), diff = deg)
  bstart <- log(sum(dta) / nx)
  b  <- rep(bstart, nx)
  
  # Perform the iterations
  for (it in 1:iter) {
    b0  <- b
    eta <- X %*% b
    gam <- exp(eta)
    mu  <- C %*% gam
    w   <- c(1 / mu, rep(lambda, nx - deg)) 
    Gam <- gam %*% rep(1, nx)
    Q   <- C %*% (Gam * X)
    z   <- c(dta - mu + Q %*% b, rep(0, nx - deg))
    Fit <- lsfit(rbind(Q, D), z, wt = w, intercept = F) 
    b   <- Fit$coef
    db  <- max(abs(b - b0))
    
    if (show)  cat(it, " ", db, "\n")
    if (db < tol) break
  }
  
  # Regression diagnostics
  R    <- t(Q) %*% diag(c(1 / mu)) %*% Q
  H    <- solve(R + lambda * t(D) %*% D, R)
  H0   <- solve(R + lambda * t(D) %*% D) # variance-covariance matrix Bayesian approach
  H1   <- H0 %*% R %*% H0 # variance-covaraince matrix sandwitch estimator
  vcov <- list(vcov_bayesian = H0, vcov_sandwitch = H1)
  
  trace <- sum(diag(H))
  ok    <- dta > 0
  dev   <- 2 * sum(dta[ok] * log(dta[ok] / mu[ok]))
  aic   <- dev + 2 * trace
  bic   <- dev + log(length(dta)) * trace
  # standard errors using sandwitch estimator 
  # (use alternatively H0 for Bayesian approach)
  s.e. <- sqrt(diag(H1))
  gof  <- list(AIC = aic, BIC = bic, deviance = dev, standard.errors = s.e.)
  
  # confidence intervals 5%
  qn    <- qnorm(1 - ci.level/2)
  lower <- exp(eta - qn*s.e.)
  upper <- exp(eta + qn*s.e.)
  fit   <- data.frame(value = gam, lower.bound = lower, upper.bound = upper)
  
  out <- list(input = input, fitted.values = fit, 
              goodness.of.fit = gof, lambda = lambda)
  return(out)
}











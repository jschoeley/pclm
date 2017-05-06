
#' @title Estimate death rates using PCLM
#' @description Estimate 1-year age-specific death rates from aggregated data 
#' using the univariate penalized composite link model.
#' @inheritParams pclm
#' @seealso \code{\link{pclm}}
#' @return The function will return 3 pclm objects: one for dta data, 
#' one for data provided in the offset argument, and one for the death rates.
#' The first two pclm objects can be obtained by fitting two separate \code{\link{pclm}}s too, 
#' using the two data sets (dta and offset).
#' 
#' @examples 
#' library(pclm)
#' 
#' # Select data
#' year = 2014
#' sex  = "Total"
#' 
#' exposure <- hmdEx[hmdEx$Year == year, ]
#' deaths <- hmdDx[hmdDx$Year == year, ]
#' 
#' # Aggregate them artificially in 5-years age classes with 85+
#' # add last bin with 0 counts to complete the histogram
#' breaks = c(seq(0, 85, by = 5), 115, 130)
#' exposure$Groups_Counts = deaths$Groups_Counts <- c(rep(1:17, each = 5), rep(18, 26))
#' 
#' aggredate_Dx <- aggregate(deaths$Total, FUN = "sum",
#'                           by = list(deaths$Groups_Counts))$x
#' aggredate_Dx <- c(aggredate_Dx, 0)
#' 
#' aggredate_Ex <- aggregate(exposure$Total, FUN = "sum",
#'                           by = list(exposure$Groups_Counts))$x
#' aggredate_Ex <- c(aggredate_Ex, 0)
#' 
#' # Fit the models ---------------------
#' model <- pclm_mx(dta = aggredate_Dx, offset = aggredate_Ex, breaks)
#' 
#' M1 = model$pclm_dta
#' M2 = model$pclm_offset
#' M3 = model$pclm_mx
#' 
#' # Plot1 - Observed and estimated death distribution
#' par(mfrow = c(1, 3))
#' plot(M1)
#' lines(0:110, deaths$Total, type = "o", col = "blue")
#' 
#' # Plot 2 - Observed and estimated distribution of population exposed to risk
#' plot(M2)
#' lines(0:110, exposure$Total, type = "o", col = "blue")
#' 
#' # Plot 3 - Age-specific death rate estimation
#' plot(M3)
#' mx_observed  <- deaths[, 'Total']/exposure[, 'Total']
#' points(mx_observed)
#' 
#' @export
pclm_mx <- function(dta, offset, breaks, lambda = NULL, deg = 2, 
                    show = FALSE, ci.level = 0.05,
                    objective.fun = 'AIC', opt.interval = c(0, 10^5)){
  
  mdl_Dx = pclm(dta = dta, breaks, offset = NULL, 
                lambda, deg, show, ci.level, 
                objective.fun, opt.interval)
  
  mdl_Ex = pclm(dta = offset, breaks, offset = NULL, 
                lambda, deg, show, ci.level,
                objective.fun, opt.interval)
  
  mdl_mx = pclm(dta = dta, breaks, offset = mdl_Ex$fitted.values$value, 
                lambda, deg, show, ci.level,
                objective.fun, opt.interval)
  
  out <- list(pclm_dta = mdl_Dx, pclm_offset = mdl_Ex, pclm_mx = mdl_mx)
  return(out)
}



% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pclm_source.R
\name{pclm}
\alias{pclm}
\title{Fit Univariate Penalized Composite Link Model}
\usage{
pclm(dta, breaks, offset = NULL, lambda = NULL, show = FALSE,
  ci.level = 0.05, objective.fun = "AIC", opt.interval = c(0, 10^5), ...)
}
\arguments{
\item{dta}{Vector of observed counts}

\item{breaks}{Breaks in the data set. For example: if we have 3 bins 
\code{[0-5), [5-10) and [10-15]} the breaks will be defined as a vector 
\code{c(0, 5, 10, 16)}}

\item{offset}{Offset term.}

\item{lambda}{Smoothing parameter. If \code{lambda = NULL} (default) an 
algorithm will find the optimal value.}

\item{show}{Logical value. Indicates whether iteration details should be shown or not.
Default: \code{FALSE}.}

\item{ci.level}{Level of significance for computing confidence intervals. 
Default: \code{0.05}.}

\item{objective.fun}{Objective function used in optimisation process.
Choose between \code{AIC} or \code{BIC}. Default: \code{AIC}}

\item{opt.interval}{Interval to be used in one-dimensional optimisation of lambda. 
Default: \code{c(0, 10^5)}.}

\item{...}{Other parameters to pe passed in \code{fit_pclm} internal function.}
}
\value{
A \code{pclm} object.
}
\description{
Fit Univariate Penalized Composite Link Model
}
\examples{
 
library(pclm)

# Select data
year = 2014
sex  = "Total"

deaths <- hmdDx[hmdDx$Year == year, ]

# Aggregate them artificially in 5-years age classes with 85+ 
# add last bin with 0 counts to complete the histogram
breaks = c(seq(0, 85, by = 5), 115, 130)
deaths$Groups_Counts <- c(rep(1:17, each = 5), rep(18, 26))
aggredate_dx <- aggregate(deaths$Total, FUN = "sum",
                          by = list(deaths$Groups_Counts))$x
aggredate_dx <-  c(aggredate_dx, 0)

# Fit model
mod = pclm(aggredate_dx, breaks)
mod
fitted.values(mod)

# Generic plot
plot(mod)
# Add real data on top
lines(0:110, deaths$Total, type = "o", col = "blue") 

}
\seealso{
\code{\link{pclm_mx}}
}

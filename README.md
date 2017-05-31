# pclm R Package - work in progress...
[![Build Status](https://travis-ci.org/mpascariu/pclm.svg?branch=master)](https://travis-ci.org/mpascariu/pclm)
[![license](https://img.shields.io/github/license/mpascariu/pclm.svg)]()

This repository contains code for estimating smooth distributions using the penalized composite link methodology. Examples applied to demographic data are provided.

## Installation

1. Make sure you have the most recent version of R
2. Run the following code in your R console 

```r
# install.packages("devtools")
devtools::install_github("mpascariu/pclm")
```


## Abstract
Ungrouping binned data can be desirable for many reasons: Bins can be too coarse to allow for accurate analysis; comparisons can be hindered when different grouping approaches are used in different histograms; and the last interval is often wide and open-ended and, thus, covers a lot of information in the tail area. Age group–specific disease incidence rates and abridged life tables are examples of binned data. We propose a versatile method for ungrouping histograms that assumes that only the underlying distribution is smooth. Because of this modest assumption, the approach is suitable for most applications. The method is based on the composite link model, with a penalty added to ensure the smoothness of the target distribution. Estimates are obtained by maximizing a penalized likelihood. This maximization is performed efficiently by a version of the iteratively reweighted least-squares algorithm. Optimal values of the smoothing parameter are chosen by minimizing Akaike’s Information Criterion. We demonstrate the performance of this method in a simulation study and provide several examples that illustrate the approach. Wide, open-ended intervals can be handled properly. The method can be extended to the estimation of rates when both the event counts and the exposures to risk are grouped.

## Reference
[Silvia Rizzi](http://findresearcher.sdu.dk:8080/portal/en/person/srizzi), [Jutta Gampe](http://www.demogr.mpg.de/en/institute/staff_directory_1899/jutta_gampe_655.htm) and Paul H. C. Eilers - (2015) - Efficient Estimation of Smooth Distributions From Coarsely Grouped Data - Am J Epidemiol  182 (2): 138-147. DOI: https://doi.org/10.1093/aje/kwv020

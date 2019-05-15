
<!-- README.md is generated from README.Rmd. Please edit that file -->
hindex
======

<!-- badges: start -->
<!-- badges: end -->
This package provides functionality to simulate the development of h-index and h-alpha values of scientists who collaborate on writing papers and to visualize the simulated data. The effect of publishing, being cited, and strategic collaborating can be simulated. 

Installation
------------

You can install the released version of hindex from [GitHub](https://github.com) using the devtools package:

``` r
# if devtools is not installed yet:
install.packages("devtools")

devtools::install_github("hindex")
```

Example
-------

``` r
library(hindex)
simdata <- simulate_hindex(runs = 2, n = 20, periods = 3)
plot_hsim(simdata, plot_hindex = TRUE)
```

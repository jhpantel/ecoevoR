
# ecoevor

<!-- badges: start -->
<!-- badges: end -->

ecoevor is a collection of computational and statistical tools to study eco-evolutionary dynamics.

## Installation

You can install the development version of ecoevor from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("jhpantel/ecoevoR")
```
``` r
# install.packages("R.rsp")
# options(timeout=10000) # depends on your download speed, package is 700mb
devtools::install_github("jhpantel/ecoevoR",build_vignettes=TRUE)
```

## Example

This is one example of runnig a simulation of Lotka-Volterra competition dynamics with an evolving trait x that impacts growth:

``` r
library(ecoevor)
sim_num <- 10
m1 <- rand_sim(n=sim_num,x0=c(0.4,0.5,0.6),evol=TRUE)
## Visualize the results
matplot(1:300, m1$pop_record[1,,7:9], type = "l", xlab="time",ylab="N")
```


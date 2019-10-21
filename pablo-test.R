##
## R vs Cpp comparison
## 20-OCT-2019
##

library(Rcpp)
library(microbenchmark)

## --------------------------------------------------------
## TEST: R/rhoOpt.R vs src/rhoOpt.cpp
source("R/rhoOpt.R")
Rcpp::sourceCpp("src/rhoOpt.cpp")

## Both versions return the same values:
xx <- runif(1000,-10, 10)
cc <- runif(1, 1, 5)
all.equal(rhoOpt(xx, cc), rhoOpt_pablo(xx, cc))
all.equal(psiOpt(xx, cc), psiOpt_pablo(xx, cc))
all.equal(derpsiOpt(xx, cc), derpsiOpt_pablo(xx, cc))
## Benchmarking:
microbenchmark(rhoOpt(xx, cc), rhoOpt_pablo(xx, cc))
microbenchmark(psiOpt(xx, cc), psiOpt_pablo(xx, cc))
microbenchmark(derpsiOpt(xx, cc), derpsiOpt_pablo(xx, cc))

## --------------------------------------------------------
## TEST: R/mscale.R vs src/mscale.cpp
source("R/mscale.R")
Rcpp::sourceCpp("src/rhoOpt.cpp")

## Both versions return the same values:
set.seed(1)
xx <- runif(1000,-10, 10)
bb <- 0.5
cc <- runif(1, 1, 5)
all.equal(Mscale(xx, bb, cc), Mscale_pablo(xx, bb, cc))
## Benchmarking:
microbenchmark(Mscale(xx, bb, cc), Mscale_pablo(xx, bb, cc))

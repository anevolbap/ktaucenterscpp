##
## R vs Cpp comparison
## 20-OCT-2019
##

library(Rcpp)
library(microbenchmark)

## --------------------------------------------------------
## TEST: R/rhoOpt.R vs src/rhoOpt.cpp
source("R/rhoOpt.R")
Rcpp::sourceCpp("src/mscale.cpp")

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
source("R/rhoOpt.R")
source("R/mscale.R")
Rcpp::sourceCpp("src/mscale.cpp")

## Both versions return the same values:
xx <- runif(1e2, -10, 10)
bb <- 0.5
cc <- runif(1, 1, 5)
all.equal(Mscale(xx, bb, cc), Mscale_pablo(xx, bb, cc))
all.equal(Mscale(xx, bb, cc), Mscale_pablo_sin(xx, bb, cc))

## Benchmarking:
m <- microbenchmark(Mscale(xx, bb, cc),
                    Mscale_pablo(xx, bb, cc),
                    Mscale_pablo_sin(xx, bb, cc),
                    mad(xx),
                    sd(xx),
                    times = 1e2)
plot(m)

m

Rprof("profileo-pablo.out")
invisible(replicate(2000, Mscale_pablo(xx, bb, cc)))
Rprof(NULL)

summaryRprof("profileo-pablo.out")

##
## R vs Cpp comparison
## 27-JAN-2020
##

library(Rcpp)
library(microbenchmark)


## Identical??
set.seed(5)
Z <- rnorm(600)
mues <- rep(c(-4, 0, 4), 200)
X <-  matrix(Z + mues, ncol = 2)
ktau_cpp <- ktaucenterscpp::ktaucenters(X, K = 4, nstart = 15)

set.seed(5)
Z <- rnorm(600)
mues <- rep(c(-4, 0, 4), 200)
X <-  matrix(Z + mues, ncol = 2)
ktau <- ktaucenters::ktaucenters(X, K = 4, nstart = 15)

identical(ktau, ktau_cpp)
all.equal(ktau, ktau_cpp)

## Performance improved?
start = Sys.time()
for (iter in 1:10){
    set.seed(iter)
    Z <- rnorm(600)
    mues <- rep(c(-4, 0, 4), 200)
    X <-  matrix(Z + mues, ncol = 2)
    aux = ktaucenterscpp::ktaucenters(X, K = 3, nstart = 10)
}




## Bottlenecks?
Rprof()
set.seed(3)
Z <- rnorm(600)
mues <- rep(c(-4, 0, 4), 200)
X <-  matrix(Z + mues, ncol = 2)
ktaucenterscpp::ktaucenters(X, K = 3, nstart = 10)
Rprof(NULL)
summaryRprof()

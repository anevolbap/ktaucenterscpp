##
## R vs Cpp comparison
## 10-FEB-2021
##

library(Rcpp)
library(microbenchmark)

library(tidyverse)

## Data
set.seed(5)
Z <- rnorm(600)
mues <- rep(c(-4, 0, 4), 200)
X <-  matrix(Z + mues, ncol = 2)
X[sample(1:300, 60), ] <- matrix(runif(40, 2 * min(X), 2 * max(X)),
                                ncol = 2, nrow = 60)
set.seed(7)
k1 = ktaucenterscpp::ktaucenters(as_tibble(X), centers = 6, n_runs = 10)
set.seed(7)

k2 = ktaucenters_original(X, centers = 4)

sort(k1$tauPath)
sort(k2$tauPath)

## ---- 
m = microbenchmark::microbenchmark(
                ktaucenterscpp::ktaucenters(X, centers = 5, n_runs = 10),
                ktaucenters::ktaucenters(X, K = 5, nstart = 10), times = 50)

m = bench::mark(check=FALSE, filter_gc = FALSE,
    ktaucenterscpp::ktaucenters(X, centers = 5, n_runs = 10),
    ktaucenters::ktaucenters(X, K = 5, nstart = 10), iterations = 50)

plot(m)


a = ktaucenterscpp::ktaucenters(X, centers = 5, n_runs = 500)

b = ktaucenters::ktaucenters(X, K = 5, nstart = 500)

## Performance improved?
start = Sys.time()
for (iter in 1:10){
    set.seed(iter)
    Z <- rnorm(600)
    mues <- rep(c(-4, 0, 4), 200)
    X <-  matrix(Z + mues, ncol = 2)
    X[sample(1:300, 60), ] <- matrix(runif(40, 2 * min(X), 2 * max(X)),
                                 ncol = 2, nrow = 60)
    aux = ktaucenters::ktaucenters(X, K = 3, nstart = 10)
}
end = Sys.time() - start

ktaucenters_original = function (X, K, centers = NULL, tolmin = 1e-06, NiterMax = 100, 
          nstart = 1, startWithKmeans = TRUE, startWithROBINPD = TRUE, 
          cutoff = 0.999) {
    
    if (!is.matrix(X)) {
        X = as.matrix(X)
    }
    init_centers <- centers
    taumin <- 1e+20
    n <- nrow(X)
    p <- ncol(X)
    aux = list()
    centers0 <- matrix(0, nrow = K, ncol = p)
    start = 1 * (!startWithKmeans)
    nstartEnd = nstart + 1 * (startWithROBINPD)
    for (trial in start:nstartEnd) {
        set.seed(1)
        if (trial == 0) {
            sal0 <- kmeans(X, centers = K, nstart = 20)
            sal0$labels <- sal0$cluster
            for (jota in 1:K) {
                if (sum(sal0$labels == jota) == 1) {
                    centers0[jota, ] <- X[sal0$labels == jota, 
                    ]
                }
                if (sum(sal0$labels == jota) > 1) {
                    centers0[jota, ] <- apply(as.matrix(X[sal0$labels == 
                                                              jota, ]), 2, mean)
                }
            }
        }
        if (trial >= 1) {
            centers0 = X[sample(1:dim(X)[1], K), ]
        }
        if ((trial == 1) & (!is.null(init_centers))) {
            centers0 = init_centers
        }
        if (trial == nstart + 1) {
            retROB <- ktaucenters::ROBINDEN(D = dist(X), data = X, k = K)
            centers0 <- X[retROB$centers, ]
        }
        centers = centers0
        ret_ktau = ktaucenters::ktaucenters_aux(X = X, K = K, centers = centers, 
                                   tolmin = tolmin, NiterMax = NiterMax)
        tauPath = ret_ktau$tauPath
        niter = ret_ktau$niter
        aux = append(aux, tauPath[niter])
        if (tauPath[niter] < taumin) {
            taumin = tauPath[niter]
            best_tauPath = tauPath
            best_ret_ktau = ret_ktau
        }
    }
    newClusters <- best_ret_ktau$cluster
    squaredi <- (best_ret_ktau$di)^2
    robustScale = Mscale(u = sqrt(squaredi), b = 0.5, c = normal_consistency_constants(p))
    outliers = c()
    value <- qchisq(cutoff, df = p)
    for (j in 1:K) {
        indices <- which(newClusters == j)
        booleansubindices <- (squaredi[indices]/(robustScale^2)) > 
            value
        outliersk <- indices[booleansubindices]
        outliers <- c(outliersk, outliers)
    }
    best_ret_ktau$outliers = outliers
    best_ret_ktau$aux = aux
    best_ret_ktau
}


# Distancias
k = 5
p = 2
n = 100
centers = matrix(runif(p * k), ncol = p)
data = matrix(rnorm(n * p), nrow = n)

bench::mark(d1 <- sqrt(apply(t(centers), 2, 
                             function(m){
                                 colSums((t(data) - m)^2)
                             })),
            d2 <- distance_to_centers(data, centers)$distance_matrix,iterations = 20000L)


d1 <- sqrt(apply(t(centers), 2, 
                 function(m){
                     colSums((t(data) - m)^2)
                 }))

d2 <- distance_to_centers(data, centers)

distances_min <- apply(d1, 1, min)

# Membership vector (when ties, first is chosen)
cluster <- max.col(-d1, ties = "first")

d2$membership 

cluster

all(d1 == d2$distance_matrix)

distances_min == d2$min_distance


# Distancias
k = 5
p = 2
n = 100
centers = matrix(runif(p * k), ncol = p)
data = matrix(rnorm(n * p), nrow = n)
bench::mark(d1 <- ktaucenters_aux(data, centers, TOLERANCE_DEFAULT, MAX_ITER_DEFAULT),
            d2 <- ktaucenters_aux_cpp(data, centers, TOLERANCE_DEFAULT, MAX_ITER_DEFAULT), 
            iterations = 100L, check = TRUE)


a = runif(100)
b = cluster
cc = replace(rowsum(a,b),list = 1:2, c(0,0))

a

a/cc[b]




cero = rep(0,K)
cero[rownames(rowsum(cluster, cluster))] = rowsum(cluster, cluster)

K=5
aux = c(cluster, seq(K))
bench::mark(table(aux),
            rowsum(rep(1,aux), aux), check = FALSE)

1L * length(aux)

table(aux)
rowsum(rep(1, length(aux)),aux)


# Distancias
k = 5
p = 2
n = 100
centers = matrix(runif(p * k), ncol = p)
data = matrix(rnorm(n * p), nrow = n)
cluster = sample(1:k , n, replace = TRUE)
weights = runif(n)

mapply(colWeightedMeans, split.data.frame(data, cluster), split(weights, cluster))

rowsum(data * weights, cluster)

dataW = crossprod(diag(weights), data)
Reduce(rbind,
       by(dataW,
          cluster,
          colSums))

set.seed(1)
k1=ktaucenterscpp::ktaucenters(X, centers = 5, n_runs = 2)
set.seed(1)
k2=ktaucenters::ktaucenters(X, K = 5, nstart = 2)

names(k1)
names(k2)

max(abs(k1$tauPath - k2$tauPath))

max(abs(k1$weights - k2$weights))

class(k2$weights)

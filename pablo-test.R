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

## ---- 
m = bench::mark(check=FALSE, filter_gc = FALSE,
    ktaucenterscpp::ktaucenters(X, centers = 5, n_runs = 20),
    ktaucenters::ktaucenters(X, K = 5, nstart = 20), 
    iterations = 50)

set.seed(6)
Z <- rnorm(600)
mues <- rep(c(-4, 0, 4), 200)
X <- matrix(Z + mues, ncol = 2)
X[sample(1:300, 60), ] <- matrix(runif(40, 2 * min(X), 2 * max(X)),
                                 ncol = 2, nrow = 60)
set.seed(6)
ktau_cpp <- ktaucenterscpp::ktaucenters(data=X, centers = 4, n_runs = 1)


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

ktau_aux_orig = function (X, K, centers, tolmin, NiterMax) {
    if (!is.matrix(centers)) {
        centers = as.matrix(centers)
    }
    if (!is.matrix(X)) {
        X = as.matrix(X)
    }
    emptyCluster <- FALSE
    emptyClusterFlag <- FALSE
    n <- nrow(X)
    p <- ncol(X)
    c1 <- constC1(p)
    b1 <- 0.5
    c2 <- constC2(p)
    b2 <- 1
    tauPath <- c()
    niter <- 0
    tol <- tolmin + 1
    repeatedCentersMatrix <- matrix(0, ncol = p, nrow = n)
    distances <- matrix(0, ncol = K, nrow = n)
    while ((niter < NiterMax) & (tol > tolmin)) {
        for (h in 1:K) {
            distances[, h] <- sqrt(apply(sweep(X, 2, centers[h, 
            ])^2, 1, sum))
        }
        distances_min <- apply(distances, 1, min)
        cluster <- apply(distances, 1, function(x) which(x == 
                                                             min(x))[1])
        ms <- Mscale(u = distances_min, b = b1, c = c1)
        dnor <- distances_min/ms
        tau <- ms * sqrt(mean(rhoOpt(dnor, cc = c2)))/sqrt(b2)
        tauPath <- c(tauPath, tau)
        Du <- mean(psiOpt(dnor, cc = c1) * dnor)
        Cu <- mean(2 * rhoOpt(dnor, cc = c2) - psiOpt(dnor, cc = c2) * 
                       dnor)
        Wni <- (Cu * psiOpt(dnor, cc = c1) + Du * psiOpt(dnor, 
                                                         cc = c2))/dnor
        Wni1 = Wni
        if (sum(distances_min == 0) > 0) {
            Wni[distances_min == 0] <- (Du * derpsiOpt(0, cc = c2) + 
                                            Cu * derpsiOpt(0, cc = c1))
        }
        Wni2 = Wni
        weights <- 0 * Wni
        for (jota in 1:K) {
            if ((sum(Wni[cluster == jota])) != 0) {
                weights[cluster == jota] = Wni[cluster == jota]/sum(Wni[cluster == 
                                                                            jota])
            }
            if ((sum(Wni[cluster == jota])) == 0) {
                mmm = length(cluster == jota)
                if (sum(weights[cluster == jota] == 0) == mmm) 
                    weights[cluster == jota] = 1/mmm
            }
        }
        XW <- 0 * X
        XW = sweep(X, 1, weights, FUN = "*")
        oldcenters <- centers
        auxx = rep(0, K)
        for (jota in 1:K) {
            auxx[jota] <- sum(cluster == jota)
            if (auxx[jota] > 0) {
                centers[jota, ] = apply(as.matrix(XW[cluster == 
                                                         jota, ]), 2, sum)
            }
        }
        if (sum(auxx > 0) != K) {
            furtherIndices = order(distances_min, decreasing = TRUE)[1:sum(auxx == 
                                                                               0)]
            centers[auxx == 0, ] = X[furtherIndices, ]
            cluster[furtherIndices] = which(auxx == 0)
            emptyClusterFlag = TRUE
        }
        tol = sqrt(sum((oldcenters - centers)^2))
        niter = niter + 1
    }
    ret = list(tauPath = tauPath, niter = niter, centers = centers, 
               cluster = cluster, emptyCluster = emptyCluster, tol = tol, 
               weights = weights, di = distances_min, 
               Wni = Wni, Wni1 = Wni1, Wni2 = Wni2,  
               emptyClusterFlag = emptyClusterFlag)
    return(ret)
}

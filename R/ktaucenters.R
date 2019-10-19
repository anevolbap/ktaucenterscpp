#' ktaucenters
#'
#' Robust Clustering algorithm based on centers, a robust and efficient version of KMeans.
#' @param X numeric matrix  of size n x p.
  #' @param K the number of cluster.
#' @param centers a matrix of size K x p containing the K initial centers,
#'  one at each matrix-row. If centers is NULL a random set of (distinct) rows in  \code{X}
#'  are chosen as the initial centres.
#' @param tolmin a tolerance parameter used for the algorithm stopping rule
#' @param NiterMax a maximun number of iterations used for the algorithm stopping rule
#' @param nstart the number of trials that the base algorithm ktaucenters_aux is run.
#' If it is greater than 1 and center is not set as NULL, a random set of (distinct) rows
#' in \code{X} will be chosen as the initial centres.
#' @param startWithKmeans  TRUE if kmean centers values is included as starting point.
#' @param startWithROBINPD TRUE if ROBINDEN estimator is included as starting point
#' @param cutoff optional argument for outliers detection - quantiles of chi-square to be used as a threshold for outliers detection, defaults to 0.999


#' @return A list including the estimated K centers and labels for the observations
##' \itemize{
##'  \item{\code{centers}}{:   matrix of size K x p, with the estimated K centers.}
##'  \item{\code{cluster}}{: array of size n x 1  integers labels between 1 and K.}
##'  \item{\code{tauPath}}{: sequence of tau scale values at each iterations.}
##'  \item{\code{Wni}}{: numeric array of size n x 1 indicating the weights associated to each observation.}
##'  \item{\code{emptyClusterFlag}}{: a boolean value. True means that in some iteration there were clusters totally empty}
##'  \item{\code{niter}}{: number of iterations untill convergence is achived or maximun number of iteration is reached}
##'  \item{\code{di}}{: distance of each observation to its assigned cluster-center}
##'  \item{\code{outliers}}{: indices observation that can be considered as outliers}
##' }

#' @examples
#' # Generate Sintetic data (three cluster well separated)
#' Z <- rnorm(600);
#' mues <- rep(c(-3, 0, 3), 200)
#' X <-  matrix(Z + mues, ncol=2)
#'
#' # Generate 60 sintetic outliers (contamination level 20%)
#' X[sample(1:300,60), ] <- matrix(runif( 40, 3 * min(X), 3 * max(X) ),
#'                                 ncol = 2, nrow = 60)
#'
#' ### Applying the algortihm ####
#'sal <- ktaucenters(
#'      X, K=3, centers=X[sample(1:300,3), ],
#'      tolmin=1e-3, NiterMax=100)
#'
#' #### plotting  the clusters####
#'
#' oldpar=par(mfrow = c(1,2))
#' 
#' plot(X,type = "n", main = "ktaucenters (Robust) \n outliers: solid black dots")
#' points(X[sal$cluster==1,],col=2);
#' points(X[sal$cluster==2,],col=3);
#' points(X[sal$cluster==3,],col=4)
#' points(X[sal$outliers,1], X[sal$outliers,2], pch=19)
#'
#' ### Applying a classical (non Robust) algortihm ###
#' sal <- kmeans(X, centers=3,nstart=100)
#'
#' ### plotting the clusters ###
#' plot(X, type ="n", main = "kmeans (Classical)")
#' points(X[sal$cluster==1,],col=2);
#' points(X[sal$cluster==2,],col=3);
#' points(X[sal$cluster==3,],col=4)
#'
#' par(oldpar)
#' @references Gonzalez, J. D., Yohai, V. J., & Zamar, R. H. (2019). 
#' Robust Clustering Using Tau-Scales. arXiv preprint arXiv:1906.08198. 
#'
#'
#'
#'
#' @importFrom stats kmeans dist qchisq
#' @export
ktaucenters=function(X,K,centers=NULL,tolmin=1e-06,NiterMax=100,nstart=1,startWithKmeans=TRUE,startWithROBINPD=TRUE,cutoff=0.999){
  if (!is.matrix(X)){
    X=as.matrix(X);
  }
  init_centers <- centers
  taumin <- 1e+20
  n <- nrow(X);
  p <- ncol(X);
  centers0 <- matrix(0, nrow=K, ncol=p)
  start=1*(!startWithKmeans); # the start value is one or zero.

  nstartEnd=nstart+ 1*(startWithROBINPD);

  for (trial in start:nstartEnd){
    # if startWithKmeans its true, start=0, then trial take the zero value.
    if (trial==0){
      sal0 <- kmeans(X, centers=K, nstart = 20)
      sal0$labels <- sal0$cluster
      for (jota in 1: K){
        # when there is a single observation, it is not possible to use apply function
        if (sum(sal0$labels==jota)==1){
          centers0[jota, ] <- X[sal0$labels==jota, ]
        }

        if (sum(sal0$labels==jota)>1){
          # as.matrix below is necessary because if p=1 function "apply" will not work.
          centers0[jota, ] <- apply(as.matrix(X[sal0$labels==jota, ]), 2, mean)
        }
      }
    }

    if (trial>=1){
      centers0=X[sample(1:dim(X)[1],K), ]
    }

    if ((trial==1) & (!is.null(init_centers))){
      centers0=init_centers;
    }
    if (trial==nstart+ 1){
      retROB <- ROBINDEN(D = dist(X), data = X, k = K);
      centers0 <- X[retROB$centers, ]
    }

    centers=centers0;
    ret_ktau=ktaucenters_aux(X = X,K = K,
                             centers = centers,tolmin = tolmin,
                             NiterMax = NiterMax)
    tauPath=ret_ktau$tauPath;
    niter=ret_ktau$niter

    if (tauPath[niter]<taumin) {
      # si la escala es menor que las otras actualizo
      taumin=tauPath[niter];
      best_tauPath=tauPath
      best_ret_ktau=ret_ktau
    }
  }

  newClusters<-best_ret_ktau$cluster;
  squaredi <-(best_ret_ktau$di)^2;
  robustScale=Mscale(u=sqrt(squaredi),b=0.5,c = normal_consistency_constants(p));
  outliers=c()
  value <- qchisq(cutoff, df = p)
  for (j in 1:K){
    indices <- which(newClusters == j)
    booleansubindices <-(squaredi[indices]/(robustScale^2)) > value
    outliersk <- indices[booleansubindices]
    outliers <- c(outliersk, outliers)
  }
  best_ret_ktau$outliers=outliers
  best_ret_ktau
}

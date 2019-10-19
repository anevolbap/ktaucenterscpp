  #' ktaucenters_aux
  #'
  #' Robust Clustering algorithm based on centers, a robust and efficient version of KMeans.
  #' @param X A matrix  of size n x p.
  #' @param K The number of clusters.
  #' @param centers  matrix of size K x p containing the K initial centers, one at each matrix-row.
  #' @param tolmin  tolerance parameter used for the algorithm stopping rule
  #' @param NiterMax a maximun number of iterations used for the algorithm stopping rule
  #' @return A list including the estimated K centers and labels for the observations
  #' \itemize{
  #'  \item{\code{centers}}{:   matrix of size K x p, with the estimated K centers.}
  #'  \item{\code{cluster}}{: array of size n x 1  integers labels between 1 and K.}
  #'  \item{\code{tauPath}}{: sequence of tau scale values at each iterations.}
  #'  \item{\code{Wni}}{: numeric array of size n x 1 indicating the weights associated to each observation.}
  #'  \item{\code{emptyClusterFlag}}{: a boolean value. True means that in some iteration there were clusters totally empty}
  #'  \item{\code{niter}}{: number of iterations untill convergence is achived or maximun number of iteration is reached}
  #'  \item{\code{di}}{distance of each observation to its assigned cluster-center}
  #'  }
  #' @examples
  #'
  #' # Generate Sintetic data (three cluster well separated)
  #' Z=rnorm(600);
  #' mues=rep(c(0,10,20),200)
  #' X= matrix(Z+mues,ncol=2)
  #'
  #' # Applying the algortihm
#' sal = ktaucenters_aux(
#' X, K=3, centers=X[sample(1:300,3), ],
#' tolmin=1e-3, NiterMax=100)
#'
#' #plot the results
#' plot(X,type="n")
#' points(X[sal$cluster==1,],col=1);
#' points(X[sal$cluster==2,],col=2);
#' points(X[sal$cluster==3,],col=3);
#'
#' @note
#' Some times, if the initial centers are wrong, the algorithm converges to a non-optimal (local) solution.
#' To avoid that, the algorithm must be run several times. This task is carried out by \code{\link{ktaucenters}}
#' @seealso \code{\link{ktaucenters}}
#' @references Gonzalez, J. D., Yohai, V. J., & Zamar, R. H. (2019). 
#' Robust Clustering Using Tau-Scales. arXiv preprint arXiv:1906.08198. 
#' @export
ktaucenters_aux=function(X,K,centers,tolmin,NiterMax){
  if (!is.matrix(centers)){
    centers=as.matrix(centers);
  }


  if (!is.matrix(X)){
    X=as.matrix(X);
  }
  emptyCluster <- FALSE
  emptyClusterFlag <- FALSE
  n <- nrow(X);
  p <- ncol(X);
  c1 <- constC1(p)
  b1 <-  .5
  c2 <- constC2(p)
  b2 <- 1
  tauPath <- c();
  niter <- 0;
  tol <- tolmin+1
  repeatedCentersMatrix <- matrix(0,ncol=p,nrow=n)
  distances <- matrix(0,ncol=K,nrow=n)
  while((niter<NiterMax) & (tol>tolmin)){
    for (h in 1:K){

#      for (iw in 1:n){
#        repeatedCentersMatrix[iw, ] <- centers[h,]
#      }
#     a cada fila de x le resto el centerside y despues calculo la distancia
#     resultado: una matriz donde en cada fila tiene di1 di2  diK, donde
#     dij es la distancia entre xi y el centerside muj
#     OLD distances[, h] <- sqrt((X - repeatedCentersMatrix)^2  %*% rep(1,p))
      distances[, h] <- sqrt(apply(sweep(X,2,centers[h, ])^2,1,sum))
    }
    #  Taking the minimun distsance from each observation to its center.
    distances_min <- apply(distances, 1, min)
    # cluster is an array of nx1, that contains integers from 1 to K
    # cluster[l] = j means  that observation x_l is assigned to cluster j.
    cluster <- apply(distances, 1, function(x) which(x==min(x))[1])
    ms <- Mscale(u=distances_min, b=b1, c=c1)
    dnor <- distances_min/ms; # normalized distance
    tau <- ms*sqrt(mean(rhoOpt(dnor,cc=c2)))/sqrt(b2)
    # vector that save tau scale values of the distances
    tauPath <- c(tauPath,tau)

    ###############################################################################
    # define weights and constants in each iteration
    Du <-  mean(psiOpt(dnor, cc=c1)*dnor);
    Cu <-  mean(2*rhoOpt(dnor, cc=c2)-psiOpt(dnor, cc=c2)*dnor)
    Wni <-  (Cu*psiOpt(dnor, cc=c1) + Du*psiOpt(dnor, cc=c2))/dnor
    ###############################################################################


    # Atenttion: when di=0. psi_1(dnor)=0 y psi_2(dnor)=0.
    # Then the weight w is undefined due to  dividing by zero.
    # Given that psi_1(0)=0, Wni can be obtained throug the derivative
    # of psi_1 in this case. that is that the following lines do.
    if(sum(distances_min==0)>0){
      Wni[distances_min==0] <- (Du * derpsiOpt(0,cc=c2) + Cu * derpsiOpt(0,cc=c1))
    }

    weights <- 0*Wni;

    for (jota in 1:K){
      if ( (sum(Wni[cluster==jota])) !=0 ){
        weights[cluster==jota]=Wni[cluster==jota]/sum(Wni[cluster==jota])
      }
      ## check this piece of code
      if ( (sum(Wni[cluster==jota]))==0 ){
        # esto significa que el algoritmo converjio.
        # entonces no  actualizo las X's.
        # pero si llega a pasar  que son cero las actualizo:
        mmm=length(cluster==jota)
        if (sum(weights[cluster==jota]==0)==mmm)
          weights[cluster==jota]=1/mmm;
      }
    }


    XW <- 0*X

    # each X row is multiplied by their corresponding weight
    # OLD OLD OLD OLD OLD OLDOLD OLD OLDOLD OLD OLD
    #for (fila in 1:n){
    #  XW[fila, ] <- X[fila, ] * weights[fila]
    #}

    # New  New New New New New New New New  New New New New New New New
    XW=sweep(X,1,weights,FUN='*')

    # New centers are named centers.
    oldcenters <- centers;
    ## sometimes a cluster has no observations, the followin code
    ## deals that situation ...
    auxx=rep(0, K)

    for (jota in 1:K){
      auxx[jota] <- sum(cluster==jota)
      if (auxx[jota]>0){
        # The "as.matrix" in the line below
        # is nedeed when the matrix has a single observation.
        # In that case R transforms it
        # into a vector and the assignment is not possible
        centers[jota,]=apply(as.matrix(XW[cluster==jota,]),2,sum)
      }
    }

    if (sum(auxx>0)!=K){
      # if not all clusters are filled,
      # ther are replaced for the furthest centers
      # this is speccially important when the number of clusters
      # K is high.
      furtherIndices= order(distances_min,decreasing=TRUE)[1:sum(auxx==0)];
      centers[auxx==0,]=X[furtherIndices, ]
      cluster[furtherIndices]=which(auxx==0)
      emptyClusterFlag=TRUE;
    }

    # condition sum(auxx>0)==K means that  all clusters are filled.
    tol=sqrt(sum((oldcenters-centers)^2));

    niter=niter+1
  }
  ret=list(tauPath=tauPath, niter=niter, centers=centers,cluster=cluster,emptyCluster=emptyCluster, tol=tol,weights=weights,di=distances_min,Wni=Wni,emptyClusterFlag=emptyClusterFlag)
  return(ret)
}







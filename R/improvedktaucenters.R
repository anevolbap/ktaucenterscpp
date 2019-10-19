#' improvedktaucenters
#'
#' Robust Clustering algorithm for non-spherical data. This function estimate clusters taking into account that clusters may have
#' different size, volume or orientation.
#' @param X numeric matrix  of size n x p.
#' @param K the number of cluster.
#' @param cutoff optional argument for getOutliers - quantiles of chi-square to be used as a threshold for outliers detection, defaults to 0.999
#' @param nstart optional the number of trials that the base algorithm ktaucenters_aux is run
#' at the first stage. #' If it is greater than 1 and center is not set as NULL,
#' a random set of (distinct) rows in x is chosen as the initial centres for each try.
#' @param INITcenters numeric matrix  of size K x p indicating the initial centers for that clusters and robust covarianze matrices will be 
#' computed, if it is set as NULL the algorithm will compute @param INITcenters from ktaucenters routine. Set to NULL by default. 
#' @return A list including the estimated K centers and clusters labels for the observations

## @details text describing parameter inputs in more detail.
#' \itemize{
#'  \item{\code{centers}}{:   matrix of size K x p, with the estimated K centers.}
#'  \item{\code{cluster}}{: array of size n x 1  integers labels between 1 and K.}
#'  \item{\code{tauPath}}{: sequence of tau scale values at each iterations.}
#'  \item{\code{Wni}}{: numeric array of size n x 1 indicating the weights associated to each observation.}
#'  \item{\code{emptyClusterFlag}}{: a boolean value. True means that in some iteration there were clusters totally empty.}
#'  \item{\code{niter}}{: number of iterations untill convergence is achived or maximun number of iteration is reached.}
#'  \item{\code{sigmas}}{: a list containing the k covariance matrices found by the procedure at its second step.}
#'  \item{\code{outliers}}{: indices observation that can be considered as outliers.}
#' }
#' @importFrom GSE GSE   getOutliers  getLocation  getScatter
#' @importFrom MASS ginv
#' @importFrom stats   mahalanobis  qchisq
## #
#' @examples
#'
#' # Generate Sintetic data (three normal cluster in two dimension)
#' # clusters have different shapes and orentation.
#' # The data is contaminated uniformly (level 20%).
#'
#' ################################################
#' #### Start data generating process ############
#' ##############################################
#'
#' # generates base clusters
#' set.seed(1)
#' Z1 <- c(rnorm(100,0),rnorm(100,0),rnorm(100,0))
#' Z2 <- rnorm(300);
#' X <-  matrix(0, ncol=2,nrow=300);
#' X[,1]=Z1;X[,2]=Z2
#' true.cluster= c(rep(1,100),rep(2,100),rep(3,100))
#'
#' # rotate, expand and translate base clusters
#' theta=pi/3;
#' aux1=matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),nrow=2)
#' aux2=sqrt(4)*diag(c(1,1/4))
#' B=aux1%*%aux2%*%t(aux1)
#' X[true.cluster==3,]=X[true.cluster==3,]%*%aux2%*%aux1 + matrix(c(5,2),byrow = TRUE,nrow=100,ncol=2)
#' X[true.cluster==2,2]=X[true.cluster==2,2]*5
#' X[true.cluster==1,2]=X[true.cluster==1,2]*0.1
#' X[true.cluster==1, ]=X[true.cluster==1,]+ matrix(c(-5,-1),byrow = TRUE,nrow=100,ncol=2)

#' ### Generate 60 sintetic outliers (contamination level 20%)
#'
#' outliers=sample(1:300,60)
#' X[outliers, ] <- matrix(runif( 40, 2 * min(X), 2 * max(X) ),
#'                                 ncol = 2, nrow = 60)
#' ###############################################
#' #### END data generating process ############
#' #############################################
#'
#' #############################################
#' ### Applying the algortihm ##################
#' #############################################
#' ret=improvedktaucenters(X,K=3,cutoff=0.999)
#'
#' #############################################
#' ### plotting results ########################
#' #############################################
#' oldpar=par(mfrow=c(2,1))
#' #' plot(X,main="actual clusters")
#' for (j in 1:3){
#'  points(X[true.cluster==j,],pch=19, col=j+1)
#' }
#' points(X[outliers,],pch=19,col=1)

#' plot(X,main="clusters estimation")
#' for (j in 1:3){
#'  points(X[ret$cluster==j,],pch=19, col=j+1)
#' }
#' points(X[ret$outliers,],pch=19)
#' 
#' par(oldpar)
#' @references Gonzalez, J. D., Yohai, V. J., & Zamar, R. H. (2019). 
#' Robust Clustering Using Tau-Scales. arXiv preprint arXiv:1906.08198. 
#' @export

improvedktaucenters=function(X,K,cutoff=0.999,nstart =5,INITcenters=NULL){

  n=dim(X)[1]
  p=dim(X)[2]
  sigmas <- vector(mode="list", length=K);
  centers=INITcenters
  #### FIRST STEP: determine the best centers  #####

  # if centers are not given, we calculate with ktaucenters routine
  if (is.null(centers)){
    ret_ktau <- ktaucenters(X, K, nstart = nstart)
    centers <- ret_ktau$centers
    newClusters<-ret_ktau$cluster
  }
  
  #### if centers are given, we  just 
  #### calculate the clusters from these centers. 
  sphericalDistanceMatrix<- matrix(0, ncol = K, nrow = n)
    if (!is.null(centers)){
    for (j in 1:K){
    sphericalDistanceMatrix[,j] <- mahalanobis(x=X,
                                           center=centers[j, ]
                                           ,cov=diag(p))
      
    }  
    newClusters <- apply(sphericalDistanceMatrix, 1,
                             function(x) which( x == min(x))[1])
  }
  
  
  newcentersaux <- 0*centers
  for (j in 1:K){
    sigmas[[j]] <- diag(p)
  }

#=========== and SECOND STEP ========
iter=0
newcenters <- centers;
while(iter<1){
    for (j in 1:K){
      Xcluster <- X[newClusters==j,];
      nk <- dim(Xcluster)[1];
      
      # dim works properly when there are two or more observations, for that reason: 
      if (sum(newClusters==j)<2){
        nk=sum(newClusters==j)
        }
      
      # if there is enough observations, we proceed to compute 
      # robust covarianze matrix 
      if(nk>(3*p)){
        #  Xcentered <- sweep(Xcluster, 2, centers[j,], FUN="-")
          sal1 <- GSE(Xcluster, tol=1e-4, maxiter=150)
          newcentersaux[j,] <- getLocation(sal1)
          sigmas[[j]] <- getScatter(sal1)
      }
      # if there is not enough observations, we don't compute 
      # robust covarianze matrix 
      if(nk<=(3*p)){
        newcentersaux[j,]= newcenters[j, ]
      }
      
    }

    mahalanobisMatrix <- matrix(0, ncol = K, nrow = n)
    for (j in 1:K){
      mahalanobisMatrix[,j] <- mahalanobis(x=X,
                                        center=newcentersaux[j, ]
                                        ,cov=ginv(sigmas[[j]]),
                                        inverted = TRUE)

    }
    newClusters_Mah <- apply(mahalanobisMatrix, 1,
                          function(x) which( x == min(x))[1])
    value <- qchisq(cutoff, df = p)

    #IMPROVED oultiers determination. 
    outliers=c()
    for (j in 1:K){
      indices <- which(newClusters_Mah == j)
      booleansubindices <- mahalanobisMatrix[indices, j] > value
      outliersk <- indices[booleansubindices]
      outliers <- c(outliersk, outliers)
    }

  # In the future: implement this lines should be worked better than current for unbalanced-cluster cases.
  # qsMatrix=matrix(0,ncol=K,nrow=n)
  #  logdet=rep(0,K);

  #  for (j in 1:K){
  #    logalpha=sum(newClusters_Mah==j)/(n-length(outliers))
  #    logdet[j]=log(det(sigmas[[j]]))
  #    qsMatrix[,j]= logalpha[j] - .5*logdet[j] -.5*mahalanobisMatrix[,j]
  #  }
  #  newClusters_qs=apply(qsMatrix,1,function(x) which(x==min(x))[1])

    newClusters <- newClusters_Mah;
    newcenters <- newcentersaux;
    iter <- iter + 1
  }

  list(cluster=newClusters, centers=newcenters, outliers=outliers, sigmas=sigmas);
}

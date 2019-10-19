#' ROBINDEN (ROBust INitialization based on inverse DENsity estimator)
#'
#' @description
#' \code{ROBINDEN} searches for k initial cluster seeds for k-means-based clustering methods.
#'
#' @param D A distance matrix calculated on \code{data}.
#' @param data A data matrix with n observations and p variables.
#' @param k The number of cluster centers to find.
#' @param mp The number of the nearest neighbors to find dense regions by LOF, the default is 10.
#'
#'
#' @return
#' \item{centers}{A numeric vector of \code{k} initial cluster centers corresponding to the k indices of observations.}
#' \item{idpoints}{A real vector containing the inverse density values of each point (observation).}
#'
#' @export
#'
#' @details
#' The centers are the observations located in the most dense region
#' and far away from each other at the same time.
#' In order to find the observations in the highly dense region,
#' ROBINPOINTDEN uses point density estimation
#' (instead of Local Outlier Factor, Breunig et al (2000)), see more details.
#'
#' @examples
# Generate Sintetic data (7 cluster well separated)
#' K=5;
#' nk=100
#' Z <- rnorm(2 * K * nk);
#' centers_aux <- -floor(K/2):floor(K/2)
#' mues <- rep(5*centers_aux,2*nk*K )
#' X <-  matrix(Z + mues, ncol=2)

#' # Generate sintetic outliers (contamination level 20%)
#' X[sample(1:(nk * K),(nk * K) * 0.2), ] <-matrix(runif((nk * K) * 0.2 * 2,
#'                                           3 * min(X), 3 * max(X)),
#'                                           ncol = 2, nrow = (nk * K) * 0.2)
#' res <- ROBINDEN(D =dist(X), data=X, k = K);
#' # plot the Initial centers found
#' plot(X)
#' points(X[res$centers,],pch=19,col=4,cex=2)

#' @note this is a slightly modified version of ROBIN algorithm implementation done by
#' Sarka Brodinova <sarka.brodinova@tuwien.ac.at>.
#' @author Juan Domingo Gonzalez <juanrst@hotmail.com>
#'
#' @references Hasan AM, et al. Robust partitional clustering by
#' outlier and density insensitive seeding. Pattern Recognition Letters, 30(11), 994-1002, 2009.
#'
#' @seealso \code{\link[dbscan]{lof}}
#'
#' @importFrom dbscan lof
#' @importFrom dplyr %>%
#' @importFrom dbscan kNN
#
ROBINDEN <- function(D,data,k,mp=10){
  # THIS FUNCTION DO THE SAME THAN ROBIN BUT TAKing INTO ACCOUNT THE density (Instead of the
  # inverse of the average relative local density (known as LOF))

  # calculates the inverse density point.
  idp <-1/denpoints(D,k=mp)
  n <- nrow(data)
  
  # Observation: 
  # Outliers have a high idp value.
  # In umbalanced cases and when K increases, all the observations from a group might be above the critRobin, 
  # So we need to increase the critRobin in order to avoid two initials centers from the same group. 
  vvv=max(0.5,0.96*(1-1.5/k));
  critRobin=sort(idp,decreasing = FALSE)[floor(vvv*n)]; # vvv=0.5 implies critRobin= median(idp). 
  
  # ORIGINAL r <-sample(n,1)
  # MODIFICATION: Start with a point whose its density is maximun
  r <-which(idp==min(idp))[1]
  # END MODIFICATION.
  id.means <- numeric(k)
  m <-1
  D <- as.matrix(D)
  while(m<=k){
    # find the observations that are far away from each other - the potential cluster seeds
    if(m<=2){
      sort.points <- sort.int(D[r,],decreasing = TRUE) %>%names()%>%as.numeric()
    }else{
      sort.points <- sort.int(apply(D[id.means[1:m],],2,min),decreasing=TRUE) %>%names() %>% as.numeric()
    }
    idp.sort.points <- idp[sort.points]
    id <- which(c(idp.sort.points<=critRobin))[1]
    
    if (all((idp.sort.points>critRobin))){
    # Juan D. Gonzalez  
    #some times id is empty because all idp.sort.points are greater than 
    # the critRobin Value, then I take the nearest point to critRobin
    id=order((idp.sort.points-critRobin),decreasing=FALSE)[1]
    }
    
    r <- sort.points[id]
    id.means[m] <- r

    m <- m+1
  }
  return(list(centers=id.means,idpoints=idp))
}
############################
#########################


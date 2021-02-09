#' Robust initialization based on inverse density estimator
#'
#' @description \code{robust_init_density} searches for k initial cluster seeds
#'     for k-means-based clustering methods.
#'
#' @param dist_matrix A distance matrix calculated on \code{data}.
#' @param data A data matrix with n rows and p columns.
#' @param k The number of cluster centers to find.
#' @param mp The number of the nearest neighbors to find dense regions
#'     by LOF, the default is 10.
#'
#' @return
#' \item{centers}{A numeric vector of \code{k} initial cluster
#'     centers corresponding to the k indices of observations.}
#' \item{idpoints}{A real vector containing the inverse density
#'     values of each point (observation).}
#'
#' @export
#'
#' @details 
#'   This function do the same as ROBIN but taking into account the 
#'   density (instead of the inverse of the average relative local
#'   density known as LOF)
#'   
#'   The centers are the observations located in the most dense
#'     region and far away from each other at the same time.  In order
#'     to find the observations in the highly dense region,
#'     ROBINPOINTDEN uses point density estimation (instead of Local
#'     Outlier Factor, Breunig et al. (2000)), see more details.
#'     
#'     Observation: Outliers have a high 'idp' value. In imbalanced cases
#'     and when K increases, all the observations from a group might be
#'     above the critRobin, So we need to increase the critRobin in order
#'     to avoid two initials centers from the same group.
#'     modification: start with a point whose density is maximum
#'     
#' @examples
# Generate Synthetic data (7 cluster well separated)
#' K=5;
#' nk=100
#' Z <- rnorm(2 * K * nk);
#' centers_aux <- -floor(K/2):floor(K/2)
#' mues <- rep(5*centers_aux,2*nk*K )
#' X <-  matrix(Z + mues, ncol=2)

#' # Generate synthetic outliers (contamination level 20%)
#' X[sample(1:(nk * K),(nk * K) * 0.2), ] <-matrix(runif((nk * K) * 0.2 * 2,
#'                                           3 * min(X), 3 * max(X)),
#'                                           ncol = 2, nrow = (nk * K) * 0.2)
#' res <- robust_init_density(dist_matrix =dist(X), data=X, k = K);
#' # plot the Initial centers found
#' plot(X)
#' points(X[res$centers,],pch=19,col=4,cex=2)

#' @note this is a slightly modified version of ROBIN algorithm
#'     implementation done by Sarka Brodinova
#'     <sarka.brodinova@tuwien.ac.at>.
#' @author Juan Domingo Gonzalez <juanrst@hotmail.com>
#'
#' @references Hasan AM, et al. Robust partitional clustering by
#'     outlier and density insensitive seeding. Pattern Recognition
#'     Letters, 30(11), 994-1002, 2009.
#'
#' @seealso \code{\link[dbscan]{lof}}
#'
#' @importFrom dbscan lof
#' @importFrom dbscan kNN
#
robust_init_density <- function(dist_matrix, data, n_clusters, 
                                mp = 10, method = "density") {
  
  id_means <- numeric(n_clusters)
  iter <- 1  
  n <- nrow(data)
  
  # Compute the inverse density points.
  idp <- density_points(dist_matrix, k = mp)
  position <- floor(max(0.5, 0.96 * (1 - 1.5 / n_clusters)) * n)
  crit_robin <- sort(idp)[position]
  
  if (method == "robin") r <- sample(n, 1)
  if (method == "density") r <- which.min(idp)
  
  dist_matrix = as.matrix(dist_matrix) # FIXME: check types 
  while(iter <= n_clusters){
    if (iter <= 2) {
      sorted_points <- order(dist_matrix[r, ], decreasing = TRUE)
    } else {
      sorted_points <- order(apply(dist_matrix[id_means[1:iter], ], 2, min), #FIXME: try Rfast package
                             decreasing = TRUE)
    }
    idp_sorted_points <- idp[sorted_points]
    
    id <- ifelse(any(idp_sorted_points <= crit_robin),
                 which(idp_sorted_points <= crit_robin)[1], # FIXME: efficient "first true"
                 which.min(idp_sorted_points - crit_robin, decreasing))
    
    r <- sorted_points[id]
    id_means[iter] <- r
    
    iter <- iter + 1
  }
  
  return(list(centers = id_means,
              id_points = idp))
}
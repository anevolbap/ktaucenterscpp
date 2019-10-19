#' denpoints
#' @description
#' \code{denpoints} Estimates the densities values of a sample.
#'
#' @param x A distance matrix calculated on \code{data} or a matrix
#' @param k The number of nearest neighbors to calculate local point density
#'
#' @return
#' \item{dpoints}{A real vector containing the density values for each point}
#'
#' @export
#'
#' @details
#'
#' For a fixed  \code{y}, density of \code{y} is  defined as the sum of  \code{distance(y,z)}
#' on all \code{z} that are the k-nearest neighbors of \code{y}
#'
#' @examples
#' #generate normal data in dimension 2
#' X=matrix(rnorm(1000),ncol=2)
#' a<- denpoints(x=X,k=4)
#'
#'
#' #### ten most isolated points
#' most_isolated=order(a)[1:10]
#'
#' ### plotting results: (most isolated points should be shown in green)
#' plot(X)
#' points(X[ most_isolated, ], pch=19,col=3)
#'
#' @author Juan Domingo Gonzalez <juanrst@hotmail.com>
#'
#' @references Hasan AM, et al. Robust partitional clustering by
#' outlier and density insensitive seeding. Pattern Recognition Letters, 30(11), 994-1002, 2009.
#' @importFrom dbscan kNN
#' @importFrom methods is

denpoints=function (x, k = 4)
{
  if (is(x, "dist"))
    n <- attr(x, "Size")
  else n <- nrow(x)
  if (is.null(n))
    stop("x needs to be a matrix or a dist object!")
  if (k < 1 || k >= n)
    stop("k has to be larger than 1 and smaller than the number of points")
  d <- kNN(x, k)
  lrd <- numeric(n)

  for (i in 1:n) lrd[i] <- 1/(sum(apply(cbind(d$dist[d$id[i,
                                                          ], k], d$dist[i, ]), 1, max))/k)
  pointDen <- lrd
  pointDen[is.nan(pointDen)] <- 1
  dpoints=pointDen
  dpoints
}



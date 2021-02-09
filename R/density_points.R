#' denpoints
#'
#' @description \code{density_points} Estimates the densities values of a
#'     sample.
#'
#' @param x A distance matrix calculated on \code{data} or a matrix
#' @param k The number of nearest neighbors to calculate local point
#'     density
#'
#' @return \item{dpoints}{A real vector containing the density values
#'     for each point}
#'
#' @export
#'
#' @details
#'
#' For a fixed \code{y}, density of \code{y} is defined as the sum of
#' \code{distance(y,z)} on all \code{z} that are the k-nearest
#' neighbors of \code{y}
#'
#' @examples
#' ## generate normal data in dimension 2
#' X <- matrix(rnorm(1000), ncol = 2)
#' a <- denpoints(x = X, k = 4)
#'
#'
#' ## ten most isolated points
#' most_isolated = order(a)[1:10]
#'
#' ## plotting results: (most isolated points should be shown in green)
#' plot(X)
#' points(X[ most_isolated, ], pch = 19, col = 3)
#'
#' @author Juan Domingo Gonzalez <juanrst@hotmail.com>
#'
#' @references Hasan AM, et al. Robust partitional clustering by
#'     outlier and density insensitive seeding. Pattern Recognition
#'     Letters, 30(11), 994-1002, 2009.
#' @importFrom dbscan kNN
#' @importFrom methods is
density_points <- function(x, k = 4) {
    
    d <- kNN(x, k)
    pointDen <- rowMeans(pmax(matrix(d$dist[t(d$id), k], ncol=k, byrow = TRUE),
                              d$dist))
    pointDen[is.nan(pointDen)] <- 1
    return(pointDen)
}
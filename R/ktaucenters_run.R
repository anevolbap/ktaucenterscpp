#' ktaucenters_run
#'
#' Robust Clustering algorithm based on centers, a robust and
#' efficient version of K-Means.
#' @param data A matrix of size n x p.
#' @param K The number of clusters.
#' @param centers matrix of size K x p containing the K initial
#'     centers, one at each matrix-row.
#' @param tolerance tolerance parameter used for the algorithm stopping
#'     rule
#' @param max_iter a maximum number of iterations used for the
#'     algorithm stopping rule
#' @return A list including the estimated K centers and labels for the
#'     observations
#'
#' \itemize{
#' \item{\code{centers}}{: matrix of size K
#'     x p, with the estimated K centers.}
#' \item{\code{cluster}}{:
#'     array of size n x 1 integers labels between 1 and K.}
#' \item{\code{tauPath}}{: sequence of tau scale values at each
#'     iterations.}
#' \item{\code{Wni}}{: numeric array of size n x 1
#'     indicating the weights associated to each observation.}
#' \item{\code{empty_cluster_flag}}{: a boolean value. True means
#'     that in some iteration there were clusters totally empty}
#' \item{\code{niter}}{: number of iterations until convergence
#'     is achived or maximum number of iteration is reached}
#' \item{\code{di}}{distance of each observation to its assigned
#'     cluster-center} }
#' @examples
#'
#' # Generate Synthetic data (three cluster well separated)
#' Z=rnorm(600);
#' mues=rep(c(0,10,20),200)
#' data= matrix(Z+mues,ncol=2)
#'
#' # Applying the algorithm
#' sal = ktaucenters_aux(
#' data, K=3, centers=data[sample(1:300,3), ],
#' tolerance=1e-3, max_iter=100)
#'
#' #plot the results
#' plot(data,type='n')
#' points(data[sal$cluster==1,],col=1);
#' points(data[sal$cluster==2,],col=2);
#' points(data[sal$cluster==3,],col=3);
#'
#' @note Some times, if the initial centers are wrong, the algorithm
#'     converges to a non-optimal (local) solution.  To avoid that,
#'     the algorithm must be run several times. This task is carried
#'     out by \code{\link{ktaucenters}}
#' @seealso \code{\link{ktaucenters}}
#' @references Gonzalez, J. D., Yohai, V. J., & Zamar, R. H. (2019).
#'     Robust Clustering Using Tau-Scales. arXiv preprint
#'     arXiv:1906.08198.
#' @export

B1_DEFAULT = 0.5
B2_DEFAULT = 1

ktaucenters_run <-
  function(data, centers, tolerance, max_iter) {
    centers = as.matrix(centers)
    data = as.matrix(data)
    
    # Initialize
    empty_cluster_flag <- FALSE
    n_clusters = ifelse(is.integer(centers), centers, nrow(centers))
    n <- nrow(data)
    p <- ncol(data)
    c1 <- constC1(p)
    c2 <- constC2(p)
    b1 <- B1_DEFAULT
    b2 <- B2_DEFAULT
    tauPath <- NULL
    iter <- 0
    tol <- tolerance + 1
    
    while ((iter < max_iter) & (tol > tolerance)) {
      
      # Step 1: (re)compute labels
      dists <- distance_to_centers(data, centers) 
      distances_min <- dists$min_distance 
      cluster <- dists$membership
      ms <- Mscale(u = distances_min, b = b1, cc = c1) # m-scale
      dnor <- distances_min / ms  # normalized distance
      tau <-
        ms * sqrt(mean(rhoOpt(dnor, cc = c2))) / sqrt(b2) # tau-scale
      tauPath <- c(tauPath, tau)
      
      # Step 2: (re)compute centers
      oldcenters <- centers
      Du <- mean(psiOpt(dnor, cc = c1) * dnor)
      Cu <- mean(2 * rhoOpt(dnor, cc = c2) - psiOpt(dnor, cc = c2) * dnor)
      Wni <- (Cu * psiOpt(dnor, cc = c1) + Du * psiOpt(dnor, cc = c2)) / dnor
      Wni[distances_min == 0] <- (Du * derpsiOpt(0, cc = c2) +
                                    Cu * derpsiOpt(0, cc = c1))
      weights <- Wni / rowsum(Wni, cluster)[cluster] # compute weights
      weights[is.na(weights)] = Wni[is.na(weights)] # FIXME: see original func
      
      # Clusters could have no observations (filled with furthest centers)
      empty_clusters <- tabulate(c(cluster, seq(n_clusters))) == 1
      centers[!empty_clusters,] <- rowsum(data * weights, cluster)
      if (any(empty_clusters)) {
        furthest_indices <- head(order(distances_min, decreasing = TRUE),
                                 sum(empty_clusters))
        centers[empty_clusters, ] <- data[furthest_indices, , drop = FALSE]
        cluster[furthest_indices] <- which(empty_clusters)
        empty_cluster_flag <- TRUE
      }
      tol <- sqrt(sum((oldcenters - centers) ^ 2))
      iter <- iter + 1
    }
    
    return (
      list(
        tauPath = tauPath,
        last_iter = iter,
        tauPath_last_iter = tauPath[iter],
        centers = centers,
        cluster = cluster,
        tol = tol,
        p = p,
        weights = weights,
        distances_min = distances_min,
        Wni = Wni,
        empty_cluster_flag = empty_cluster_flag
      )
    )
  }
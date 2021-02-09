#' ktaucenters_aux
#'
#' Robust Clustering algorithm based on centers, a robust and
#' efficient version of K-Means.
#' @param X A matrix of size n x p.
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
#' \item{\code{emptyClusterFlag}}{: a boolean value. True means
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
#' X= matrix(Z+mues,ncol=2)
#'
#' # Applying the algorithm
#' sal = ktaucenters_aux(
#' X, K=3, centers=X[sample(1:300,3), ],
#' tolerance=1e-3, max_iter=100)
#'
#' #plot the results
#' plot(X,type='n')
#' points(X[sal$cluster==1,],col=1);
#' points(X[sal$cluster==2,],col=2);
#' points(X[sal$cluster==3,],col=3);
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
ktaucenters_aux <- function(X, K, centers, tolerance, max_iter) {

    # Sanitize
    if (!is.matrix(centers)) centers = as.matrix(centers)
    if (!is.matrix(X)) X = as.matrix(X)

    emptyCluster <- FALSE
    emptyClusterFlag <- FALSE

    # Initialize stuff
    n <- nrow(X)
    p <- ncol(X)
    c1 <- constC1(p)
    b1 <- 0.5 # FIXME: avoid hard coded constants
    c2 <- constC2(p)
    b2 <- 1 # FIXME: avoid hard coded constants
    tauPath <- c()
    niter <- 0
    tol <- tolerance + 1
    # repeatedCentersMatrix <- matrix(0L, ncol = p, nrow = n)
    # distances <- matrix(0L, ncol = K, nrow = n)
    
    while ((niter < max_iter) & (tol > tolerance)) {
        
        # Step 1: recompute labels
        
        # Compute distances from observations to each center
        distances <- sqrt(apply(t(centers), 2, 
                                function(m){
                                    colSums((t(X) - m)^2)
                                }))

        ## Find the closest center for each observation
        distances_min <- apply(distances, 1, min)
        
        # Membership vector (when ties, first is chosen)
        cluster <- max.col(-distances, ties = "first")
        
        # Tau-scale
        ms <- Mscale(u = distances_min, b = b1, cc = c1)
        dnor <- distances_min / ms  # normalized distance
        tau <- ms * sqrt(mean(rhoOpt(dnor, cc = c2))) / sqrt(b2)
        tauPath <- c(tauPath, tau)

        # Compute constants in each iteration
        Du <- mean(psiOpt(dnor, cc = c1) * dnor)
        Cu <- mean(2 * rhoOpt(dnor, cc = c2) - psiOpt(dnor, cc = c2) * dnor)
        Wni <- (Cu * psiOpt(dnor, cc = c1) + Du * psiOpt(dnor, cc = c2))/dnor

        # FIXME: di?
        # Attention: when di=0, psi_1(dnor)=0 and psi_2(dnor)=0.
        # Then weight w is undefined due to dividing by zero.
        # Given that psi_1(0)=0, Wni can be obtained
        # through the derivative of psi_1 in this case:
        if (sum(distances_min == 0) > 0) {
            Wni[distances_min == 0] <- (Du * derpsiOpt(0, cc = c2) +
                                        Cu * derpsiOpt(0, cc = c1))
        }

        # Compute weights
        # weights <- numeric(length(Wni))
        # for (jota in 1:K) {
        #     if ((sum(Wni[cluster == jota])) != 0) {
        #         weights[cluster == jota] = Wni[cluster == jota]/sum(Wni[cluster ==
        #           jota])
        #     }
        #     # FIXME: [1] check this piece of code 
        #     if ((sum(Wni[cluster == jota])) == 0) { 
        #         # Esto significa que el algoritmo convergio. Entonces
        #         # no actualizo las X's. Pero si llega a pasar que son
        #         # cero las actualizo:
        #         mmm = length(cluster == jota)
        #         if (sum(weights[cluster == jota] == 0) == mmm)
        #           weights[cluster == jota] = 1/mmm
        #     }
        # }

        # FIXME: check [1]
        weights_aux = lapply(split(Wni, cluster), 
                     function(x) {
                         if (sum(x)) {
                             ret = x / sum(x)}
                         else {
                             ret = ret #1/length(x) 
                         }
                     })
        weights = unsplit(weights_aux, cluster)
        
        # Weight observations
        XW <- crossprod(diag(weights), X)

        # Update centers
        oldcenters <- centers
        
        # Step 2: recompute centers 
        
        # Sometimes a cluster has no observations:
        obs_per_cluster <- table(c(cluster, seq(K))) - 1
        nonempty_cluster <- obs_per_cluster > 0
        centers[nonempty_cluster, ] = Reduce(rbind, 
                                             by(XW, 
                                                cluster,
                                                colSums)) ## FIXME
        
        # if not all clusters are filled, they are replaced for
        # the furthest centers. This is specially important when
        # the number of clusters K is high.
        if (any(!nonempty_cluster)) {
            furtherIndices = head(order(distances_min,
                                        decreasing = TRUE), 
                                  sum(!nonempty_cluster))
            centers[obs_per_cluster == 0, ] = X[furtherIndices, ]
            cluster[furtherIndices] = which(!nonempty_cluster)
            emptyClusterFlag = TRUE
        }

        ## Exit condition
        tol <- sqrt(sum((oldcenters - centers)^2))
        
        ## Next iteration...
        niter <- niter + 1
    }

    return (list(tauPath = tauPath,
                 niter = niter,
                 tauPath_niter = tauPath[niter],
                 centers = centers,
                 cluster = cluster,
                 emptyCluster = emptyCluster,
                 tol = tol,
                 p = p,
                 weights = weights,
                 di = distances_min,
                 Wni = Wni,
                 emptyClusterFlag = emptyClusterFlag))
}
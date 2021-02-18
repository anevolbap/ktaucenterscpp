#' ktaucenters
#'
#' Robust Clustering algorithm based on centers, a robust and
#' efficient version of K-Means.
#' @param X numeric matrix of size n x p.
#' @param K the number of cluster.
#' @param centers a matrix of size K x p containing the K initial
#'     centers, one at each matrix-row. If centers is NULL a random
#'     set of (distinct) rows in \code{X} are chosen as the initial
#'     centres.
#' @param tolerance a tolerance parameter used for the algorithm stopping
#'     rule
#' @param max_iter a maximum number of iterations used for the
#'     algorithm stopping rule
#' @param n_runs the number of trials that the base algorithm
#'     ktaucenters_aux is run.  If it is greater than 1 and center is
#'     not set as NULL, a random set of (distinct) rows in \code{X}
#'     will be chosen as the initial centers.
#' @param startWithKmeans TRUE if kmean centers values is included as
#'     starting point.
#' @param startWithROBINPD TRUE if ROBINDEN estimator is included as
#'     starting point
#' @param flag_outliers optional argument for outliers detection - quantiles
#'     of chi-square to be used as a threshold for outliers detection,
#'     defaults to 0.999

#' @return A list including the estimated K centers and labels for the observations
##' \itemize{
##'  \item{\code{centers}}{:   matrix of size K x p, with the estimated K centers.}
##'  \item{\code{cluster}}{: array of size n x 1  integers labels between 1 and K.}
##'  \item{\code{tauPath}}{: sequence of tau scale values at each iterations.}
##'  \item{\code{Wni}}{: numeric array of size n x 1 indicating the weights
##' associated to each observation.}
##'  \item{\code{emptyClusterFlag}}{: a boolean value. True means that in some
##' iteration there were clusters totally empty}
##'  \item{\code{niter}}{: number of iterations until convergence is achieved
##' or maximun number of iteration is reached}
##'  \item{\code{di}}{: distance of each observation to its assigned cluster-center}
##'  \item{\code{outliers}}{: indices observation that can be considered as outliers}
##' }

#' @examples
#' # Generate Sinthetic data (three cluster well separated)
#' Z <- rnorm(600);
#' mus <- rep(c(-3, 0, 3), 200)
#' X <-  matrix(Z + mus, ncol=2)
#'
#' # Generate 60 sinthetic outliers (contamination level 20%)
#' X[sample(1:300,60), ] <- matrix(runif(40, 3 * min(X), 3 * max(X)),
#'                                 ncol = 2, nrow = 60)
#'
#' ### Applying the algorithm ####
#'sal <- ktaucenters(
#'      X, K=3, centers=X[sample(1:300,3), ],
#'      tolerance=1e-3, max_iter=100)
#'
#' ### plotting the clusters ###
#'
#' oldpar = par(mfrow = c(1, 2))
#'
#' plot(X, type = 'n', main = 'ktaucenters (Robust) \n outliers: solid black dots')
#' points(X[sal$cluster == 1, ], col = 2);
#' points(X[sal$cluster == 2, ], col = 3);
#' points(X[sal$cluster == 3, ], col = 4)
#' points(X[sal$outliers, 1], X[sal$outliers, 2], pch = 19)
#'
#' ### Applying a classical (non Robust) algortihm ###
#' sal <- kmeans(X, centers = 3, n_runs = 100)
#'
#' ### plotting the clusters ###
#' plot(X, type = 'n', main = 'kmeans (Classical)')
#' points(X[sal$cluster == 1, ], col = 2);
#' points(X[sal$cluster == 2, ], col = 3);
#' points(X[sal$cluster == 3, ], col = 4)
#'
#' par(oldpar)
#' @references Gonzalez, J. D., Yohai, V. J., & Zamar, R. H. (2019).
#' Robust Clustering Using Tau-Scales. arXiv preprint arXiv:1906.08198.
#'
#' @importFrom stats kmeans dist qchisq
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar% foreach
#' @export

MAX_ITER_DEFAULT = 100L
N_RUNS_DEFAULT = 1L
TOLERANCE_DEFAULT = 1e-6
N_CORES_DEFAULT = 1L
INIT_CENTERS_DEFAULT = list(quote(init_kmeans), quote(init_robin))
CUTOFF_DEFAULT = 0.999
MSCALE_BP_DEFAULT = 0.5

ktaucenters <- function(data,
                        centers,
                        tolerance = TOLERANCE_DEFAULT,
                        max_iter = MAX_ITER_DEFAULT,
                        n_runs = N_RUNS_DEFAULT,
                        init_centers = INIT_CENTERS_DEFAULT,
                        flag_outliers = outliers_tau_cutoff(CUTOFF_DEFAULT,
                                                            MSCALE_BP_DEFAULT),
                        n_cores = N_CORES_DEFAULT) {
    data = as.matrix(data)
    n_clusters = ifelse(is.list(centers), length(centers), centers)
    
    # Set up center initialization
    add_init_custom = NULL
    if (is.list(centers))
        add_init_custom = quote(init_custom)
    add_init_random = replicate(n_runs, init_random, simplify = FALSE)
    init_centers = append(init_centers, c(add_init_random, add_init_custom), 1)
    
    # Runs
    start_centers = lapply(init_centers, function(x)
        eval(x)(data, n_clusters))
    ktau_runs = lapply(
        start_centers,
        ktaucenters_run,
        data = data,
        tolerance = tolerance,
        max_iter = max_iter
    )
    
    # Pick best run
    best = which.min(lapply(ktau_runs, function(x)
        x$tauPath_last_iter))
    best_ktau = ktau_runs[[best]]
    
    # Outlier detection
    best_ktau = flag_outliers(best_ktau)
    
    return(best_ktau)
}
test_that("ktaucenters", {
  set.seed(6)
  Z <- rnorm(600)
  mues <- rep(c(-4, 0, 4), 200)
  X <- matrix(Z + mues, ncol = 2)
  X[sample(1:300, 60), ] <- matrix(runif(40, 2 * min(X), 2 * max(X)),
                                   ncol = 2, nrow = 60)
  set.seed(6)
  ktau_cpp <- ktaucenterscpp::ktaucenters(data=X, centers = 4, n_runs = 1)
  ktau_cpp$outliers <- sort(ktau_cpp$outliers)
  ktau_cpp$p = NULL
  ktau_cpp$niter = ktau_cpp$last_iter
  ktau_cpp$last_iter = NULL
  ktau_cpp$tauPath_last_iter = NULL
  ktau_cpp$di = ktau_cpp$distances_min
  ktau_cpp$distances_min = NULL
  ktau_cpp$emptyCluster = FALSE

  set.seed(6)
  ktau <- ktaucenters::ktaucenters(X, K = 4, nstart = 1)  
  ktau$outliers <- sort(ktau$outliers)
  
  expect_equal(ktau_cpp[names(ktau)], ktau, 1e-5)
  
})


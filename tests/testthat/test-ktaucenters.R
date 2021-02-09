test_that("ktaucenters", {
  set.seed(6)
  Z <- rnorm(600)
  mues <- rep(c(-4, 0, 4), 200)
  X <- matrix(Z + mues, ncol = 2)
  X[sample(1:300, 60), ] <- matrix(runif(40, 2 * min(X), 2 * max(X)),
                                   ncol = 2, nrow = 60)
  set.seed(6)
  ktau_cpp <- ktaucenterscpp::ktaucenters(X, K = 4, nstart = 1)
  ktau_cpp$outliers <- sort(ktau_cpp$outliers)
  ktau_cpp$p = NULL
  ktau_cpp$tauPath_niter = NULL
  set.seed(6)
  ktau <- ktaucenters::ktaucenters(X, K = 4, nstart = 1)  
  ktau$outliers <- sort(ktau$outliers)
  expect_equal(ktau_cpp, ktau)
  
})
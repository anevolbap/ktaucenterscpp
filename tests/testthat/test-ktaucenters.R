set.seed(6)
Z <- rnorm(600)
mues <- rep(c(-4, 0, 4), 200)
X <- matrix(Z + mues, ncol = 2)
X[sample(1:300, 60), ] <- matrix(runif(40, 2 * min(X), 2 * max(X)),
                                 ncol = 2, nrow = 60)

test_that("ktaucenters-centers", {
  set.seed(6)
  ktau_cpp <- ktaucenterscpp::ktaucenters(data=X, centers = 4, n_runs = 1)
  set.seed(6)
  ktau <- ktaucenters::ktaucenters(X, K = 4, nstart = 1)
  
  expect_equal(ktau_cpp$centers, ktau$centers, tolerance = 1e-6)
})

test_that("ktaucenters-outliers", {
  set.seed(6)
  ktau_cpp <- ktaucenterscpp::ktaucenters(data=X, centers = 4, n_runs = 1)
  set.seed(6)
  ktau <- ktaucenters::ktaucenters(X, K = 4, nstart = 1)
  
  expect_equal(sort(ktau_cpp$outliers), sort(ktau$outliers))
})

test_that("ktaucenters-Wni", {
  set.seed(6)
  ktau_cpp <- ktaucenterscpp::ktaucenters(data=X, centers = 4, n_runs = 1)
  set.seed(6)
  ktau <- ktaucenters::ktaucenters(X, K = 4, nstart = 1)
  
  expect_equal(sort(ktau_cpp$Wni), sort(ktau$Wni), tolerance = 1e-6)
})

test_that("ktaucenters-weights", {
  set.seed(6)
  ktau_cpp <- ktaucenterscpp::ktaucenters(data=X, centers = 4, n_runs = 1)
  set.seed(6)
  ktau <- ktaucenters::ktaucenters(X, K = 4, nstart = 1)
  
  expect_equal(sort(ktau_cpp$weights), sort(ktau$weights), tolerance = 1e-6)
})
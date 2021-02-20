test_that("robust-init-density", {
  set.seed(6)
  Z <- rnorm(600)
  mues <- rep(c(-4, 0, 4), 200)
  X <- matrix(Z + mues, ncol = 2)
  X[sample(1:300, 60), ] <- matrix(runif(40, 2 * min(X), 2 * max(X)),
                                   ncol = 2, nrow = 60)
  k = 5
  D = dist(X)
  
  set.seed(6)
  rob_cpp <- ktaucenterscpp::robust_init_density(D, X, k)
  rob_cpp$id.means = rob_cpp$id_means
  rob_cpp$id_means = NULL
  rob_cpp$idpoints = rob_cpp$id_points
  rob_cpp$id_points = NULL
  
  set.seed(6)
  rob <- ktaucenters::ROBINDEN(D, X, k = k)
  
  expect_equal(rob, rob_cpp)
})
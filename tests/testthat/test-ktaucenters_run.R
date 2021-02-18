test_that("ktaucenters_run", {
    set.seed(6)
    Z <- rnorm(600)
    mues <- rep(c(-4, 0, 4), 200)
    X <- matrix(Z + mues, ncol = 2)
    X[sample(1:300, 60), ] <- matrix(runif(40, 2 * min(X), 2 * max(X)),
                                     ncol = 2, nrow = 60)
    set.seed(6)
    kaux_cpp = ktaucenterscpp::ktaucenters_run(X, centers = diag(4)[, 1:2],
                              tolerance=1e-7,max_iter=1000)

    set.seed(6)
    kaux = ktaucenters::ktaucenters_aux(X, K = 4, centers = diag(4)[, 1:2],
                              tolmin=1e-7, NiterMax=1000)
    
    expected_output = c(0.002063637, 0.013249235, 0.012024533,
                        0.006088092, 0.013249235, 0.012024533)


    expect_equal(kaux_cpp$weights, kaux$weights)
})

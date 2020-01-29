test_that("denpoints", {

    set.seed(1)
    m <- matrix(rnorm(1e3), ncol = 2)
    expected_output <- c(5.990406, 4.711414, 4.134972, 4.651841, 5.580436, 3.178011)

    expect_equal(round(head(denpoints(m)), 6),
                 round(expected_output, 6))
})

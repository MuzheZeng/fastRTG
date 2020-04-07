test_that("multiplication works", {
  n=10
  G = array(abs(rnorm(27)),dim = rep(3,3))
  Pi = c(0.3,0.3,0.4)
  sampleSBM = sbmt(n, Pi, G, sparsity = 0.01)
  expect_equal(length(sampleSBM@vals), 10, tolerance = 10)
})

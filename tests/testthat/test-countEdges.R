test_that("multiplication works", {
  set.seed(123)
  n = 10
  k = 3
  X = matrix(rnorm(n*k),n,k)
  X = list(X,X,X)
  G = rTensor::as.tensor(array(rexp(k^3),dim = rep(k,3)))
  res = countEdges(X, G)
  expect_equal(res[2], 0.08266169,tolerance = 1e-5)
})

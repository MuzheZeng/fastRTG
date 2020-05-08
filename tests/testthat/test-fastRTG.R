test_that("sparsity of fastRTG", {
  G = array(rgamma(120,1,1), dim = c(2,3,4,5))
  X = list()
  X[[1]] = matrix(abs(rnorm(20)),10,2)
  X[[2]] = matrix(rgamma(60,1,1),20,3)
  X[[3]] = matrix(rf(120,3,4),30,4)
  X[[4]] = matrix(1,40, 5)
  sampleTensor <- fastRTG(X, G, sparsity = 0.01, returnParameters = TRUE)
  expect_equal(sum(sampleTensor$tnsr@vals)/prod(sampleTensor$tnsr@dims), 0.01, tolerance = 0.005)
})

test_that("output datatype", {
  G = array(rgamma(24,1,1), dim = c(2,3,4))
  X = list()
  X[[1]] = matrix(abs(rnorm(20)),10,2)
  X[[2]] = matrix(rgamma(60,1,1),20,3)
  X[[3]] = matrix(rf(120,3,4),30,4)
  sampleTensor <- fastRTG(X, G, sparsity = 0.1, returnParameters = TRUE)
  expect_equal(class(sampleTensor$tnsr)[1], "sptensor")
  expect_equal(class(sampleTensor$core)[1], "Tensor")
  expect_true(class(sampleTensor$list_mat)[1] %in% c("list", "matrix"))
})


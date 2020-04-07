#' Count the edges information of an initial
#' random tensor graph.
#'
#' This function counts
#'
#' @param X A list of matrices
#' @param G The core tensor. A rTensor::Tensor object.
#'
#' @return a 2-dim vector describing sum of tensor values and average value.
#' @export
#'
#' @examples
#' n = 10
#' k = 3
#' X = matrix(rnorm(n*k),n,k)
#' X = list(X,X,X)
#' G = rTensor::as.tensor(array(rexp(k^3),dim = rep(k,3)))
#' res = countEdges(X, G)
countEdges <- function(X, G) {
  M = length(attr(G, "modes"))
  if(length(X) == 1) {
    tmp = list(X[[1]])
    for(i in 1:(M-1)){
      tmp[[i+1]] = X[[1]]
    }
    X = tmp
  } else if (length(X) != M) {
    stop("Numbers of modes doesn't match!")
  }
  Cx = list()
  for(i in 1:M){
    Cx[[i]] = diag(colSums(X[[i]]), nrow = ncol(X[[i]]), ncol = ncol(X[[i]]))
  }
  em = sum(attr(rTensor::ttl(G, Cx, 1:M),"data"))
  avgDeg = em/prod(as.numeric(lapply(X, function(a) nrow(a)[1])))
  return(c(em, avgDeg))
}

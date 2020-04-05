#' Count the edges information of an initial
#' random tensor graph.
#'
#' @param X A list of matrices
#' @param G The core tensor
#'
#' @return a 2-dim vector describing average degree
#' @export
#'
#' @examples
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
  em = sum(attr(ttl(G, Cx, 1:M),"data"))
  avDeg = em/nrow(X[[1]])
  return(c(em, avDeg))
}

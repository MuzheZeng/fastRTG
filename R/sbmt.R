#' Sampling Stochstic Tensor Block Models.
#'
#' Suppose the tensor \eqn{A\in \mathbb{R}^{n_1\times n_2...\times n_m} has expectation structure:
#' \deqn{\mathbb{E}(A) = [G;Z_1,Z_2,Z_3,...Z_m].}
#'
#' @param n
#' @param Pi
#' @param G
#' @param PoissonEdges
#' @param avgDeg
#' @param returnParameters
#' @param parametersOnly
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
sbmt <- function(n, Pi, G, PoissonEdges = FALSE, avgDeg = NULL, returnParameters = FALSE, parametersOnly = FALSE, ...) {
  K = length(Pi)
  if (K != dim(G)[1] || length(unique(dim(G))) != 1) {
    stop("Core tensor should be super-diagonal with every dimension match length of Pi!")
  }
  Z = sample(K, n, replace = TRUE, prob = Pi)
  Z = sort(Z)
  X = model.matrix(~factor(as.character(Z), levels = as.character(1:K)) - 1)

  if (length(avgDeg) == 0) {
    return(list(tensor=fastRTG(list(X), G, PoissonEdges = PoissonEdges, ...), Z = X))
  } else {
    eDbar = countEdges(list(X), G)[2]
    G = G * avgDeg/eDbar
  }

  if (!PoissonEdges) {
    if (max(G@data) >= 1) {
      warning("This combination of B and avgDeg has led to probabilities that exceed 1.
              Suggestion:  Either diminish avgDeg or enable poisson edges.")
    } else {
      G@data = -log(1-G@data)
    }
  }

  if (parametersOnly) return(list(X = list(X), core=G))
  return(list(tensor=fastRTG(X = list(X),G, avgDeg = avgDeg, PoissonEdges = PoissonEdges, returnParameters = returnParameters, ...),Z=X))
}

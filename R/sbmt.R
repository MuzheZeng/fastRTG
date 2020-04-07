#' Sampling Stochstic Tensor Block Models.
#'
#' Suppose the binary (Poisson) tensor A has expectation structure:
#' \deqn{E(A) = [G;Z_1,Z_2,Z_3,...Z_m],}
#' where each row of Z_i's has exactly one non-zero element with value equals to 1.
#'
#'
#' @param n a vector of positive integers specifying the dimensions sizes of each mode.
#' @param Pi  a vector with same size as the number of clusters specifying the cluster weights.
#' @param G the core tensor in a multi-dimensional array.
#' @param PoissonEdges  boolean indicator. If TRUE, elements of A allows multiple same edges. IF FALSE, A should be a binary tensor.
#' @param avgDeg  an integer specifying the expected degree.
#' @param returnParameters  logical. Return the parameters or not.
#' @param parametersOnly  logical. Only return the parameters or not.
#' @param ... other parameters.
#'
#' @return A list of three items: The sampled random tensor. The latent factors Zs. The ground truth core tensor G.
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
  return(fastRTG(X = list(X),G, avgDeg = avgDeg, PoissonEdges = PoissonEdges, returnParameters = returnParameters, ...))
}

#' Title
#'
#' @param X a list of matrices X_i's. Each X_i is a \eqn{n_i by k_i} matrix.
#' n_i is the dimension size in the i-th mode. k_i is the number of latent
#' factors in the i-th mode. If X only has one matrix, the default reads it
#' as a 3-mode tensor with same X_1 in all three modes.
#' @param G the core tensor. Should not contain negative values.
#' @param avgDeg  specifies the expected degree.
#' @param PoissonEdges  boolean indicator. Allow poisson edges if TRUE, otherwise only binary edges.
#' @param returnParameters  return parameter list or not.
#'
#' @return
#' @export
#'
#' @examples
fastRTG <- function(X, G, avgDeg = NULL, PoissonEdges = TRUE, returnParameters = FALSE) {
  M = length(attr(G, "modes"))
  if(length(X) == 1) {
    X = list(X[[1]], X[[1]], X[[1]])
  } else if (length(X) != M) {
    stop("Numbers of modes doesn't match!")
  }

  for(i in 1:M) {
    if (sum(X[[i]]<0)) {
      stop(paste("Negative entries found in",i,"th matrix.", sep=" "))
    }
  }
  if(sum(G@data<0)) {
    stop("Negative entries found in core tensor.")
  }

  # if (length(simple)>0) if(simple) {
  #   selfLoops = directed = PoissonEdges = FALSE
  # }

  N = as.numeric(lapply(X, nrow))
  K = G@modes

  if(length(avgDeg)>0) {
    eDbar = countEdges(X, G)[2]
    G = G * avgDeg / eDbar
  }

  Cx = list()
  for(i in 1:M){
    Cx[[i]] = diag(colSums(X[[i]]), nrow = ncol(X[[i]]), ncol = ncol(X[[i]]))
  }

  Gt = ttl(G, Cx, 1:M)
  m = rpois(n = 1, lambda = sum(Gt@data))

  if (m == 0) {
    if (returnEdgeList) return(matrix(NA, nrow = 0, ncol = m))
    subs = list(rep(1,m), N)
    vals = 0
    dims = N
    tau = sptensor(subs, vals, dims)
    if(returnParameters){
      out = list(tnsr = tau, list_mat = X, core = G)
    }
    if(!returnParameters){
      out = tau
    }
    return(out)
  }


  tabUVW = attr(as.tensor(array(rmultinom(n=1, size=m, prob=Gt@data), dim=K)), "data")
  cumsumUVW = c(0,as.numeric(cumsum(tabUVW)))

  edges = matrix(NA, nrow = m, ncol = M)
  for(k in 1:prod(K)){
    curm = tabUVW[k]

    if(curm > 0) {
      index = matrix(0, ncol = M)
      index[1] = k %% K[1]
      index[1] = ifelse(index[1] == 0, K[1], index[1])
      for(i in 2:M){
        index[i] = ceiling(k/prod(K[1:(i-1)])) %% K[i]
        index[i] = ifelse(index[i] == 0, K[i], index[i])
      }

      for(j in 1:M){
        edges[(cumsumUVW[k]+1):cumsumUVW[k+1],j] = sample(1:N[j], curm, replace=TRUE, prob = X[[j]][,index[j]])
      }
    }
  }

  if (dim(edges)[1] == 0) {
    if (returnEdgeList) return(matrix(NA, nrow = 0, ncol = M))
    subs = list(rep(1,M), N)
    vals = 0
    dims = N
    tau = sptensor(subs, vals, dims)
    if(returnParameters){
      out = list(tnsr = tau, list_mat = X, core = G)
    } else{
      out = tau
    }
    return(out)
  }

  elcount <- table(apply(edges, 1, paste, collapse = ","))
  nms <- names(elcount)
  nms <- strsplit(nms, ",")
  nms <- lapply(nms, as.numeric)
  subs <- as.data.frame(nms)
  vals = elcount
  if (PoissonEdges) {
    tau = sptensor(subs, vals, N)
  }else {
    vals[vals>1] = 1
    tau = sptensor(subs, vals, N)
  }

  if (returnParameters) {
    return(list(tnsr = tau, list_mat = X, core = G))
  } else {
    return(tau)
  }

}

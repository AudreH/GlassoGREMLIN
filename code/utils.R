#' @importFrom igraph as_adj
#' @example 
#' Lambda <- matrix(c(1, 2, 3, 2, 1, 4, 3, 4, 1), 3, 3)
#' p <- 20
#' alpha <- c(1/2, 1/4, 1/4)
#' mySBM <- rSBMLaplace(p, Lambda, alpha)
#' cl <- V(mySBM)$membership
#' Omega <- mySBM %>% as_adj(attr = "weight")
#' X <- rMVTnorm(100, mySBM)
rMVTnorm <- function(n, SBMLaplace, mu = rep(0, ncol(Omega))) {
  Omega <- as_adj(SBMLaplace, attr = "weight")
  diag(Omega) <- colSums(Omega)
  Sigma <- as.matrix(chol2inv(chol(Omega)))
  res <- rmvnorm(n, mean = mu, sigma = Sigma)
  res
}

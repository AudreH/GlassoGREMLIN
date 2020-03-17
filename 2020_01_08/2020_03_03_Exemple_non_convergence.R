setwd("/home/hulot/Documents/packages_R/GlassoGREMLIN/2020_01_08/")
# .libPaths("/projet/gabi/save/ahulot/R_lib/")
# setwd("/projet/gabi/save/ahulot/Documents/GRGL/")

## ---- PACKAGES : -----

library(Matrix)
library(GREMLIN)
library(aricode)
library(glassoFast)
library(blockmodels)

library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggpubr)

library(MASS)

library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

library(igraph)

# ---- SOURCE FONCTIONS : -----

source("Code_Sophie/functions_plot.R")

# source("Functions.R") # CVglasso, pas forcement utile ici.
# source("Functions_boucles.R")

## ----- CODE JULIEN : (SIMULATION DONNEES) ----
#' rSBM_Laplace
#'
#' Generate a symmetric matrix under the stochastic block model
#' with a Laplace distribution as the emission law of the edges.
#'
#' @param p number of nodes
#' @param Lambda a QxQ matrix of connectivity parameters
#' @param alpha a size-Q vector of mixture
#' @return an igraph object: undirected weigthed graph with a attribute "membership" for the vertices
#' @import Matrix igraph
#' @importFrom rmutil rlaplace

rSBM_Laplace <- function(p, Lambda, alpha = c(1), Mean = NULL) {
  
  Q <- length(alpha)
  Z <- t(rmultinom(p, 1, prob = alpha))
  A <- matrix(0, p, p)
  subset <- upper.tri(A)
  lambda <- (Z %*% Lambda %*% t(Z))[subset]
  if(is.null(Mean)){Mean = matrix(0, ncol(Lambda), ncol(Lambda))}
  means <- (Z %*% Mean %*% t(Z))[subset]
  x <- rmutil::rlaplace(n=p*(p-1)/2, m = means, s = lambda)
  A[subset] <- x
  A <- A + t(A)
  diag(A) <- colSums(abs(A)) # diagonale dominante pour inversible
  mySBM <- graph_from_adjacency_matrix(A, weighted = TRUE, mode = "undirected")
  vertex_attr(mySBM, "membership") <- Z %*% 1:Q
  mySBM
}

Sigma_From_Omega = function(Omega){
  Sigma <- as.matrix(chol2inv(chol(Omega)))
  return(Sigma)
}

rmvnorm_SBM <- function(n, sbm) {
  Omega <- as_adj(sbm, attr = "weight")
  Sigma = Sigma_From_Omega(Omega)
  res <- mvtnorm::rmvnorm(n, sigma = Sigma)
  res
}

## ---- SETTINGS POUR SIMULATION : ----

n_ind = 5000
n_feat = 30
nb_groupes = 3
n_iter = 10 # Nombre max d'iterations, pas forcement atteint

alpha. = rep(1/nb_groupes, nb_groupes)
round(alpha., 2)

Lambda. = matrix(0.001, nb_groupes, nb_groupes)
diag(Lambda.) = c(0.01, 0.10, 0.01)

Means. = matrix(0, nb_groupes, nb_groupes)

set.seed(1992)
Sim_SBM = rSBM_Laplace(n_feat, Lambda., alpha., Means.)
dat = rmvnorm_SBM(n = n_ind, Sim_SBM)

classification = V(Sim_SBM)$membership

# ---- Visualisation données simulées : ----

Omega_sim = as_adjacency_matrix(Sim_SBM, attr = "weight", sparse = FALSE)
Sigma_sim = Sigma_From_Omega(Omega_sim)

g_sig_sim = image(Matrix(Sigma_sim)[order(classification), order(classification)])
g_om_sim = image(Matrix(Omega_sim)[order(classification), order(classification)])
g_cov = image(Matrix(cov(dat))[order(classification), order(classification)])
g_inv = image(Matrix(solve(cov(dat)))[order(classification), order(classification)])

ggarrange(g_sig_sim, g_om_sim, g_cov, g_inv, labels = c("Sigma Sim", "Omega sim", "Cov dat", "Inv Cov"), nrow = 2, ncol = 2)

# ---- GREMLIN : (multipartiteBM error) ----

gl1 = glassoFast(S = cov(dat), rho = matrix(0.01, ncol = ncol(cov(dat)), nrow = nrow(cov(dat))))

g_w = image(Matrix(gl1$w[order(classification), order(classification)]))
g_wi = image(Matrix(gl1$wi[order(classification), order(classification)]))
ggarrange(g_sig_sim, g_om_sim, g_w, g_wi, labels = c("Sigma Sim", "Omega sim", "GL cov", "GL Inv"), nrow = 2, ncol = 2)

Net = defineNetwork(gl1$w, typeInter = "adj", rowFG = "Features", colFG = "Features")

gr1 = multipartiteBM(list_Net = list(Net), namesFG = "Features", v_distrib = "gaussian", v_Kmin = 2, v_Kmax = 3, initBM = FALSE)
gr1 = multipartiteBM(list_Net = list(Net), namesFG = "Features", v_distrib = "gaussian", v_Kmin = 2, v_Kmax = 3, initBM = FALSE, maxiterVE = 10000)

# ---- GREMLIN : (multipartiteBMFixedModel error) : ----

gl1 = glassoFast(S = cov(dat), rho = matrix(0.1, ncol = ncol(cov(dat)), nrow = nrow(cov(dat))))

g_w = image(Matrix(gl1$w[order(classification), order(classification)]))
g_wi = image(Matrix(gl1$wi[order(classification), order(classification)]))
ggarrange(g_sig_sim, g_om_sim, g_w, g_wi, labels = c("Sigma Sim", "Omega sim", "GL cov", "GL Inv"), nrow = 2, ncol = 2)

Net = defineNetwork(gl1$w, typeInter = "adj", rowFG = "Features", colFG = "Features")
gr1 = multipartiteBMFixedModel(list_Net = list(Net), namesFG = "Features", v_distrib = "gaussian", v_K = 3)


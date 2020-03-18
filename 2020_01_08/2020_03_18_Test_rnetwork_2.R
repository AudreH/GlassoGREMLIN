# Changement fonction rnetWork pour test

rm(list = ls())
setwd("/home/hulot/Documents/packages_R/GlassoGREMLIN/2020_01_08/")

## ---- PACKAGES : -----

library(Matrix)
library(GREMLIN)
library(aricode)
library(glassoFast)
library(blockmodels)
library(simone)

library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggpubr)

library(MASS)

library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

library(igraph)
library(ggraph)

library(parallel)

# ---- SOURCE FONCTIONS : -----

source("Code_Sophie/functions_plot.R")
source("Function_graph.R")
source("Functions_boucles.R")

# ---- Fonctions simone: ----

# Laplacian, reprise simone
laplacian = function (M) 
{
  n <- nrow(M)
  diag(M) <- 1
  D <- colSums(M)
  L <- diag(rep(1, n)) - diag(D^(-1/2)) %*% M %*% diag(D^(-1/2))
  L[is.nan(L)] <- 0
  return(L)
}

# Changement pour signed & ajout rnorm
rNetwork2 =   function (p, pi, alpha = c(1), sdmat = NULL, directed = FALSE, name = "a network",
                         signed = TRUE)
{
  A <- matrix(0, p, p)
  Q <- length(alpha)
  pi <- as.matrix(pi)
  fixed.edges <- FALSE
  
  if(is.null(sd)) sdmat = matrix(1, ncol = length(alpha), nrow = length(alpha))
  
  if (max(pi) > 1) {
    if (sum(pi)%%floor(sum(pi)) == 0) {
      fixed.edges <- TRUE
    }
  }
      pi.bounds <- floor((p * alpha) %*% t(p * alpha))
      
    if (any(pi > pi.bounds)) {
      cat("\nMatrix Pi is too dense: here is the adjusted version\n")
      pi[pi > pi.bounds] <- pi.bounds[pi > pi.bounds]
      print(pi)
    }
  
  cl.levels <- factor(1:ncol(pi))
  if (is.null(rownames(pi))) {
    colnames(pi) <- rownames(pi) <- 1:ncol(pi)
  }
  names(alpha) <- cl.levels

    clusters <- factor(sample(cl.levels, p, prob = alpha,
                              replace = TRUE), levels = cl.levels)
  
  couple <- expand.grid(clusters, clusters)
  if (directed) {
    possible.edges <- rep(TRUE, nrow(couple))
  }
  else {
    up <- upper.tri(A, diag = FALSE)
    possible.edges <- rep(TRUE, nrow(couple)) & up
  }
  for (q in 1:Q) {
    for (l in 1:Q) {
      ql <- which(couple[, 1] == q & couple[, 2] == l &
                    possible.edges)
        A[ql] <- rbinom(length(ql), 1, pi[q, l]) * rnorm(length(ql), 0, sd = sdmat[q,l])
    }
  }
  if (!directed) {
    A <- (A | t(A)) * 1
  }
  dimnames(A) <- list(clusters, clusters)
  if (signed) {
    switchsign <- matrix(sample(c(1, -1), p * p, replace = TRUE),
                         p, p)
  }
  else {
    switchsign <- matrix(1, p, p) ## Modif here
  }
  eig = eigen(A)
  A = eig$vectors %*% diag(abs(eig$values)) %*% t(eig$vectors)
  Theta = -switchsign * t(switchsign) * A
  # Theta = A
  # if (!directed) {
  #   L <- laplacian(A)
  #   Theta <- -switchsign * t(switchsign) * L
  #   diag(Theta) <- rep(1, p)
  # }
  # else {
  #   Theta <- A
  #   corr <- runif(p * p)
  #   Theta[Theta != 0] <- corr[rank(corr) > sum(A == 0)] *
  #     ((switchsign * A)[A != 0])
  #   if (!(max(Mod(eigen(Theta)$values)) < 1)) {
  #     Theta <- Theta/(max(Mod(eigen(Theta)$values)) + 1e-09)
  #   }
  # }
  nodes <- as.character(paste("g", 1:p, sep = ""))
  dimnames(Theta) <- list(nodes, nodes)
  return(structure(list(A = A, Theta = Theta, directed = directed,
                        clusters = clusters, name = name), class = "simone.network"))
}


# ---- TESTS RNETWORK : -----
p = 30
nb_groupes = 3
pi = matrix(0, nb_groupes, nb_groupes)
diag(pi) = c(1, 1, 1)
alpha = rep(1/nb_groupes, nb_groupes)

sdmat =  matrix(0, nb_groupes, nb_groupes)
diag(sdmat) = c(0.5, 0.5, 0.5)

set.seed(1992)
net_1 = rNetwork2(p = p,
                 pi = pi,
                 alpha    = alpha,
                 sdmat = sdmat,
                 directed = FALSE, 
                 name     = "net_1")

Omega = net_1$Theta
Sigma = solve(net_1$Theta)
Sigma_chol = chol2inv(net_1$Theta) # verif que ca marche, mais laplacien = ok de toute facon
classification = net_1$clusters

col2 = list("Cl" = c("1" = "blue", "2" = "azure3", "3" = "darkorchid3", "4" = "darksalmon", "5" = "darkolivegreen3"))
ha = HeatmapAnnotation(df = data.frame("Cl" = as.factor(classification)), col = col2)
ha_r = HeatmapAnnotation(df = data.frame("Cl" = as.factor(classification)), which = "row", col = col2)

draw(Heatmap(Omega - diag(diag(Omega)), top_annotation = ha, right_annotation = ha_r, name = "Omega",
             row_order = order(classification), column_order = order(classification)), merge_legend = TRUE)
draw(Heatmap(Sigma - diag(diag(Sigma)), top_annotation = ha, right_annotation = ha_r, name = "Sigma",
             row_order = order(classification), column_order = order(classification)), merge_legend = TRUE)


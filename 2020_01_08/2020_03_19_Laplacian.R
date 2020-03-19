rm(list = ls())
library(igraph)
library(corrplot)

p  <- 120
sizes <- c(50, 30, 30, 10)
K <- length(sizes)
PI <- diag(.4, K) + .05 # connectivity matrix
mySBM <- igraph::sample_sbm(p, PI, sizes)
mySBM <- set_vertex_attr(mySBM, "cluster", value = rep(1:K, sizes))
par(mfrow= c(1,1))
plot(mySBM, vertex.color = V(mySBM)$cluster)

## random Gaussian weights
mu <- 2; sigma <- 1
mySBM <- set_edge_attr(mySBM, "weight", value = rnorm(gsize(mySBM), mu , sigma))

epsilon <- 1e-1
Theta <- laplacian_matrix(mySBM, normalized = TRUE)
Theta <- Theta + diag(epsilon, p, p)
Sigma <- solve(Theta)

par(mfrow = c(2,2))
Theta_plot <- as.matrix(Theta); diag(Theta_plot) <- NA
Sigma_plot <- as.matrix(Sigma); diag(Sigma_plot) <- NA
corrplot(Theta_plot, method = "color", is.corr = FALSE, tl.pos = "n")
corrplot(Sigma_plot, method = "color", is.corr = FALSE, tl.pos = "n")
hist(Theta_plot)
hist(Sigma_plot)

source("utils_laplace.R")
source("VEM_SBMLaplace.R")

Lambda <- matrix(c(1, 2, 3, 2, 1, 4, 3, 4, 1), 3, 3)
p <- 200
alpha <- c(1/2, 1/4, 1/4)

mySBM <- rSBMLaplace(p, Lambda, alpha)

cl <- V(mySBM)$membership
Omega <- mySBM %>% as_adj(attr = "weight")
image(Omega[order(cl), order(cl)])

print(aricode::ARI(kmeans(Omega, centers = 3, nstart = 10)$cl , cl))

out <- VEM_SBM_laplace(mySBM, 3)

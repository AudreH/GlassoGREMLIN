source("utils_laplace.R")
source("VEM_SBMLaplace.R")

Lambda <- matrix(c(1, 2, 3, 2, 1, 4, 3, 4, 1), 3, 3)
p <- 200
alpha <- c(1/2, 1/4, 1/4)

mySBM <- rSBMLaplace(p, Lambda, alpha)
Omega <- mySBM %>% as_adj(attr = "weight")



# HOW to charge the diagonal  : sol 1 : using the principal minors...  numerical pob . 
# eta <- 1;
# Omega[1,1] = rexp(1,eta)
# Omega <- as.matrix(Omega)
# for (i in 2:nrow(Omega)){
#       U <- rexp(1,eta)
#       Binf <-  0 
#       Omegaii  = as.matrix(Omega[1:i,1:i])
#       for (j in 1:(i-1)){
#         Binf <- Binf + (-1)^(i + j)*Omegaii[i,j]*det(as.matrix(Omegaii[-i,-j]))
#         print(c(j,Binf))
#       }
#       
#       
#       Omega[i,i] <- Binf / ((-1)^{i+i}*det(as.matrix(Omegaii[-i,-i]))) + U  
# }
# 

# HOW to charge the diagonal  : sol 2 : using the charge the diagonal  : very large values on the diagonal 
eta <- 1;
diag(Omega) = 0
for (i in 1:nrow(Omega)){
  eps <- rexp(1)  
  Omega[i,i]  = eps + sum(abs(Omega[i,]))
}
min(eigen(Omega)$values)





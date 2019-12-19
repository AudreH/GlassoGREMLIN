require(R.utils)

# ---------------- GLGR_FUNC : ----

glgr_func = function(X, s, rho,
                     thr,
                     maxit,
                     penalize.diagonal,
                     start = "warm",
                     w.init, wi.init,
                     trace = FALSE,
                     nb_feat, 
                     blocks_names,
                     v_Kmin, 
                     v_Kmax, 
                     verbose = FALSE,
                     # maxiterVE = 100,
                     nbCores = 1,
                     V = NULL
                  
){
  # cat("I'm here ! \n")
  
  initial.guess = NULL
  Sigma0 = NULL
  if(start == "warm"){initial.guess = wi.init ; Sigma0 = w.init}
  
  Rho = NULL
  if(is.null(V)){
    Rho = matrix(rho, ncol(X), ncol(X))
  }else{
    Rho = lam*V
  }
  
  GLASSO = graphical.lasso(X = X, 
                           Sigma = Sigma0,
                           Rho = Rho, 
                           initial.guess = initial.guess,
                            eps = 1e-6, maxIt = maxit) 
  
  if(any(is.na(GLASSO$Theta))){
    GLASSO = graphical.lasso(X = X, 
                    Rho = matrix(rho, ncol(X), ncol(X)), 
                    initial.guess = initial.guess,
                    eps = 1e-6, maxIt = maxit) 
  }
  
  A_GLASSO <- ifelse(GLASSO$Theta!=0 & row(GLASSO$Theta)!=col(GLASSO$Theta),1,0)
  
  GREM = gremlin_list(A_GLASSO,  nb_feat = nb_feat, 
                      blocks_names = blocks_names,
                      v_Kmin = v_Kmin, v_Kmax = v_Kmax, 
                      verbose = verbose, 
                      # maxiterVE = maxiterVE,
                      nbCores = nbCores
                      )

  GLASSO = graphical.lasso(X = X,
                           Sigma = GLASSO$Sigma,
                           Rho = rho*(1/GREM$matrix_reconstruct), 
                           initial.guess = GLASSO$Theta,
                           eps = 1e-6, maxIt = maxit) 
  
  if(any(is.na(GLASSO$Theta))){
    GLASSO = graphical.lasso(X = X,
                             # Sigma = GLASSO$Sigma,
                             Rho = rho*(1/GREM$matrix_reconstruct), 
                             initial.guess = GLASSO$Theta,
                             eps = 1e-6, maxIt = maxit) 
  }
  
  ## Note sure that's useful
  A_GLASSO2 = ifelse(GLASSO$Theta!=0 & row(GLASSO$Theta)!=col(GLASSO$Theta),1,0)
  diff = sum((A_GLASSO - A_GLASSO)!=0)
  iter = 1
  while(diff>0.05*length(A_GLASSO) & iter < 10){
    
    A_GLASSO <- ifelse(GLASSO$Theta!=0 & row(GLASSO$Theta)!=col(GLASSO$Theta),1,0)
    
    GREM = gremlin_list(A_GLASSO,  nb_feat = nb_feat, 
                        blocks_names = blocks_names,
                        v_Kmin = v_Kmin, v_Kmax = v_Kmax, 
                        verbose = verbose, 
                        # maxiterVE = maxiterVE,
                        nbCores = nbCores
    )
    
    GLASSO = graphical.lasso(X = X,
                             Sigma = GLASSO$Sigma,
                             Rho = rho*(1/GREM$matrix_reconstruct), 
                             initial.guess = GLASSO$Theta,
                             eps = 1e-6, maxIt = maxit) 
    
    if(any(is.na(GLASSO$Theta))){
      GLASSO = graphical.lasso(X = X,
                               # Sigma = GLASSO$Sigma,
                               Rho = rho*(1/GREM$matrix_reconstruct), 
                               initial.guess = GLASSO$Theta,
                               eps = 1e-6, maxIt = maxit) 
    }
    
    A_GLASSO2 = ifelse(GLASSO$Theta!=0 & row(GLASSO$Theta)!=col(GLASSO$Theta),1,0)
    diff = sum((A_GLASSO - A_GLASSO)!=0)
    iter = iter+1
  }
  
  GLASSO$wi = GLASSO$Theta
  GLASSO$w = GLASSO$Sigma
  
  return(list(GLASSO = GLASSO, GREM = GREM))
}


# GLASSO = glasso(s = S.train, rho = lam_, thr = tol,
#                 maxit = maxit, penalize.diagonal = diagonal,
#                 start = start, w.init = init, wi.init = initOmega,
#                 trace = FALSE, ...)
# 
# A_GLASSO <- ifelse(GLASSO$wi!=0 & row(GLASSO$wi)!=col(GLASSO$wi),1,0)
# 
# GREM = gremlin_list(A_GLASSO,  nb_feat = nb_feat,
#                     blocks_names = blocks_names,
#                     v_Kmin = v_Kmin, v_Kmax = v_Kmax,
#                     verbose = verbose)
# 
# GLASSO = glasso(s = S.train, rho = lam_*(1/GREM$matrix_reconstruct), thr = tol,
#                 maxit = maxit, penalize.diagonal = diagonal,
#                 start = start,
#                 w.init = GLASSO$w, wi.init = GLASSO$wi,
#                 # w.init = init, wi.init = initOmega,
#                 trace = FALSE, ...)

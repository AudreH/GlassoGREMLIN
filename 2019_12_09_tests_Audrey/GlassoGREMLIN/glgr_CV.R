glgr_CV = function(X = NULL, S = NULL, nlam = 10, lam.min.ratio = 0.01, 
                   lam = NULL, diagonal = FALSE, path = FALSE, tol = 1e-04, 
                   maxit = 10000, adjmaxit = NULL, K = 5, 
                   crit.cv = c("loglik", "AIC", "BIC"),
                   start = c("warm", "cold"), 
                   cores = 1, 
                   trace = c("progress", "print", "none"),
                   V = NULL, # MODIF matrix of edge penalization. Multiplier of lambda.
                   ## ARGUMENTS GREMLIN
                   nb_feat, # nombre de features par blocks (decoupage des matrices en reseaux, important)
                   blocks_names = names(nb_feat), # nom des blocks
                   v_Kmin = rep(2, length(nb_feat)), 
                   v_Kmax =  v_Kmin+3, 
                   verbose = FALSE, # de base il vaut mieux que ce soit faux pour de la crossvalidation
                   nbcores_grem = 1, 
                   # maxiterVE = 100,
                   ...){
  
  
  GLASSO = CVgraphical.lasso(X = X, S = S, nlam = nlam, lam.min.ratio = lam.min.ratio, 
                             lam = lam, diagonal = diagonal, path = path, tol = tol, maxit = maxit, 
                             adjmaxit = adjmaxit, K = K, 
                             crit.cv = crit.cv, 
                             start = start, cores = cores, trace = trace, V = V)
  
  A_GLASSO <- ifelse(GLASSO$Omega!=0 & row(GLASSO$Omega)!=col(GLASSO$Omega),1,0)
  
  resGREMLIN =  gremlin_list(A_GLASSO,  nb_feat = nb_feat, 
                             blocks_names = blocks_names,
                             v_Kmin = v_Kmin, v_Kmax = v_Kmax, 
                             verbose = verbose, 
                             # maxiterVE = maxiterVE,
                             nbCores = nbcores_grem
  )
  
  GLASSO = CVgraphical.lasso(X = X, S = GLASSO$w, 
                             nlam = nlam, lam.min.ratio = lam.min.ratio, lam = lam, 
                             diagonal = diagonal, path = path,
                             tol = tol, maxit = maxit, 
                             adjmaxit = adjmaxit, K = K, 
                             crit.cv = crit.cv, 
                             start = start, cores = cores, trace = trace, 
                             V = 1/(resGREMLIN$matrix_reconstruct))
  
  A_GLASSO2 = ifelse(GLASSO$Omega!=0 & row(GLASSO$Omega)!=col(GLASSO$Omega),1,0)
  
  diff = sum(abs(A_GLASSO - A_GLASSO2))/length(A_GLASSO)
  cat("Diff : ", diff, "\n")
  if(diff>0.05){
    cat("Diff between adjacency matrices >5%, starting iterations \n")
    iter = 1
    A_GLASSO =  list()
    A_GLASSO[[1]] =  ifelse(GLASSO$Omega!=0 & row(GLASSO$Omega)!=col(GLASSO$Omega),1,0)
    
    while(diff>0.05 & iter <= 10){
      cat(" * iter ", iter, " - ")
      resGREMLIN =  gremlin_list(A_GLASSO[[iter]],  nb_feat = nb_feat, 
                                 blocks_names = blocks_names,
                                 v_Kmin = v_Kmin, v_Kmax = v_Kmax, 
                                 verbose = verbose, 
                                 # maxiterVE = maxiterVE,
                                 nbCores = nbcores_grem
      )
      
      GLASSO = CVgraphical.lasso(X = X, S = GLASSO$w, 
                                 nlam = nlam, lam.min.ratio = lam.min.ratio, lam = lam, 
                                 diagonal = diagonal, path = path, 
                                 tol = tol, maxit = maxit, 
                                 adjmaxit = adjmaxit, K = K, 
                                 crit.cv = crit.cv, 
                                 start = start, cores = cores, trace = trace, 
                                 V = 1/(resGREMLIN$matrix_reconstruct))
      
      A_GLASSO[[iter+1]] =  ifelse(GLASSO$Omega!=0 & row(GLASSO$Omega)!=col(GLASSO$Omega),1,0)
      diff = sum(abs((A_GLASSO[[iter+1]] - A_GLASSO[[iter]])))/length(A_GLASSO[[iter+1]])
      cat(" Diff = ", diff, "\n")
      
      iter = iter + 1
    }
  }
  
  return(list(resGREMLIN = resGREMLIN, 
         GLASSO = GLASSO))
}
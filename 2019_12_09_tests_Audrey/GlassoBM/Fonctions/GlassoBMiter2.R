# ---- glBM_CV : ----

glBM_CV = function(X = NULL, S = NULL, nlam = 10, lam.min.ratio = 0.01, 
                   lam = NULL, diagonal = FALSE, path = FALSE, tol = 1e-04, 
                   maxit = 10000, adjmaxit = NULL, K = 5, 
                   solver = "shooting",
                   crit.cv = c("loglik", "AIC", "BIC"),
                   start = c("warm", "cold"), 
                   cores = 1, 
                   trace = c("progress", "print", "none"),
                   V = NULL, # MODIF matrix of edge penalization. Multiplier of lambda.
                   ## ARGUMENTS GREMLIN
                   membership_type,
                   verbosity=0,
                   autosave='',
                   plotting=character(0),
                   exploration_factor=1.5,
                   explore_min=4,
                   explore_max=Inf,
                   ncores=detectCores(),
                   thre.iter = 10^-6, # seuil pour commencer les iterations... 
                   n.iter = 10 # pas les iterations pour le glasso, mais les iterations pour l'alternance Glasso BM
                   ){
  
  
  GLASSO = CVgraphical.lasso(X = X, S = S, nlam = nlam, lam.min.ratio = lam.min.ratio, 
                             lam = lam, diagonal = diagonal, path = path, tol = tol, maxit = maxit, 
                             adjmaxit = adjmaxit, K = K, 
                             crit.cv = crit.cv, 
                             start = start, cores = cores, trace = trace, V = V, solver = solver)
  
  old_sig = GLASSO$Sigma
  
  BM_step = BM_gaussian_step(
    verbosity = verbosity,
    membership_type = membership_type, 
    explore_min = explore_min, 
    explore_max = explore_max,
    exploration_factor = exploration_factor,
    autosave = autosave,
    plotting = plotting,
    adj = GLASSO$Sigma,
    ncores = ncores)
  
  GLASSO = CVgraphical.lasso(X = X, S = GLASSO$Theta, nlam = nlam, lam.min.ratio = lam.min.ratio, 
                    lam = lam, diagonal = diagonal, path = path, tol = tol, maxit = maxit, 
                    adjmaxit = adjmaxit, K = K, 
                    crit.cv = crit.cv, 
                    start = start, cores = cores, trace = trace, V = BM_step$mat_penalty, solver = solver)
  
  first_estimate = GLASSO
  first_BM = BM_step
  
  diff = sum(abs(old_sig-GLASSO$Sigma))
  print(diff)
  diff_vect = c(diff)
  STOP = FALSE
  iter = 1
  
  # if(diff>thre.iter & diff<10*length(GLASSO$Sigma)) cat("Starting iterations because of the difference\n")
  # else if(diff>10*length(GLASSO$Sigma)){
  #   STOP = TRUE
  #   cat("It seems the first estimate is too high to even start, stoping for safety...\n")
  # }
  
  
  while(diff_vect[length(diff_vect)]>thre.iter & iter < n.iter & STOP == FALSE){
    
    old_sig = GLASSO$Sigma
    
    BM_step = BM_gaussian_step(
      verbosity = verbosity,
      membership_type = membership_type, 
      explore_min = explore_min, 
      explore_max = explore_max,
      exploration_factor = exploration_factor,
      autosave = autosave,
      plotting = plotting,
      adj = GLASSO$Sigma, ncores = ncores)
    
    
    GLASSO = CVgraphical.lasso(X = X, S = GLASSO$Theta, nlam = nlam, lam.min.ratio = lam.min.ratio, 
                               lam = lam, diagonal = diagonal, path = path, tol = tol, maxit = maxit, 
                               adjmaxit = adjmaxit, K = K, 
                               crit.cv = crit.cv, 
                               start = start, cores = cores, trace = trace, V = BM_step$mat_penalty, solver = solver)
    
    diff_vect = c(diff_vect, sum(abs(old_sig-GLASSO$Sigma)))
    print(diff_vect[length(diff_vect)])
    
    # if(diff_vect[length(diff_vect)]>10*length(GLASSO$Sigma)){
    #   STOP = TRUE
    #   cat("It seems the algorithm is not converging, stoping for safety... Keeping the first estimate as output \n")
    #   GLASSO = first_estimate
    #   BM_step = first_BM
    # }
    
    iter = iter+1
  }
  
  GLASSO$wi = GLASSO$Theta
  GLASSO$w = GLASSO$Sigma
  
  return(list(BM_step = BM_step, 
              GLASSO = GLASSO))
}
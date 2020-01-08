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
                   ## ARGUMENTS Block Model
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
  
  listGLASSO = NULL
  GLASSO = CVgraphical.lasso(X = X, S = S, nlam = nlam, lam.min.ratio = lam.min.ratio, 
                             lam = lam, diagonal = diagonal, path = path, tol = tol, maxit = maxit, 
                             adjmaxit = adjmaxit, K = K, 
                             crit.cv = crit.cv, 
                             start = start, cores = cores, trace = trace, V = V, solver = solver)
  listGLASSO[[length(listGLASSO)+1]] = GLASSO
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
  listGLASSO[[length(listGLASSO)+1]] = GLASSO
  first_estimate = GLASSO
  first_BM = BM_step
  
  diff = sum(abs(old_sig-GLASSO$Sigma))
  print(diff)
  diff_vect = c(diff)
  STOP = FALSE
  iter = 1
  
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
    listGLASSO[[length(listGLASSO)+1]] = GLASSO
    diff_vect = c(diff_vect, sum(abs(old_sig-GLASSO$Sigma)))
    print(diff_vect[length(diff_vect)])
    
    iter = iter+1
  }
  
  GLASSO$wi = GLASSO$Theta
  GLASSO$w = GLASSO$Sigma
  mat_penalty = BM_step$mat_penalty
  if(!diagonal) diag(mat_penalty) = 0
  
  return(list(BM_step = BM_step, 
              GLASSO = GLASSO, 
              mat_penalty = mat_penalty,
              listGLASSO = listGLASSO))
}
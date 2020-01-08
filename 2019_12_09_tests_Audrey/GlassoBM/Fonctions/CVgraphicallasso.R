require(glasso)
require(parallel)
require(dplyr)
require(doParallel)

## Matt Galloway
# Modifications 2/12/19 : ajout d'un parametre V, matrice de proportions, a multiplier par rho, pour que la regularisation soit
# differente d'une arete a l'autre. V sera la matrice des probas de connexions (inverse terme a terme pour que les aretes ayant
# une plus forte proba de connexion soient moins regularisee ?)
# Le reste de la fonction n'a pas ete change (pour l'instant)
# Modifications 6/12/19 et 9/12/19 : changement pour graphical.lasso (simone) au lieu de glasso (fortran, qui bug)

###############################################################################################################
# ----- CVgraphical.lasso FUNC : -----

CVgraphical.lasso = function(X = NULL, 
                             S = NULL, # Can be an initial guess of S.
                             nlam = 10, # Longueur de sequence de lambdas que l'on veut tester
                             lam.min.ratio = 0.01, # donne le lambda de debut de sequence en fonction du lambda max
                             lam = NULL, # Sequence de lambdas a tester
                             diagonal = FALSE,  # Penalisation de la diagonale ou non
                             path = FALSE, 
                             tol = 1e-04, 
                             maxit = 10000, # Maximum iteration pour le glasso.
                             adjmaxit = NULL,
                             start = c("warm", "cold"), # use initial guess for starting the graphical lasso
                             K = 5, # Nombre de folds pour la validation croisee
                             crit.cv = c("loglik",   "AIC", "BIC"), # Critere pour la validation croisee
                             cores = 1, # Pour le parallele ou  non
                             trace = c("progress", "print", "none"), # verbose, du plus causant au silencieux.
                             V = NULL, # MODIF matrix of edge penalization. Multiplier of lambda.
                             solver = "shooting",
                             ...) {
  
  ### checks
  if (is.null(X) && is.null(S)) {
    stop("Must provide entry for X or S!")
  }
  if (!all(lam > 0)) {
    stop("lam must be positive!")
  }
  if (!(all(c(tol, maxit, adjmaxit, K, cores) > 0))) {
    stop("Entry must be positive!")
  }
  if (!(all(sapply(c(tol, maxit, adjmaxit, K, cores, nlam, 
                     lam.min.ratio), length) <= 1))) {
    stop("Entry must be single value!")
  }
  if (all(c(maxit, adjmaxit, K, cores)%%1 != 0)) {
    stop("Entry must be an integer!")
  }
  if (cores < 1) {
    stop("Number of cores must be positive!")
  }
  if (cores > 1 && path) {
    cat("Parallelization not possible when producing solution path. Setting cores = 1...\n\n")
    cores = 1
  }
  K = ifelse(path, 1, K)
  if (cores > K) {
    cat("Number of cores exceeds K... setting cores = K\n\n")
    cores = K
  }
  if (is.null(adjmaxit)) {
    adjmaxit = maxit
  }
  
  # match values
  crit.cv = match.arg(crit.cv)
  start = match.arg(start)
  trace = match.arg(trace)
  call = match.call()
  MIN.error = AVG.error = CV.error = NULL
  n = ifelse(is.null(X), nrow(S), nrow(X))
  
  ###  compute sample covariance matrix, if necessary
  if (is.null(S)) {
    S = (nrow(X) - 1)/nrow(X) * cov(X)
  }
  if(!is.null(V)){if(ncol(V)!=ncol(S) | nrow(V)!=nrow(S)) stop("V must have the same dimensions as S")} ### MODIF
  
  Sminus = S
  diag(Sminus) = 0
  
  ### compute grid of lam values, if necessary
  if (is.null(lam)) {
    if (!((lam.min.ratio <= 1) && (lam.min.ratio > 0))) {
      cat("lam.min.ratio must be in (0, 1]... setting to 1e-2!\n\n")
      lam.min.ratio = 0.01
    }
    if (!((nlam > 0) && (nlam%%1 == 0))) {
      cat("nlam must be a positive integer... setting to 10!\n\n")
      nlam = 10
    }
    
    # calculate lam.max and lam.min
    lam.max = max(abs(Sminus))
    lam.min = lam.min.ratio * lam.max
    
    # calculate grid of lambda values
    lam = 10^seq(log10(lam.min), log10(lam.max), length = nlam)
    
  } else {
    # sort lambda values
    lam = sort(lam)
    
  }
  
  ### perform cross validation, if necessary
  if ((length(lam) > 1) & (!is.null(X) || path)) {
    
    # run CV in parallel?
    if (cores > 1) {
      
      # execute CVP
      GLASSO = CVPgl(X = X, lam = lam, diagonal = diagonal, 
                     tol = tol, maxit = maxit, adjmaxit = adjmaxit, 
                     K = K, crit.cv = crit.cv, cores = cores, 
                     trace = trace, V = V, solver = solver, ...) 
      MIN.error = GLASSO$min.error
      AVG.error = GLASSO$avg.error
      CV.error = GLASSO$cv.error
      
    } else {
      
      # execute CV_ADMMc
      if (is.null(X)) {
        X = matrix(0)
      }
      GLASSO = CVgl(X = X, lam = lam, diagonal = diagonal, 
                    path = path, tol = tol, maxit = maxit, adjmaxit = adjmaxit, 
                    K = K, crit.cv = crit.cv, trace = trace, 
                    V = V, solver = solver, ...)
      MIN.error = GLASSO$min.error
      AVG.error = GLASSO$avg.error
      CV.error = GLASSO$cv.error
      Path = GLASSO$path
      
    }
    
    # print warning if lam on boundary
    if ((GLASSO$lam == lam[1]) && (length(lam) != 1) && 
        !path) {
      cat("\nOptimal tuning parameter on boundary... consider providing a smaller lam value or decreasing lam.min.ratio!")
    }
    
    # specify initial estimate for Sigma
    if (diagonal) {
      
      # simply force init to be positive definite final diagonal
      # elements will be increased by lam
      init = S + GLASSO$lam
      
    } else {
      
      # provide estimate that is pd and dual feasible
      alpha = min(c(GLASSO$lam/max(abs(Sminus)), 1))
      init = (1 - alpha) * S
      diag(init) = diag(S)
      
    }
    
    # compute final estimate at best tuning parameters
    lam_ = GLASSO$lam
    
    initial.guess = NULL
    if(start == "warm"){
      initial.guess = diag(ncol(S))
    }
    
    Rho = NULL
    if(is.null(V)){
      Rho = matrix(lam_, ncol(X), ncol(X))
    }else{
      Rho = lam_*V
    }
    if(!diagonal) diag(Rho) = 0
    GLASSO = graphical.lasso(X = X, 
                             Rho = Rho, 
                             initial.guess = initial.guess,
                             eps = 1e-6, maxIt = maxit, solver = solver) 
    
    GLASSO$wi = GLASSO$Theta
    GLASSO$w = GLASSO$Sigma
    GLASSO$lam = lam_
    GLASSO$Rho = Rho
    
  } else { # No cross_validation case
    
    # execute ADMM_sigmac
    if (length(lam) > 1) {
      stop("Must set specify X, set path = TRUE, or provide single value for lam.")
    }
    
    # specify initial estimate for Sigma
    if (diagonal) {
      
      # simply force init to be positive definite final diagonal
      # elements will be increased by lam
      init = S + lam
      
    } else {
      
      # provide estimate that is pd and dual feasible
      alpha = min(c(lam/max(abs(Sminus)), 1))
      init = (1 - alpha) * S
      diag(init) = diag(S)
      
    }
    
    initial.guess = NULL
    Sigma0 = NULL
    if(start == "warm"){
      initial.guess = diag(ncol(S))
      Sigma0 = init
    }
    
    Rho = NULL
    if(is.null(V)){
      Rho = matrix(lam, ncol(X), ncol(X))
    }else{
      Rho = lam*V
    }
    if(!diagonal) diag(Rho) = 0
    GLASSO = graphical.lasso(X = X, Rho = Rho, eps = 1e-6, 
                             maxIt = maxit, solver = solver, initial.guess = initial.guess) 
    
    GLASSO$wi = GLASSO$Theta
    GLASSO$w = GLASSO$Sigma
    GLASSO$lam = lam
    GLASSO$Rho = Rho
  }
  
  
  # option to penalize diagonal
  if (diagonal) { C = 1
  } else {
    C = 1 - diag(ncol(S))
  }
  
  # compute penalized loglik
  # if(is.null(V)){ ### MODIF 
  #   loglik = (-n/2) * (sum(GLASSO$wi * S) - determinant(GLASSO$wi, 
  #                                                       logarithm = TRUE)$modulus[1] + GLASSO$lam * sum(abs(C * GLASSO$wi)))
  # }else{
  #   loglik = (-n/2) * (sum(GLASSO$wi * S) - determinant(GLASSO$wi, 
  #                                                       logarithm = TRUE)$modulus[1] +  GLASSO$lam * sum(V * abs(C * GLASSO$wi)))
  # }
  
    loglik = (-n/2) * (sum(GLASSO$wi * S) - determinant(GLASSO$wi,
                                                        logarithm = TRUE)$modulus[1] + sum(abs(GLASSO$Rho * GLASSO$wi)))
  
  
  # return values
  tuning = matrix(c(log10(GLASSO$lam), GLASSO$lam), ncol = 2)
  colnames(tuning) = c("log10(lam)", "lam")
  if (!path) {
    Path = NULL
  }
  
  returns = list(Call = call, Iterations = GLASSO$niter, 
                 Tuning = tuning, Lambdas = lam, maxit = maxit, Omega = GLASSO$wi, 
                 Sigma = GLASSO$w, Path = Path, Loglik = loglik, MIN.error = MIN.error, 
                 AVG.error = AVG.error, CV.error = CV.error, GLASSO = GLASSO)
  
  class(returns) = "CVglasso"
  return(returns)
  
}

#############################################################################################################"

# ----- CVgl FUNC : -----
#' @title Cross Validation
#' 
#' 

CVgl = function(X = NULL, S = NULL, lam = 10^seq(-2, 2, 0.2), 
                diagonal = FALSE, path = FALSE, tol = 1e-04, maxit = 10000, 
                adjmaxit = NULL, K = 5, crit.cv = c("loglik", "AIC", "BIC"), 
                start = c("warm", "cold"), # use initial guess for starting graphical lasso
                cores = 1, trace = c("progress", "print", "none"),
                solver = "shooting",
                V = NULL, ### MODIF
                ...) {
  
  ### match values
  crit.cv = match.arg(crit.cv)
  start = match.arg(start)
  trace = match.arg(trace)
  lam = sort(lam)
  
  ### initialize
  Path = NULL
  initmaxit = maxit
  # S.train = S.valid = S
  CV_errors = array(0, c(length(lam), K))
  
  ### set progress bar
  if (trace == "progress") {
    progress = txtProgressBar(max = K * length(lam), style = 3)
  }
  
  ### no need to create folds if K = 1
  if (K == 1) {
    
    # set sample size
    n = nrow(S)
    
    # initialize Path, if necessary
    if (path) {
      Path = array(0, c(ncol(S), ncol(S), length(lam)))
    }
    
  } else {
    
    # designate folds and shuffle -- ensures randomized folds
    n = nrow(X)
    ind = sample(n)
    
  }
  
  ### parse data into folds and perform CV
  for (k in 1:K) {
    if (K > 1) {
      
      # cat(paste0("Fold ", k, "\n"))
      
      # training set
      leave.out = ind[(1 + floor((k - 1) * n/K)):floor(k * n/K)]
      X.train = X[-leave.out, , drop = FALSE]
      X_bar = apply(X.train, 2, mean)
      X.train = scale(X.train, center = X_bar, scale = FALSE)
      
      # validation set
      X.valid = X[leave.out, , drop = FALSE]
      X.valid = scale(X.valid, center = X_bar, scale = FALSE)
      n = nrow(X.valid)
      
      # sample covariances
      S.train = crossprod(X.train)/(dim(X.train)[1])
      S.valid = crossprod(X.valid)/(dim(X.valid)[1])
      
    }
    
    # re-initialize values for each fold
    maxit = initmaxit
    init = S.train
    initOmega = diag(ncol(S.train))
    
    ### loop over all tuning parameters
    for (i in 1:length(lam)) {
      
      # set temporary tuning parameter
      lam_ = lam[i]
      
      # print(lam_)
      
      # compute the penalized likelihood precision matrix estimator
      
      initial.guess = NULL
      if(start == "warm"){initial.guess = initOmega}
      
      Rho = NULL
      if(is.null(V)){
        Rho = matrix(lam_, ncol(X), ncol(X))
      }else{
        Rho = lam_*V
      }
      
      if(!diagonal) diag(Rho) = 0
      
      GLASSO = graphical.lasso(X = X, 
                               Rho = Rho, 
                               initial.guess = initial.guess,
                               eps = 1e-6, maxIt = maxit,
                               solver = solver) 
      
      GLASSO$wi = GLASSO$Theta
      GLASSO$w = GLASSO$Sigma
      GLASSO$lam = lam_
      
      # if (start == "warm") {
      #   
      #   # option to save initial values for warm starts
      #   init = GLASSO$w
      #   initOmega = GLASSO$wi
      #   maxit = adjmaxit
      #   
      # }
      
      # compute the observed negative validation loglikelihood
      # (close enoug)
      CV_errors[i, k] = (nrow(X)/2) * (sum(GLASSO$wi * S.valid) - determinant(GLASSO$wi, logarithm = TRUE)$modulus[1])
      
      # update for crit.cv, if necessary
      if (crit.cv == "AIC") {
        CV_errors[i, k] = CV_errors[i, k] + sum(GLASSO$wi !=   0)
      }
      if (crit.cv == "BIC") {
        CV_errors[i, k] = CV_errors[i, k] + sum(GLASSO$wi !=   0) * log(nrow(X))/2
      }
      
      # save estimate if path = TRUE
      if (path) {
        Path[, , i] = GLASSO$wi
      }
      
      # update progress bar
      if (trace == "progress") {
        setTxtProgressBar(progress, i + (k - 1) * length(lam))
        
        # if not quiet, then print progress lambda
      } else if (trace == "print") {
        cat("\nFinished lam = ", paste(lam_, sep = ""), "\n")
      }
    }
    
    # if not quiet, then print progress kfold
    if (trace == "print") {
      cat("\nFinished fold ", paste(k, sep = ""), "\n")
    }
  }
  
  # determine optimal tuning parameters
  AVG = apply(CV_errors, 1, mean)
  best_lam = lam[which.min(AVG)]
  error = min(AVG)
  
  
  # return best lam and alpha values
  return(list(lam = best_lam, path = Path, min.error = error, 
              avg.error = AVG, cv.error = CV_errors))
  
}

#############################################################################################################"

# ----- CVgl FUNC : -----
CVPgl = function(X = NULL,
                 # S = NULL, 
                 lam = 10^seq(-2, 2, 0.2), diagonal = FALSE, 
                 tol = 1e-04, maxit = 10000, adjmaxit = NULL, K = 5,
                 crit.cv = c("loglik", "AIC", "BIC"), 
                 start = c("warm", "cold"), cores = 1, 
                 trace = c("progress", "print", "none"), 
                 V = NULL, ### MODIF
                 solver = "shooting",
                 ...) {
  
  # match values
  crit.cv = match.arg(crit.cv)
  start = match.arg(start)
  trace = match.arg(trace)
  lam = sort(lam)
  
  # make cluster and register cluster
  num_cores = detectCores()
  if (cores > num_cores) {
    cat("\nOnly detected", paste(num_cores, "cores...", 
                                 sep = " "))
  }
  if (cores > K) {
    cat("\nNumber of cores exceeds K... setting cores = K")
    cores = K
  }
  
  cluster = makeCluster(cores)
  clusterExport(cl=cluster, varlist=c("graphical.lasso"))
  registerDoParallel(cluster)
  
  # use cluster for each fold in CV
  n = nrow(X)
  ind = sample(n)
  k = NULL
  CV = foreach(k = 1:K, 
               .packages = c("CVglasso", "glasso","igraph", "network", "simone"),
               .combine = "cbind", 
               .inorder = FALSE) %dopar% {
                 
                 # set progress bar
                 if (trace == "progress") {
                   progress = txtProgressBar(max = length(lam), style = 3)
                 }
                 
                 # training set
                 leave.out = ind[(1 + floor((k - 1) * n/K)):floor(k *  n/K)]
                 X.train = X[-leave.out, , drop = FALSE]
                 X_bar = apply(X.train, 2, mean)
                 X.train = scale(X.train, center = X_bar, scale = FALSE)
                 
                 # validation set
                 X.valid = X[leave.out, , drop = FALSE]
                 X.valid = scale(X.valid, center = X_bar, scale = FALSE)
                 
                 # sample covariances
                 S.train = crossprod(X.train)/(dim(X.train)[1])
                 S.valid = crossprod(X.valid)/(dim(X.valid)[1])
                 
                 # initial estimates
                 # init = S.train
                 initOmega = diag(ncol(S.train))
                 CV_error = array(0, length(lam))
                 
                 # # initial sigma
                 # if (!diagonal) {
                 #   
                 #   # provide estimate that is pd and dual feasible
                 #   Sminus = S.train
                 #   diag(Sminus) = 0
                 #   alpha = min(c(lam[1]/max(abs(Sminus)), 1))
                 #   init = (1 - alpha) * S.train
                 #   diag(init) = diag(S.train)
                 #   
                 # }
                 
                 # loop over all tuning parameters
                 for (i in 1:length(lam)) {
                   
                   # set temporary tuning parameter
                   lam_ = lam[i]
                   
                   # # update diagonal elements of init, if necessary
                   # if (diagonal) {
                   #   diag(init) = diag(S.train) + lam_
                   # }
                   
                   # compute the penalized likelihood precision matrix
                   # estimator
                   initial.guess = NULL
                   if(start == "warm"){initial.guess = initOmega}
                   
                   Rho = NULL
                   if(is.null(V)){
                     Rho = matrix(lam_, ncol(X), ncol(X))
                   }else{
                     Rho = lam_*V
                   }
                   
                   if(!diagonal) diag(Rho) = 0
                   GLASSO = graphical.lasso(X = X, 
                                            Rho = Rho, 
                                            initial.guess = initial.guess, solver = solver, 
                                            eps = 1e-6, maxIt = maxit) 
                   
                   GLASSO$wi = GLASSO$Theta
                   GLASSO$w = GLASSO$Sigma
                   GLASSO$lam = lam_
                   
                   
                   # if (start == "warm") {
                   #   
                   #   # option to save initial values for warm starts
                   #   init = GLASSO$w
                   #   initOmega = GLASSO$wi
                   #   maxit = adjmaxit
                   #   
                   # }
                   
                   # compute the observed negative validation loglikelihood
                   # (close enoug)
                   CV_error[i] = (nrow(X)/2) * (sum(GLASSO$wi * S.valid) - 
                                                  determinant(GLASSO$wi, logarithm = TRUE)$modulus[1])
                   
                   # update for crit.cv, if necessary
                   if (crit.cv == "AIC") {
                     CV_error[i] = CV_error[i] + sum(GLASSO$wi != 0)
                   }
                   if (crit.cv == "BIC") {
                     CV_error[i] = CV_error[i] + sum(GLASSO$wi != 0) * log(nrow(X))/2
                   }
                   
                   # update progress bar
                   if (trace == "progress") {
                     setTxtProgressBar(progress, i)
                     
                     # if not quiet, then print progress lambda
                   } else if (trace == "print") {
                     cat("\nFinished lam = ", paste(lam_, sep = ""))
                   }
                 }
                 
                 # return CV errors
                 return(CV_error)
                 
               }
  
  # determine optimal tuning parameters
  AVG = apply(CV, 1, mean)
  best_lam = lam[which.min(AVG)]
  error = min(AVG)
  
  # stop cluster
  stopCluster(cluster)
  
  # return best lam and alpha values
  return(list(lam = best_lam, min.error = error, avg.error = AVG, 
              cv.error = CV))
  
}

##-----------------------------------------------------------------------------------
# ----- PRINT FUNC : -----

print.CVglasso = function(x, ...) {
  
  # print warning if maxit reached
  if (x$maxit <= x$Iterations) {
    cat("\nMaximum iterations reached...!")
  }
  
  # print call
  cat("\n\nCall: ", paste(deparse(x$Call), sep = "\n", collapse = "\n"), 
      "\n", sep = "")
  
  # print iterations
  cat("\nIterations:\n")
  print.default(x$Iterations, quote = FALSE)
  
  # print optimal tuning parameters
  cat("\nTuning parameter:\n")
  print.default(round(x$Tuning, 3), print.gap = 2L, quote = FALSE)
  
  # print loglik
  cat("\nLog-likelihood: ", paste(round(x$Loglik, 5), sep = "\n", 
                                  collapse = "\n"), "\n", sep = "")
  
  # print Omega if dim <= 10
  if (nrow(x$Omega) <= 10) {
    cat("\nOmega:\n")
    print.default(round(x$Omega, 5))
  } else {
    cat("\n(...output suppressed due to large dimension!)\n")
  }
  
}

# ----- PLOT FUNC : -----
plot.CVglasso = function(x, type = c("line", "heatmap"), footnote = TRUE, 
                         main = "Cross-Validation Errors",
                         ...) {
  
  # check
  type = match.arg(type)
  Means = NULL
  if (is.null(x$CV.error)) {
    stop("No cross validation errors to plot!")
  }
  
  if (type == "line") {
    
    # gather values to plot
    cv = cbind(expand.grid(lam = x$Lambdas, alpha = 0), 
               Errors = as.data.frame.table(x$CV.error)$Freq)
    
    # produce line graph
    graph = ggplot(summarise(group_by(cv, lam), Means = mean(Errors)), 
                   aes(log10(lam), Means)) + 
      geom_jitter(width = 0.2, color = "navy blue") + theme_minimal() + geom_line(color = "red") + 
      labs(title = main, y = "Error") + 
      geom_vline(xintercept = x$Tuning[1], linetype = "dotted")
    
  } else {
    
    # augment values for heat map (helps visually)
    lam = x$Lambdas
    cv = expand.grid(lam = lam, alpha = 0)
    Errors = 1/(c(x$AVG.error) + abs(min(x$AVG.error)) + 
                  1)
    cv = cbind(cv, Errors)
    
    # design color palette
    bluetowhite <- c("#000E29", "white")
    
    # produce ggplot heat map
    graph = ggplot(cv, aes(alpha, log10(lam))) + geom_raster(aes(fill = Errors)) + 
      scale_fill_gradientn(colours = colorRampPalette(bluetowhite)(2), 
                           guide = "none") + theme_minimal() + labs(title = "Heatmap of Cross-Validation Errors") + 
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
            axis.ticks.x = element_blank())
    
  }
  
  if (footnote) {
    
    # produce with footnote
    graph + labs(caption = paste("**Optimal: log10(lam) = ", 
                                 round(x$Tuning[1], 3), sep = ""))
    
  } else {
    
    # produce without footnote
    graph
    
  }
}


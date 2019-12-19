require(R.utils)

require(parallel)
require(doParallel)
require(dplyr)

## Matt Galloway
# Modifications 2/12/19 : ajout d'un parametre V, matrice de proportions, a multiplier par rho, pour que la regularisation soit
# differente d'une arete a l'autre. V sera la matrice des probas de connexions (inverse terme a terme pour que les aretes ayant
# une plus forte proba de connexion soient moins regularisee ?)
# Le reste de la fonction n'a pas ete change (pour l'instant)
# Modifications 10/12/19 : on change Gremlin pour BM simple. 


###############################################################################################################
# ----- CVglasso FUNC : -----

CVglassoBM = function(X = NULL, S = NULL, nlam = 10, lam.min.ratio = 0.01, 
                      lam = NULL, diagonal = FALSE, path = FALSE, tol = 1e-04, 
                      maxit = 10000, adjmaxit = NULL, K = 5, 
                      crit.cv = c("loglik", "AIC", "BIC"),
                      # start = c("warm", "cold"), 
                      cores = 1, 
                      trace = c("progress", "print", "none"),
                      solver = "shooting",
                      ## ARGUMENTS BM
                      membership_type = 'SBM', 
                      verbosity=0,
                      autosave='',
                      plotting=character(0),
                      exploration_factor=1.5,
                      explore_min=4,
                      explore_max=Inf,
                      ncores=detectCores(),
                      ## ARGUMENTS Alernance BM/Glasso
                      thre.iter = 10^-6, # seuil pour commencer les iterations... 
                      n.iter = 10, # pas les iterations pour le glasso, mais les iterations pour l'alternance Glasso BM
                      ...) {
  
  # checks
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
  # start = match.arg(start)
  trace = match.arg(trace)
  call = match.call()
  MIN.error = AVG.error = CV.error = NULL
  n = ifelse(is.null(X), nrow(S), nrow(X))
  
  # compute sample covariance matrix, if necessary
  if (is.null(S)) {
    S = (nrow(X) - 1)/nrow(X) * cov(X)
  }
  
  Sminus = S
  diag(Sminus) = 0
  
  # compute grid of lam values, if necessary
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
  
  # perform cross validation, if necessary
  if ((length(lam) > 1) & (!is.null(X) || path)) {
    
    # run CV in parallel?
    if (cores > 1) {
      
      # execute CVP ##### A CHANGER
      GLASSO = CVP_glBM(X = X, S = S, lam = lam, diagonal = diagonal, 
                        path = path, tol = tol, adjmaxit = adjmaxit, 
                        cores = cores,
                        K = K, crit.cv = crit.cv,  
                        maxit = maxit,
                        trace = trace,
                        solver = solver,
                        membership_type = membership_type, 
                        verbosity=verbosity,
                        autosave=autosave,
                        plotting=plotting,
                        exploration_factor=exploration_factor,
                        explore_min=explore_min,
                        explore_max=explore_max,
                        ncores=ncores,
                        thre.iter = thre.iter, # seuil pour commencer les iterations... 
                        n.iter = n.iter, # pas les iterations pour le glasso, mais les iterations pour l'alternance Glasso BM
                        ...) ## MODIF
      MIN.error = GLASSO$min.error
      AVG.error = GLASSO$avg.error
      CV.error = GLASSO$cv.error
      
    } else {
      
      # execute CV_ADMMc
      if (is.null(X)) {
        X = matrix(0)
      }
      # cat("I'm starting the cross-validation \n")
      GLASSO = CV_glBM(X = X, S = S, lam = lam, diagonal = diagonal, 
                       path = path, tol = tol, adjmaxit = adjmaxit, 
                       K = K, crit.cv = crit.cv,  
                       maxit = maxit,
                       trace = trace,
                       solver = solver,
                       membership_type = membership_type, 
                       verbosity=verbosity,
                       autosave=autosave,
                       plotting=plotting,
                       exploration_factor=exploration_factor,
                       explore_min=explore_min,
                       explore_max=explore_max,
                       ncores=ncores,
                       thre.iter = thre.iter, # seuil pour commencer les iterations... 
                       n.iter = n.iter, # pas les iterations pour le glasso, mais les iterations pour l'alternance Glasso BM
                       ...) ## MODIF
      MIN.error = GLASSO$min.error
      AVG.error = GLASSO$avg.error
      CV.error = GLASSO$cv.error
      Path = GLASSO$path
      
      # cat("\n I finished the CV glasso ! \n")
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
    
    res = BMGlassoIter(X = X, rho = lam_,
                       thr = tol,
                       maxit = maxit,
                       penalize.diagonal = diagonal,
                       trace = trace,
                       solver = solver,
                       membership_type = membership_type, 
                       verbosity=verbosity,
                       autosave=autosave,
                       plotting=plotting,
                       exploration_factor=exploration_factor,
                       explore_min=explore_min,
                       explore_max=explore_max,
                       ncores=ncores,
                       thre.iter = thre.iter, # seuil pour commencer les iterations... 
                       n.iter = n.iter # pas les iterations pour le glasso, mais les iterations pour l'alternance Glasso BM
    )
    
    
    GLASSO = res$GLASSO
    GLASSO$lam = lam_
    V = res$BM_step$mat_penalty
    
  } else {
    
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
    
    res = BMGlassoIter(X = X, rho = lam_,
                       thr = tol,
                       maxit = maxit,
                       penalize.diagonal = diagonal,
                       trace = trace,
                       solver = solver,
                       membership_type = membership_type, 
                       verbosity=verbosity,
                       autosave=autosave,
                       plotting=plotting,
                       exploration_factor=exploration_factor,
                       explore_min=explore_min,
                       explore_max=explore_max,
                       ncores=ncores,
                       thre.iter = thre.iter, # seuil pour commencer les iterations... 
                       n.iter = n.iter # pas les iterations pour le glasso, mais les iterations pour l'alternance Glasso BM
    )
    
    
    
    GLASSO = res$GLASSO
    
    GLASSO$lam = lam
    V = res$BM_step$mat_penalty
    
  }
  
  
  # option to penalize diagonal
  if (diagonal) {
    C = 1
  } else {
    C = 1 - diag(ncol(S))
  }
  
  # compute penalized loglik
  ### MODIF 
  loglik = (-n/2) * (sum(GLASSO$wi * S) - determinant(GLASSO$wi,
                                                      logarithm = TRUE)$modulus[1] +
                       GLASSO$lam * sum(V * abs(C * GLASSO$wi)))
  
  # return values
  tuning = matrix(c(log10(GLASSO$lam), GLASSO$lam), ncol = 2)
  colnames(tuning) = c("log10(lam)", "lam")
  if (!path) {
    Path = NULL
  }
  
  returns = list(Call = call,
                 res = res,
                 Iterations = GLASSO$niter,
                 Tuning = tuning, Lambdas = lam, 
                 maxit = maxit,
                 Omega = GLASSO$wi,
                 Sigma = GLASSO$w,
                 V = V,
                 Path = Path, Loglik = loglik, MIN.error = MIN.error, 
                 AVG.error = AVG.error, CV.error = CV.error)
  
  class(returns) = "CVglasso"
  return(returns)
  
}


#############################################################################################################"

# ----- CV FUNC : -----
#' @title Parallel Cross Validation
#' @description Parallel implementation of cross validation.
#'
#' @param X option to provide a nxp data matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
#' @param S option to provide a pxp sample covariance matrix (denominator n). If argument is \code{NULL} and \code{X} is provided instead then \code{S} will be computed automatically.
#' @param lam positive tuning parameters for elastic net penalty. If a vector of parameters is provided, they should be in increasing order. Defaults to grid of values \code{10^seq(-2, 2, 0.2)}.
#' @param diagonal option to penalize the diagonal elements of the estimated precision matrix (\eqn{\Omega}). Defaults to \code{FALSE}.
#' @param path option to return the regularization path. This option should be used with extreme care if the dimension is large. If set to TRUE, cores must be set to 1 and errors and optimal tuning parameters will based on the full sample. Defaults to FALSE.
#' @param tol convergence tolerance. Iterations will stop when the average absolute difference in parameter estimates in less than \code{tol} times multiple. Defaults to 1e-4.
#' @param maxit maximum number of iterations. Defaults to 1e4.
#' @param adjmaxit adjusted maximum number of iterations. During cross validation this option allows the user to adjust the maximum number of iterations after the first \code{lam} tuning parameter has converged. This option is intended to be paired with \code{warm} starts and allows for 'one-step' estimators. Defaults to NULL.
#' @param K specify the number of folds for cross validation.
#' @param crit.cv cross validation criterion (\code{loglik}, \code{AIC}, or \code{BIC}). Defaults to \code{loglik}.
#' @param start specify \code{warm} or \code{cold} start for cross validation. Default is \code{warm}.
#' @param cores option to run CV in parallel. Defaults to \code{cores = 1}.
#' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.
#' @param ... additional arguments to pass to \code{glasso}.
#' 
#' @return returns list of returns which includes:
#' \item{lam}{optimal tuning parameter.}
#' \item{min.error}{minimum average cross validation error (cv.crit) for optimal parameters.}
#' \item{avg.error}{average cross validation error (cv.crit) across all folds.}
#' \item{cv.error}{cross validation errors (cv.crit).}
#' 
#' @keywords internal

# we define the CV function
CV_glBM = function(X = NULL, S = NULL, lam = 10^seq(-2, 2, 0.2), 
                   diagonal = FALSE, path = FALSE, tol = 1e-04, maxit = 10000, 
                   adjmaxit = NULL, K = 5, crit.cv = c("loglik", "AIC", "BIC"), 
                   # start = c("warm", "cold"), 
                   cores = 1, trace = c("progress", "print", "none"),
                   solver = "shooting",
                   ## ARGUMENTS BM
                   membership_type,
                   verbosity=0,
                   autosave='',
                   plotting=character(0),
                   exploration_factor=1.5,
                   explore_min=4,
                   explore_max=Inf,
                   ncores=detectCores(),
                   ## ARGUMENTS Alernance BM/Glasso
                   thre.iter = 10^-6, # seuil pour commencer les iterations... 
                   n.iter = 10, # pas les iterations pour le glasso, mais les iterations pour l'alternance Glasso BM
                   ...) {
  
  # match values
  crit.cv = match.arg(crit.cv)
  # start = match.arg(start)
  trace = match.arg(trace)
  lam = sort(lam)
  
  # initialize
  Path = NULL
  initmaxit = maxit
  S.train = S.valid = S
  CV_errors = array(0, c(length(lam), K))
  
  # set progress bar
  if (trace == "progress") {
    progress = txtProgressBar(max = K * length(lam), style = 3)
  }
  
  # no need to create folds if K = 1
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
  
  # parse data into folds and perform CV
  for (k in 1:K) {
    cat(paste0("* Fold ", k, "\n")) ### TO DELETE
    if (K > 1) {
      
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
    
    # initial sigma
    if (!diagonal) {
      
      # provide estimate that is pd and dual feasible
      Sminus = S.train
      diag(Sminus) = 0
      alpha = min(c(lam[1]/max(abs(Sminus)), 1))
      init = (1 - alpha) * S.train
      diag(init) = diag(S.train)
      
    }
    
    
    # loop over all tuning parameters
    for (i in 1:length(lam)) {
      
      # set temporary tuning parameter
      lam_ = lam[i]
      
      print(lam_) ### TO DELETE
      
      # update diagonal elements of init, if necessary
      if (diagonal) {
        diag(init) = diag(S.train) + lam_
      }
      
      # compute the penalized likelihood precision matrix estimator
      ### MODIF
      res = BMGlassoIter(X = X.train, rho = lam_,
                         thr = tol,
                         maxit = maxit,
                         penalize.diagonal = diagonal,
                         trace = trace,
                         solver = solver,
                         membership_type = membership_type, 
                         verbosity=verbosity,
                         autosave=autosave,
                         plotting=plotting,
                         exploration_factor=exploration_factor,
                         explore_min=explore_min,
                         explore_max=explore_max,
                         ncores=ncores,
                         thre.iter = thre.iter, # seuil pour commencer les iterations... 
                         n.iter = n.iter # pas les iterations pour le glasso, mais les iterations pour l'alternance Glasso BM
      )
      
      GLASSO = res$GLASSO
      
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
        setTxtProgressBar(progress, i + (k - 1) * 
                            length(lam))
        
        # if not quiet, then print progress lambda
      } else if (trace == "print") {
        cat("\nFinished lam = ", paste(lam_, sep = ""))
      }
    }
    
    # if not quiet, then print progress kfold
    if (trace == "print") {
      cat("\nFinished fold ", paste(k, sep = ""))
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

# ----- CVP FUNC : ----- 
# La meme que CV FUNC mais en parallele

CVP_glBM = function(X = NULL, S = NULL, lam = 10^seq(-2, 2, 0.2), 
                    diagonal = FALSE, path = FALSE, tol = 1e-04, maxit = 10000, 
                    adjmaxit = NULL, K = 5, crit.cv = c("loglik", "AIC", "BIC"), 
                    # start = c("warm", "cold"), 
                    cores = 1,
                    trace = c("progress", "print", "none"),
                    solver = "shooting",
                    ## ARGUMENTS BM
                    membership_type,
                    verbosity=0,
                    autosave='',
                    plotting=character(0),
                    exploration_factor=1.5,
                    explore_min=4,
                    explore_max=Inf,
                    ncores=detectCores(),
                    ## ARGUMENTS Alernance BM/Glasso
                    thre.iter = 10^-6, # seuil pour commencer les iterations... 
                    n.iter = 10, # pas les iterations pour le glasso, mais les iterations pour l'alternance Glasso BM
                    ...) {
  
  # match values
  crit.cv = match.arg(crit.cv)
  # start = match.arg(start)
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
  clusterExport(cl=cluster, varlist=c("BMGlassoIter", "BM_gaussian_step", "graphical.lasso"))
  registerDoParallel(cluster)
  
  # use cluster for each fold in CV
  n = nrow(X)
  ind = sample(n)
  k = NULL
  CV = foreach(k = 1:K, 
               .packages = c("CVglasso", "glasso", "GREMLIN", "igraph", "network", "simone", "blockmodels"),
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
                 init = S.train
                 initOmega = diag(ncol(S.train))
                 CV_error = array(0, length(lam))
                 
                 # initial sigma
                 if (!diagonal) {
                   
                   # provide estimate that is pd and dual feasible
                   Sminus = S.train
                   diag(Sminus) = 0
                   alpha = min(c(lam[1]/max(abs(Sminus)), 1))
                   init = (1 - alpha) * S.train
                   diag(init) = diag(S.train)
                   
                 }
                 
                 # loop over all tuning parameters
                 for (i in 1:length(lam)) {
                   
                   # set temporary tuning parameter
                   lam_ = lam[i]
                   
                   # update diagonal elements of init, if necessary
                   if (diagonal) {
                     diag(init) = diag(S.train) + lam_
                   }
                   
                   # compute the penalized likelihood precision matrix
                   # estimator
                   
                   
                   res = BMGlassoIter(X = X.train, rho = lam_,
                                      thr = tol,
                                      maxit = maxit,
                                      penalize.diagonal = diagonal,
                                      trace = trace,
                                      solver = solver,
                                      membership_type = membership_type, 
                                      verbosity=verbosity,
                                      autosave=autosave,
                                      plotting=plotting,
                                      exploration_factor=exploration_factor,
                                      explore_min=explore_min,
                                      explore_max=explore_max,
                                      ncores=ncores,
                                      thre.iter = thre.iter, # seuil pour commencer les iterations... 
                                      n.iter = n.iter # pas les iterations pour le glasso, mais les iterations pour l'alternance Glasso BM
                   )
                   GLASSO = res$GLASSO
                   
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
                   aes(log10(lam), Means)) + geom_jitter(width = 0.2, 
                                                         color = "navy blue") + theme_minimal() + geom_line(color = "red") + 
      labs(title = "Cross-Validation Errors", y = "Error") + 
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





# ---------------- BMGlassoIter : ----

BMGlassoIter = function(X, rho,
                        thr,
                        maxit = 25,
                        penalize.diagonal,
                        trace = FALSE,
                        solver = "shooting",
                        membership_type, 
                        # adj, 
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
  
  Rho = matrix(rho, ncol(X), ncol(X))
  if(!penalize.diagonal) diag(Rho) = 0
  
  GLASSO = graphical.lasso(X = X, 
                           Rho = Rho, 
                           eps = 1e-6, 
                           maxIt = maxit, 
                           solver = solver) 
  old_sig = GLASSO$Sigma
  
    BM_step = BM_gaussian_step(
      verbosity = verbosity,
      membership_type = membership_type, 
      explore_min = explore_min, 
      explore_max = explore_max,
      exploration_factor = exploration_factor,
      autosave = autosave,
      plotting = plotting,
      adj = GLASSO$Sigma)
    
    Rho = rho*BM_step$mat_penalty
    if(!penalize.diagonal) diag(Rho) = 0
    
    GLASSO = graphical.lasso(X = X,
                             Rho = Rho, 
                             initial.guess = GLASSO$Theta,
                             eps = 1e-6, maxIt = maxit, solver= solver) 
    
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
    #   cat("It seems the first estimate is too high to even start, stopping for safety...\n")
    # }
    # 
    
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
        adj = GLASSO$Sigma)
      
      Rho = rho*BM_step$mat_penalty
      if(!penalize.diagonal) diag(Rho) = 0
      
      GLASSO = graphical.lasso(X = X,
                               Rho = Rho, 
                               initial.guess = GLASSO$Theta,
                               eps = 1e-6, maxIt = maxit, solver= solver) 
      
      diff_vect = c(diff_vect, sum(abs(old_sig-GLASSO$Sigma)))
      print(diff_vect[length(diff_vect)])
      
      # if(diff_vect[length(diff_vect)]>length(GLASSO$Sigma)){
      #   STOP = TRUE
      #   cat("It seems the algorithm is not converging, stoping for safety... Keeping the first estimate as output \n")
      #   GLASSO = first_estimate
      #   BM_step = first_BM
      # }
      
      iter = iter+1
    }
    
    GLASSO$wi = GLASSO$Theta
    GLASSO$w = GLASSO$Sigma
  
  return(list(GLASSO = GLASSO, BM_step = BM_step))
}
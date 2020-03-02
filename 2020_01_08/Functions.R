require(glassoFast)
require(dplyr)

# Reprise fonction codees par Matt Galloway dans CVGlasso, pour y mettre du glassoFast a la place du glasso "normal".
# avantages glassoFast : n'a pas le probleme de non-convergence du glasso.
# A changer : glassoFast n'a pas d'argument X, on ne peut y rentrer que la matrice de variance/covariance de sdonnées. 
# Pas genant (c'est ce que l'on cherche à estimer.)

#' @author Matt Galloway \email{gall0441@@umn.edu}

#---- CVglasso ----

#' @param X option to provide a nxp data matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
#' @param S option to provide a pxp sample covariance matrix (denominator n). If argument is \code{NULL} and \code{X} is provided instead then \code{S} will be computed automatically.
#' @param nlam number of \code{lam} tuning parameters for penalty term generated from \code{lam.min.ratio} and \code{lam.max} (automatically generated). Defaults to 10.
#' @param lam.min.ratio smallest \code{lam} value provided as a fraction of \code{lam.max}. The function will automatically generate \code{nlam} tuning parameters from \code{lam.min.ratio*lam.max} to \code{lam.max} in log10 scale. \code{lam.max} is calculated to be the smallest \code{lam} such that all off-diagonal entries in \code{Omega} are equal to zero (\code{alpha} = 1). Defaults to 1e-2.
#' @param lam option to provide positive tuning parameters for penalty term. This will cause \code{nlam} and \code{lam.min.ratio} to be disregarded. If a vector of parameters is provided, they should be in increasing order. Defaults to NULL.
#' @param diagonal option to penalize the diagonal elements of the estimated precision matrix (\eqn{\Omega}). Defaults to \code{FALSE}.
#' @param path option to return the regularization path. This option should be used with extreme care if the dimension is large. If set to TRUE, cores must be set to 1 and errors and optimal tuning parameters will based on the full sample. Defaults to FALSE.
#' @param tol convergence tolerance. Iterations will stop when the average absolute difference in parameter estimates in less than \code{tol} times multiple. Defaults to 1e-4.
#' @param maxIt maximum number of iterations. Defaults to 1e4.
#' @param adjmaxIt adjusted maximum number of iterations. During cross validation this option allows the user to adjust the maximum number of iterations after the first \code{lam} tuning parameter has converged. This option is intended to be paired with \code{warm} starts and allows for 'one-step' estimators. Defaults to NULL.
#' @param K specify the number of folds for cross validation.
#' @param crit.cv cross validation criterion (\code{loglik}, \code{AIC}, or \code{BIC}). Defaults to \code{loglik}.
#' @param start specify \code{warm} or \code{cold} start for cross validation. Default is \code{warm}.
#' @param cores option to run CV in parallel. Defaults to \code{cores = 1}.
#' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.
#' @param ... additional arguments to pass to \code{glasso}.
#' 
#' @author Matt Galloway \email{gall0441@@umn.edu}
#' 

CVglasso = function(X = NULL, S = NULL, nlam = 10, lam.min.ratio = 0.01, 
                    lam = NULL, diagonal = FALSE, path = FALSE, tol = 1e-04, 
                    maxIt = 10000, adjmaxIt = NULL, K = 5, 
                    crit.cv = c("loglik", "AIC", "BIC"), 
                    start = c("warm", "cold"), cores = 1, 
                    trace = c("progress", "print", "none"),
                    lam_prior = matrix(1, ncol(X), ncol(X)),
                    ...) {
  
  # checks
  if (is.null(X) && is.null(S)) {
    stop("Must provide entry for X or S!")
  }
  if (!all(lam > 0)) {
    stop("lam must be positive!")
  }
  if (!(all(c(tol, maxIt, adjmaxIt, K, cores) > 0))) {
    stop("Entry must be positive!")
  }
  if (!(all(sapply(c(tol, maxIt, adjmaxIt, K, cores, nlam, 
                     lam.min.ratio), length) <= 1))) {
    stop("Entry must be single value!")
  }
  if (all(c(maxIt, adjmaxIt, K, cores)%%1 != 0)) {
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
  if (is.null(adjmaxIt)) {
    adjmaxIt = maxIt
  }
  
  # match values
  crit.cv = match.arg(crit.cv)
  start = match.arg(start)
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
  
  ### MODIFICATION
  
  # perform cross validation, if necessary
  if ((length(lam) > 1) & (!is.null(X) || path)) {
    
    # run CV in parallel?
    if (cores > 1) {
      
      # execute CVP
      GLASSO = CVP(X = X, lam = lam, diagonal = diagonal, 
                   tol = tol, maxIt = maxIt, adjmaxIt = adjmaxIt, 
                   K = K, crit.cv = crit.cv, start = start, cores = cores, 
                   trace = trace, lam_prior = lam_prior, ...)
      MIN.error = GLASSO$min.error
      AVG.error = GLASSO$avg.error
      CV.error = GLASSO$cv.error
      
    } else {
      
      # execute CV_ADMMc
      if (is.null(X)) {
        X = matrix(0)
      }
      GLASSO = CV(X = X, S = S, lam = lam, diagonal = diagonal, 
                  path = path, tol = tol, maxIt = maxIt, adjmaxIt = adjmaxIt, 
                  K = K, crit.cv = crit.cv, start = start, trace = trace, lam_prior = lam_prior,
                  ...)
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
      # init = S + GLASSO$lam
      
      init = S
      diag(init) = diag(S) + GLASSO$lam

    } else {

      # provide estimate that is pd and dual feasible
      alpha = min(c(GLASSO$lam/max(abs(Sminus)), 1))
      init = (1 - alpha) * S
      diag(init) = diag(S)

    }
    
    # compute final estimate at best tuning parameters
    lam_ = GLASSO$lam
    # if(!diagonal) diag(lam_prior) = 0
    # GLASSO = glassoFast(S = S, rho = lam_*lam_prior, thr = tol, maxIt = maxIt, 
    #                     start = "warm", 
    #                     w.init = init, wi.init = diag(ncol(S)), trace = FALSE, 
    #                     ...)
    GLASSO = glassoFast(S = S, rho = lam_*lam_prior, thr = tol, maxIt = maxIt, trace = FALSE, 
                        ...)
    GLASSO$lam = lam_
    GLASSO$lam_prior = lam_prior
    
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
    
    # GLASSO = glassoFast(S = S, rho = lam*lam_prior, thr = tol, maxIt = maxIt, 
                        # start = "warm", 
                        # w.init = init, wi.init = diag(ncol(S)), trace = FALSE, 
                        # ...)
    GLASSO = glassoFast(S = S, rho = lam_*lam_prior, thr = tol, maxIt = maxIt, trace = FALSE, 
                        ...)
    GLASSO$lam = lam
    GLASSO$lam_prior = lam_prior
  }
  
  
  # option to penalize diagonal
  if (diagonal) {
    C = 1
  } else {
    C = 1 - diag(ncol(S))
  }
  
  # compute penalized loglik
  loglik = (-n/2) * (sum(GLASSO$wi * S) - determinant(GLASSO$wi, logarithm = TRUE)$modulus[1] + (GLASSO$lam) * sum(lam_prior * abs(C * GLASSO$wi)))
  
  
  # return values
  tuning = matrix(c(log10(GLASSO$lam), GLASSO$lam), ncol = 2)
  colnames(tuning) = c("log10(lam)", "lam")
  if (!path) {
    Path = NULL
  }
  
  returns = list(Call = call, Iterations = GLASSO$niter, lam_prior = lam_prior,
                 Tuning = tuning, Lambdas = lam, maxIt = maxIt, Omega = GLASSO$wi, 
                 Sigma = GLASSO$w, Path = Path, Loglik = loglik, MIN.error = MIN.error, 
                 AVG.error = AVG.error, CV.error = CV.error, GLASSO = GLASSO)
  
  class(returns) = "CVglasso"
  return(returns)
  
}

# ---- CV -----
#' @param X option to provide a nxp data matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
#' @param S option to provide a pxp sample covariance matrix (denominator n). If argument is \code{NULL} and \code{X} is provided instead then \code{S} will be computed automatically.
#' @param lam positive tuning parameters for elastic net penalty. If a vector of parameters is provided, they should be in increasing order. Defaults to grid of values \code{10^seq(-2, 2, 0.2)}.
#' @param diagonal option to penalize the diagonal elements of the estimated precision matrix (\eqn{\Omega}). Defaults to \code{FALSE}.
#' @param path option to return the regularization path. This option should be used with extreme care if the dimension is large. If set to TRUE, cores must be set to 1 and errors and optimal tuning parameters will based on the full sample. Defaults to FALSE.
#' @param tol convergence tolerance. Iterations will stop when the average absolute difference in parameter estimates in less than \code{tol} times multiple. Defaults to 1e-4.
#' @param maxIt maximum number of iterations. Defaults to 1e4.
#' @param adjmaxIt adjusted maximum number of iterations. During cross validation this option allows the user to adjust the maximum number of iterations after the first \code{lam} tuning parameter has converged. This option is intended to be paired with \code{warm} starts and allows for 'one-step' estimators. Defaults to NULL.
#' @param K specify the number of folds for cross validation.
#' @param crit.cv cross validation criterion (\code{loglik}, \code{AIC}, or \code{BIC}). Defaults to \code{loglik}.
#' @param start specify \code{warm} or \code{cold} start for cross validation. Default is \code{warm}.
#' @param cores option to run CV in parallel. Defaults to \code{cores = 1}.
#' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.
#' @param ... additional arguments to pass to \code{glasso}.


# we define the CV function
CV = function(X = NULL, S = NULL, lam = 10^seq(-2, 2, 0.2), 
              diagonal = FALSE, path = FALSE, tol = 1e-04, maxIt = 10000, 
              adjmaxIt = NULL, K = 5, crit.cv = c("loglik", "AIC", "BIC"), 
              start = c("warm", "cold"), cores = 1, 
              trace = c("progress", "print", "none"),
              lam_prior = NULL, ...) {
  
  # match values
  crit.cv = match.arg(crit.cv)
  start = match.arg(start)
  trace = match.arg(trace)
  # lam = sort(lam)
  
  # initialize
  Path = NULL
  initmaxIt = maxIt
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
    if (K > 1) {
      
      # training set
      leave.out = ind[(1 + floor((k - 1) * n/K)):floor(k * 
                                                         n/K)]
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
    maxIt = initmaxIt
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
      
      # update diagonal elements of init, if necessary
      if (diagonal) {
        diag(init) = diag(S.train) + lam_
      }
      
      # compute the penalized likelihood precision matrix
      # estimator
      # if(!diagonal) diag(lam_prior) = 0
        
      GLASSO = glassoFast(S = S.train, rho = lam_*lam_prior, thr = tol, 
                          maxIt = maxIt, 
                          start = "warm", w.init = init, wi.init = initOmega, 
                          trace = FALSE, ...)
      
      if (start == "warm") {
        
        # option to save initial values for warm starts
        init = GLASSO$w
        initOmega = GLASSO$wi
        maxIt = adjmaxIt
        
      }
      
      # compute the observed negative validation loglikelihood
      # (close enoug)
      CV_errors[i, k] = (nrow(X)/2) * (sum(GLASSO$wi * S.valid) - determinant(GLASSO$wi, logarithm = TRUE)$modulus[1])
      
      # update for crit.cv, if necessary
      if (crit.cv == "AIC") {
        CV_errors[i, k] = CV_errors[i, k] + sum(GLASSO$wi != 0)
      }
      if (crit.cv == "BIC") {
        CV_errors[i, k] = CV_errors[i, k] + sum(GLASSO$wi != 0) * log(nrow(X))/2
      }
      
      # save estimate if path = TRUE
      if (path) {
        Path[, , i] = GLASSO$wi
      }
      
      # update progress bar
      if (trace == "progress") {
        setTxtProgressBar(progress, i + (k - 1) *  length(lam))
        
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


# ---- CVP -----

#' @param X nxp data matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
#' @param lam positive tuning parameters for elastic net penalty. If a vector of parameters is provided, they should be in increasing order. Defaults to grid of values \code{10^seq(-2, 2, 0.2)}.
#' @param diagonal option to penalize the diagonal elements of the estimated precision matrix (\eqn{\Omega}). Defaults to \code{FALSE}.
#' @param tol convergence tolerance. Iterations will stop when the average absolute difference in parameter estimates in less than \code{tol} times multiple. Defaults to 1e-4.
#' @param maxIt maximum number of iterations. Defaults to 1e4.
#' @param adjmaxIt adjusted maximum number of iterations. During cross validation this option allows the user to adjust the maximum number of iterations after the first \code{lam} tuning parameter has converged. This option is intended to be paired with \code{warm} starts and allows for 'one-step' estimators. Defaults to NULL.
#' @param K specify the number of folds for cross validation.
#' @param crit.cv cross validation criterion (\code{loglik}, \code{AIC}, or \code{BIC}). Defaults to \code{loglik}.
#' @param start specify \code{warm} or \code{cold} start for cross validation. Default is \code{warm}.
#' @param cores option to run CV in parallel. Defaults to \code{cores = 1}.
#' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.
#' @param ... additional arguments to pass to \code{glasso}.

# we define the CVP function
CVP = function(X = NULL, lam = 10^seq(-2, 2, 0.2), diagonal = FALSE, 
               tol = 1e-04, maxIt = 10000, adjmaxIt = NULL, K = 5, crit.cv = c("loglik", "AIC", "BIC"), 
               start = c("warm", "cold"), cores = 1, 
               trace = c("progress", "print", "none"), 
               lam_prior = NULL,...) {
  
  # match values
  crit.cv = match.arg(crit.cv)
  start = match.arg(start)
  trace = match.arg(trace)
  # lam = sort(lam)
  
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
  registerDoParallel(cluster)
  
  # use cluster for each fold in CV
  n = nrow(X)
  ind = sample(n)
  k = NULL
  CV = foreach(k = 1:K, .packages = "CVglasso", .combine = "cbind", 
               .inorder = FALSE) %dopar% {
                 
                 # set progress bar
                 if (trace == "progress") {
                   progress = txtProgressBar(max = length(lam), style = 3)
                 }
                 
                 # training set
                 leave.out = ind[(1 + floor((k - 1) * n/K)):floor(k * 
                                                                    n/K)]
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
                   # if (diagonal) {
                     diag(init) = diag(S.train) + lam_
                   # }
                   
                   # compute the penalized likelihood precision matrix
                   # estimator
                   # if(!diagonal) diag(lam_prior) = 0
                   
                   GLASSO = glassoFast(S = S.train, rho = lam_*lam_prior, thr = tol, 
                                       maxIt = maxIt,
                                       start = "warm", w.init = init, wi.init = initOmega, 
                                       trace = FALSE, ...)
                   
                   if (start == "warm") {
                     
                     # option to save initial values for warm starts
                     init = GLASSO$w
                     initOmega = GLASSO$wi
                     maxIt = adjmaxIt
                     
                   }
                   
                   # compute the observed negative validation loglikelihood
                   # (close enoug)
                   CV_error[i] = (nrow(X)/2) * (sum(GLASSO$wi * S.valid) - determinant(GLASSO$wi, logarithm = TRUE)$modulus[1])
                   
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

# ---- print.CVglasso -----

print.CVglasso = function(x, ...) {
  
  # print warning if maxIt reached
  if (x$maxIt <= x$Iterations) {
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

# ---- plot.CVglasso -----

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
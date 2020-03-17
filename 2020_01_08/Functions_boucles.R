############# FONCTION BOUCLE GREMLIN #############################

### Pas de simulation
### Une seule table

loopGLGR_nosim = function(lambda_seq, dat,  n_iter = 25){
  listLambda = list()
  
  for(lam in lambda_seq){
    cat("** Lambda : ", lam, "\n")
    
    list_res_glasso = list()
    list_res_gremlin = list()
    list_clusters = list()
    list_diff_iter_prec = list()
    list_link = list()
    
    gl0 = gl1 = gr0 = NULL
    
    gl0 = glassoFast(S = cov(dat), rho = matrix(lam, ncol = ncol(cov(dat)), nrow = nrow(cov(dat))))
    list_res_glasso[[length(list_res_glasso)+1]] = gl0
    
    Net = defineNetwork(gl0$w, typeInter = "adj", rowFG = "Features", colFG = "Features")
    gr0 = multipartiteBM(list_Net = list(Net), namesFG = "Features",
                         v_distrib = "gaussian", v_Kmin = 2, v_Kmax = 5, verbose = FALSE, initBM = FALSE,  nbCores = 4)
    
    list_res_gremlin[[length(list_res_gremlin)+1]] = gr0
    
    clustr = gr0$fittedModel[[1]]$paramEstim$Z[[1]]
    
    list_clusters[[length(list_clusters)+1]] = length(unique(clustr))
    
    distrib_groups = gr0$fittedModel[[1]]$paramEstim$list_theta$FeaturesFeatures
    set.seed(1992)
    link1 = abs(distrib_groups$mean)/max(abs(distrib_groups$mean))
    link = link1[clustr, clustr]
    list_link[[length(list_link)+1]] = 1-link
    
    gl1 = glassoFast(S = cov(dat), rho = lam*(1-link), start = "warm", wi.init = gl0$wi, w.init = gl0$w)
    # gl1 = glassoFast(S = cov(dat), rho = (1-link), start = "warm", wi.init = gl0$wi, w.init = gl0$w)
    
    list_res_glasso[[length(list_res_glasso)+1]] = gl1
    list_diff_iter_prec[[length(list_diff_iter_prec)+1]] = sum(abs(list_res_glasso[[length(list_res_glasso)]]$wi - list_res_glasso[[length(list_res_glasso) - 1 ]]$wi))
    
    ##### WHILE #####
    i = 1
    # while(list_diff_iter_prec[[length(list_diff_iter_prec)]]>=10^-3 && i<n_iter){
      while(list_diff_iter_prec[[length(list_diff_iter_prec)]]>=10^-6 && i<n_iter){
      # while(i<n_iter){ # Force a faire plusieurs iterations meme si convergence avant.
      cat("** Iter ", i, " - ")
      # gr = NULL
      Net = defineNetwork(list_res_glasso[[length(list_res_glasso)]]$w , typeInter = "adj", rowFG = "Features", colFG = "Features")
      gr =  tryCatch({ 
        multipartiteBM(list_Net = list(Net), namesFG = "Features",
                       v_distrib = "gaussian", v_Kmin = 2, v_Kmax = 5, verbose = FALSE, initBM = FALSE,  # initBM ne change rien.
                       # maxiterVE = 200, maxiterVEM = 200, 
                       nbCores = 4)
        
        
      }, error = function(e){ NULL }) # parfois error : `*tmp*`[[k]] : indice hors limites. Pas de convergence nulle part ?
      
      
      list_res_gremlin[[length(list_res_gremlin)+1]] = gr
      
      if(!is.null(gr)){ # convergence
        clustr = gr$fittedModel[[1]]$paramEstim$Z[[1]]
        
        list_clusters[[length(list_clusters)+1]] = length(unique(clustr))
        cat(list_clusters[[length(list_clusters)]], " - ")
        
        distrib_groups = gr$fittedModel[[1]]$paramEstim$list_theta$FeaturesFeatures
        link1 = abs(distrib_groups$mean)/max(abs(distrib_groups$mean))
        
      }else{ # non convergence de toute ce qui a ete teste ?
        distrib_groups =  list("mean" = matrix(0, 3, 3))
        link1 = distrib_groups$mean
        
        set.seed(i*1992)
        clustr = sample(c(1,2,3), size = ncol(list_res_glasso[[length(list_res_glasso)]]$wi), replace = TRUE)
        
        list_clusters[[length(list_clusters)+1]] = NA
        cat("DID NOT CONVERGE - ")
        
        list_NID[[length(list_NID)+1]] = NA
        
      }
      
      link = link1[clustr, clustr]
      list_link[[length(list_link)+1]] = 1-link
      
      if(!any(is.na(list_res_glasso[[length(list_res_glasso)]]$wi))){
        gl = glassoFast(S = cov(dat), rho = (1-link)*lam, start = "warm", wi.init = list_res_glasso[[length(list_res_glasso)]]$wi, w.init = list_res_glasso[[length(list_res_glasso)]]$w)
        # gl = glassoFast(S = cov(dat), rho = (1-link), start = "warm", wi.init = list_res_glasso[[length(list_res_glasso)]]$wi, w.init = list_res_glasso[[length(list_res_glasso)]]$w)
      }else{
        gl = glassoFast(S = cov(dat), rho = (1-link)*lam, start = "cold")
        # gl = glassoFast(S = cov(dat), rho = (1-link), start = "cold")
      }
      
      list_res_glasso[[length(list_res_glasso)+1]] = gl
      list_diff_iter_prec[[length(list_diff_iter_prec)+1]] = sum(abs(list_res_glasso[[length(list_res_glasso)]]$wi - list_res_glasso[[length(list_res_glasso) - 1 ]]$wi))
      i = i+1
      cat( " -- ")
    }
    cat("\n")
    
    listLambda[[length(listLambda) +1 ]] = list(
      "glasso" = list_res_glasso,
      "gremlin" = list_res_gremlin ,
      "clusters" = list_clusters,
      "diff_prev_iter" = list_diff_iter_prec,
      "link" = list_link
    )
  }
  return(listLambda)
}

############# FONCTION BOUCLE GREMLIN #############################

loopGLGR = function(lambda_seq, dat, Omega_sim, n_iter, classification){
  listLambda = list()
  
  for(lam in lambda_seq){
    cat("** Lambda : ", lam, "\n")
    
    list_res_glasso = list()
    list_res_gremlin = list()
    list_clusters = list()
    list_NID = list()
    list_diff_sim = list()
    list_diff_iter_prec = list()
    list_link = list()
    
    gl0 = gl1 = gr0 = NULL
    
    gl0 = glassoFast(S = cov(dat), rho = matrix(lam, ncol = ncol(cov(dat)), nrow = nrow(cov(dat))))
    list_res_glasso[[length(list_res_glasso)+1]] = gl0
    list_diff_sim[[length(list_diff_sim)+1]] = sum(abs(list_res_glasso[[length(list_res_glasso)]]$wi - Omega_sim))
    
    Net = defineNetwork(gl0$w, typeInter = "adj", rowFG = "Features", colFG = "Features")
    gr0 = multipartiteBM(list_Net = list(Net), namesFG = "Features",
                         v_distrib = "gaussian", v_Kmin = 2, v_Kmax = 5, verbose = FALSE, initBM = FALSE,  nbCores = 4)
    
    list_res_gremlin[[length(list_res_gremlin)+1]] = gr0
    
    clustr = gr0$fittedModel[[1]]$paramEstim$Z[[1]]
    
    list_clusters[[length(list_clusters)+1]] = length(unique(clustr))
    list_NID[[length(list_NID)+1]] = NID(clustr, classification)
    
    distrib_groups = gr0$fittedModel[[1]]$paramEstim$list_theta$FeaturesFeatures
    set.seed(1992)
    link1 = abs(distrib_groups$mean)/max(abs(distrib_groups$mean))
    link = link1[clustr, clustr]
    list_link[[length(list_link)+1]] = 1-link
    
    gl1 = glassoFast(S = cov(dat), rho = lam*(1-link), start = "warm", wi.init = gl0$wi, w.init = gl0$w)
    # gl1 = glassoFast(S = cov(dat), rho = (1-link), start = "warm", wi.init = gl0$wi, w.init = gl0$w)
    
    list_res_glasso[[length(list_res_glasso)+1]] = gl1
    list_diff_sim[[length(list_diff_sim)+1]] = sum(abs(list_res_glasso[[length(list_res_glasso)]]$wi - Omega_sim))
    list_diff_iter_prec[[length(list_diff_iter_prec)+1]] = sum(abs(list_res_glasso[[length(list_res_glasso)]]$wi - list_res_glasso[[length(list_res_glasso) - 1 ]]$wi))
    
    ##### WHILE #####
    i = 1
    # while(list_diff_iter_prec[[length(list_diff_iter_prec)]]>=10^-3 && i<n_iter){
    while(list_diff_iter_prec[[length(list_diff_iter_prec)]]>=10^-6 && i<n_iter){
    # while(i<n_iter){ # Force a faire plusieurs iterations meme si convergence avant.
      cat("** Iter ", i, " - ")
      # gr = NULL
      Net = defineNetwork(list_res_glasso[[length(list_res_glasso)]]$w , typeInter = "adj", rowFG = "Features", colFG = "Features")
      gr =  tryCatch({ 
        multipartiteBM(list_Net = list(Net), namesFG = "Features",
                                      v_distrib = "gaussian", v_Kmin = 2, v_Kmax = 5, verbose = FALSE, initBM = FALSE,  # initBM ne change rien.
                                      # maxiterVE = 200, maxiterVEM = 200, 
                                      nbCores = 4)
        
        
      }, error = function(e){ NULL }) # parfois error : `*tmp*`[[k]] : indice hors limites. Pas de convergence nulle part ?
      
      
      list_res_gremlin[[length(list_res_gremlin)+1]] = gr
      
      if(!is.null(gr)){ # convergence
        clustr = gr$fittedModel[[1]]$paramEstim$Z[[1]]
        
        list_clusters[[length(list_clusters)+1]] = length(unique(clustr))
        cat(list_clusters[[length(list_clusters)]], " - ")
        
        list_NID[[length(list_NID)+1]] = NID(clustr, classification)
        
        cat(round(list_NID[[length(list_NID)]],2)) 
        
        distrib_groups = gr$fittedModel[[1]]$paramEstim$list_theta$FeaturesFeatures
        link1 = abs(distrib_groups$mean)/max(abs(distrib_groups$mean))
        
      }else{ # non convergence de toute ce qui a ete teste ?
        distrib_groups =  list("mean" = matrix(0, 3, 3))
        link1 = distrib_groups$mean
        
        set.seed(i*1992)
        clustr = sample(c(1,2,3), size = ncol(list_res_glasso[[length(list_res_glasso)]]$wi), replace = TRUE)
        
        list_clusters[[length(list_clusters)+1]] = NA
        cat("DID NOT CONVERGE - ")
        
        list_NID[[length(list_NID)+1]] = NA
        
      }
      
      link = link1[clustr, clustr]
      list_link[[length(list_link)+1]] = 1-link
      
      if(!any(is.na(list_res_glasso[[length(list_res_glasso)]]$wi))){
        gl = glassoFast(S = cov(dat), rho = (1-link)*lam, start = "warm", wi.init = list_res_glasso[[length(list_res_glasso)]]$wi, w.init = list_res_glasso[[length(list_res_glasso)]]$w)
        # gl = glassoFast(S = cov(dat), rho = (1-link), start = "warm", wi.init = list_res_glasso[[length(list_res_glasso)]]$wi, w.init = list_res_glasso[[length(list_res_glasso)]]$w)
      }else{
        gl = glassoFast(S = cov(dat), rho = (1-link)*lam, start = "cold")
        # gl = glassoFast(S = cov(dat), rho = (1-link), start = "cold")
      }
      
      list_res_glasso[[length(list_res_glasso)+1]] = gl
      list_diff_sim[[length(list_diff_sim)+1]] = sum(abs(list_res_glasso[[length(list_res_glasso)]]$wi - Omega_sim))
      list_diff_iter_prec[[length(list_diff_iter_prec)+1]] = sum(abs(list_res_glasso[[length(list_res_glasso)]]$wi - list_res_glasso[[length(list_res_glasso) - 1 ]]$wi))
      i = i+1
      cat( " -- ")
    }
    cat("\n")
    
    listLambda[[length(listLambda) +1 ]] = list(
      "glasso" = list_res_glasso,
      "gremlin" = list_res_gremlin ,
      "clusters" = list_clusters,
      "NID" = list_NID,
      "diff_sim" = list_diff_sim,
      "diff_prev_iter" = list_diff_iter_prec,
      "link" = list_link
    )
  }
  return(listLambda)
}

############# FONCTION BM STEP #######################################
# ---- BM step : ----

BM_gaussian_step = function(
  membership_type, 
  adj, 
  verbosity=0,
  autosave='',
  plotting='',
  exploration_factor=1.5,
  explore_min=4,
  explore_max=Inf,
  ncores=detectCores()){
  
  if(identical((adj - diag(diag(adj))), matrix(0, ncol(adj), ncol(adj)))){
    cat("Glasso estimated a diagonal matrix, no groups to find here. Skipping next steps.\n")
    
    BM_model = BM_gaussian(
      membership_type = membership_type, 
      explore_min = explore_min, 
      explore_max = explore_max,
      exploration_factor = exploration_factor,
      autosave = autosave,
      plotting = plotting,
      verbosity = verbosity,
      adj = adj,
      ncores = ncores)
    
    mat_penalty = matrix(1, ncol(adj), ncol(adj))
    
    
  }else{
    BM_model = BM_gaussian(
      membership_type = membership_type, 
      explore_min = explore_min, 
      explore_max = explore_max,
      exploration_factor = exploration_factor,
      autosave = autosave,
      plotting = plotting,
      verbosity = verbosity,
      adj = adj,
      ncores = ncores)
    
    BM_model$estimate()
    if(length(BM_model$memberships)>0){
      groups = apply(BM_model$memberships[[which.max(BM_model$ICL)]]$Z, 1, which.max)
      
      mu_mod = BM_model$model_parameters[[which.max(BM_model$ICL)]]$mu
      penalisation_groupe = 1-abs(mu_mod)/max(abs(mu_mod)) 
      mat_penalty = penalisation_groupe[groups, groups]
    }else{
      mat_penalty = matrix(1, ncol(adj), ncol(adj))
    }
  }
  
  return(list(BM_model= BM_model, mat_penalty = mat_penalty))
}

############# FONCTION BOUCLE BlockMODEL #############################

loopGLBM = function(lambda_seq, dat, Omega_sim, n_iter, classification){
  listLambda = list()
  
  for(lam in lambda_seq){
    cat("** Lambda : ", lam, "\n")
  
    # lam = lambda_seq[1]
      
    list_res_glasso = list()
    list_res_BM = list()
    list_clusters = list()
    list_NID = list()
    list_diff_sim = list()
    list_diff_iter_prec = list()
    list_link = list()
    
    gl0 = gl1 = gr0 = NULL
    
    gl0 = glassoFast(S = cov(dat), rho = matrix(lam, ncol = ncol(cov(dat)), nrow = nrow(cov(dat))))
    list_res_glasso[[length(list_res_glasso)+1]] = gl0
    list_diff_sim[[length(list_diff_sim)+1]] = sum(abs(list_res_glasso[[length(list_res_glasso)]]$wi - Omega_sim))
    
    bm0 = BM_gaussian_step("SBM", adj = gl0$w, verbosity = 0, exploration_factor = 1.5, explore_min = 2, explore_max = 5, ncores = 4)
    
    list_res_BM[[length(list_res_BM)+1]] = bm0
    
    clustr = apply(bm0$BM_model$memberships[[which.max(bm0$BM_model$ICL)]]$Z, 1, which.max)
    
    mu_mod = bm0$BM_model$model_parameters[[which.max(bm0$BM_model$ICL)]]$mu
    penalisation_groupe = 1-abs(mu_mod)/max(abs(mu_mod)) 
    mat_penalty = penalisation_groupe[clustr, clustr]
    
    list_clusters[[length(list_clusters)+1]] = length(unique(clustr))
    list_NID[[length(list_NID)+1]] = NID(clustr, classification)
    link = 1-list_res_BM[[length(list_res_BM)]]$mat_penalty
    
    gl1 = glassoFast(S = cov(dat), rho = lam*(1-link), start = "warm", wi.init = gl0$wi, w.init = gl0$w)
    # gl1 = glassoFast(S = cov(dat), rho = (1-link), start = "warm", wi.init = gl0$wi, w.init = gl0$w)
    
    list_res_glasso[[length(list_res_glasso)+1]] = gl1
    list_diff_sim[[length(list_diff_sim)+1]] = sum(abs(list_res_glasso[[length(list_res_glasso)]]$wi - Omega_sim))
    list_diff_iter_prec[[length(list_diff_iter_prec)+1]] = sum(abs(list_res_glasso[[length(list_res_glasso)]]$wi - list_res_glasso[[length(list_res_glasso) - 1 ]]$wi))
    
    ##### WHILE #####
    i = 1
    while(list_diff_iter_prec[[length(list_diff_iter_prec)]]>10^-3 && i<n_iter){
    # while(list_diff_iter_prec[[length(list_diff_iter_prec)]]>=10^-6 && i<n_iter){ 
    # while(i<n_iter){ # Force a faire plusieurs iterations meme si convergence avant.
      cat("** Iter ", i, " - ")
      
      bm = BM_gaussian_step("SBM", adj = list_res_glasso[[length(list_res_glasso)]]$w, verbosity = 0, exploration_factor = 1.5, explore_min = 2, explore_max = 5, ncores = 4)
      
      list_res_BM[[length(list_res_BM)+1]] = bm
      
      clustr = apply(bm$BM_model$memberships[[which.max(bm$BM_model$ICL)]]$Z, 1, which.max)
      
      mu_mod = bm$BM_model$model_parameters[[which.max(bm$BM_model$ICL)]]$mu
      penalisation_groupe = 1-abs(mu_mod)/max(abs(mu_mod)) 
      mat_penalty = penalisation_groupe[clustr, clustr]
      
      list_clusters[[length(list_clusters)+1]] = length(unique(clustr))
      list_NID[[length(list_NID)+1]] = NID(clustr, classification)
      
      link = 1-list_res_BM[[length(list_res_BM)]]$mat_penalty
      list_link[[length(list_link)+1]] = 1-link
      
      if(!any(is.na(list_res_glasso[[length(list_res_glasso)]]$wi))){
        gl = glassoFast(S = cov(dat), rho = (1-link)*lam, start = "warm", wi.init = list_res_glasso[[length(list_res_glasso)]]$wi, w.init = list_res_glasso[[length(list_res_glasso)]]$w)
        # gl = glassoFast(S = cov(dat), rho = (1-link), start = "warm", wi.init = list_res_glasso[[length(list_res_glasso)]]$wi, w.init = list_res_glasso[[length(list_res_glasso)]]$w)
      }else{
        gl = glassoFast(S = cov(dat), rho = (1-link)*lam, start = "cold")
        # gl = glassoFast(S = cov(dat), rho = (1-link), start = "cold")
      }
      
      list_res_glasso[[length(list_res_glasso)+1]] = gl
      list_diff_sim[[length(list_diff_sim)+1]] = sum(abs(list_res_glasso[[length(list_res_glasso)]]$wi - Omega_sim))
      list_diff_iter_prec[[length(list_diff_iter_prec)+1]] = sum(abs(list_res_glasso[[length(list_res_glasso)]]$wi - list_res_glasso[[length(list_res_glasso) - 1 ]]$wi))
      i = i+1
      cat( " -- ")
    }
    cat("\n")
    
    listLambda[[length(listLambda) +1 ]] = list(
      "glasso" = list_res_glasso,
      "BM" = list_res_BM ,
      "clusters" = list_clusters,
      "NID" = list_NID,
      "diff_sim" = list_diff_sim,
      "diff_prev_iter" = list_diff_iter_prec,
      "link" = list_link
    )
  }
  return(listLambda)
}

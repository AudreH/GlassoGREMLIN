##### PACKAGES ####

library(glassoFast)
library(GREMLIN)
library(blockmodels)
library(parallel)
library(reshape2)
library(aricode)
library(R.utils)

library(ggalluvial)
library(tidyverse)
library(reshape2)

############# FONCTION BOUCLE GREMLIN #############################

loopGLGR = function(lambda_seq, dat = NULL, Sigma = NULL, Omega_sim, n_iter, classification, threshold =10^-6,
                    timeout = Inf, diagonal = TRUE, min_pen = 0, ncores = detectCores()) {
  listLambda = list()
  
  names(lambda_seq) = lambda_seq
  listLambda = mclapply(lambda_seq, FUN = function(lam){

    # for(lam in lambda_seq){
    start_time = Sys.time()
    cat("** Lambda : ", lam, "\n")
    
    # lam = lambda_seq[1]
    # # lam = 10^-4
    # Sigma = Sigma
    # dat = as.matrix(X)
    # n_iter = 10
    # Omega_sim = Omega
    # ncores = 2
    # diagonal = TRUE

    if(!is.null(dat)) Sigma = cov(dat)
    
    n = nrow(Sigma)
    if(!is.null(dat)) n = nrow(dat)
    
    C = matrix(1, nrow(Sigma), ncol(Sigma)) 
    if(!diagonal) C = C - diag(1, nrow(Sigma))
    
    
    list_res_glasso = list()
    list_res_gremlin = list()
    list_clusters = list()
    list_classif = list()
    list_NID = list()
    list_diff_sim = list()
    list_diff_iter_prec = list()
    list_link = list()
    list_logLik = list()
    list_AIC = list()
    list_BIC = list()
    list_ICL = list()
    exit_status1 = "Did all iterations"
    exit_status = 1
    
    gl0 = gl1 = gr0 = NULL
    
    Rho = C*matrix(lam, ncol = ncol(Sigma), nrow = nrow(Sigma))
    gl0 = glassoFast(S = Sigma, rho = Rho)
    
    
    list_res_glasso[[length(list_res_glasso)+1]] = gl0
    list_diff_sim[[length(list_diff_sim)+1]] = sum(abs(gl0$wi - Omega_sim))
    list_logLik[[length(list_logLik)+1]] = (n/2) * (sum(gl0$wi*Sigma) - determinant(gl0$wi, logarithm = TRUE)$modulus[1] +
                                                       sum(Rho*C*abs(gl0$wi)) 
    )
    list_AIC[[length(list_AIC)+1]] = list_logLik[[length(list_logLik)]] + sum(gl0$wi!=0)
    list_BIC[[length(list_BIC)+1]] = list_logLik[[length(list_logLik)]] + sum(gl0$wi!=0)*log(n/2)
    
    Net = defineNetwork(gl0$w, typeInter = "adj", rowFG = "Features", colFG = "Features")
    gr0 = tryCatch(expr = { 
      withTimeout( { multipartiteBM(list_Net = list(Net), namesFG = "Features",
                                    v_distrib = "gaussian", v_Kmin = 2, v_Kmax = 5, verbose = FALSE, 
                                    initBM = FALSE, nbCores = 2) },
                   timeout = timeout) },
      error = function(e){ NULL })
    
    if(!is.null(gr0)){
      list_res_gremlin[[length(list_res_gremlin)+1]] = gr0
      
      clustr = gr0$fittedModel[[1]]$paramEstim$Z[[1]]
      list_classif[[length(list_classif) +1 ]] = clustr
      list_clusters[[length(list_clusters)+1]] = length(unique(clustr))
      list_NID[[length(list_NID)+1]] = NID(clustr, classification)
      
      list_ICL[[length(list_ICL)+1]] = gr0$fittedModel[[1]]$ICL
      
      distrib_groups = gr0$fittedModel[[1]]$paramEstim$list_theta$FeaturesFeatures
      set.seed(1992)
      denom = ifelse(max(abs(distrib_groups$mean))!=0, max(abs(distrib_groups$mean)), 1)
      link1 = abs(distrib_groups$mean)/denom
      link = link1[clustr, clustr]
      list_link[[length(list_link)+1]] = 1-link
      
      # Mettre 0 en pénalisation, j'ai l'impression que ça fait des bugs.
      # 10^-5 minimum admissible ?
      Rho = pmax(matrix(min_pen, ncol(link), ncol(link)), lam*(1-link))*C
      gl1 = glassoFast(S = Sigma, rho = Rho, start = "warm", wi.init = gl0$wi, w.init = gl0$w)
      # gl1 = glassoFast(S = Sigma, rho = (1-link), start = "warm", wi.init = gl0$wi, w.init = gl0$w)
      
      # Pourquoi ça donne des NAn ?
      
      if(!any(is.na(gl1$wi))){
        list_res_glasso[[length(list_res_glasso)+1]] = gl1
        list_diff_sim[[length(list_diff_sim)+1]] = sum(abs(gl1$wi - Omega_sim))
        list_diff_iter_prec[[length(list_diff_iter_prec)+1]] = sum(abs(gl1$wi - list_res_glasso[[length(list_res_glasso) - 1 ]]$wi))
        list_logLik[[length(list_logLik)+1]] = (n/2) * (sum(gl1$wi*Sigma) - determinant(gl1$wi, logarithm = TRUE)$modulus[1] +
                                                           sum(Rho*C*abs(gl1$wi)) 
        )
        list_AIC[[length(list_AIC)+1]] = list_logLik[[length(list_logLik)]] + sum(gl1$wi*C!=0)
        list_BIC[[length(list_BIC)+1]] = list_logLik[[length(list_logLik)]] + sum(gl1$wi*C!=0)*log(n/2)
        
        ##### WHILE #####
        i = 0
        while(list_diff_iter_prec[[length(list_diff_iter_prec)]]>=threshold && i<n_iter){
          i = i+1
          cat("** Iter ", i, " - ")
          # gr = NULL
          Net = defineNetwork(list_res_glasso[[length(list_res_glasso)]]$w , typeInter = "adj", rowFG = "Features", colFG = "Features")
          gr =  tryCatch({ 
            withTimeout( {multipartiteBM(list_Net = list(Net), namesFG = "Features",
                                         v_distrib = "gaussian", v_Kmin = 2, v_Kmax = 5, verbose = FALSE, initBM = FALSE,  # initBM ne change rien.
                                         # maxiterVE = 200, maxiterVEM = 200, 
                                         nbCores = 4)},
                         timeout = timeout) 
          }, error = function(e){ NULL }) # parfois error : `*tmp*`[[k]] : indice hors limites. Pas de convergence nulle part ?
          
          
          list_res_gremlin[[length(list_res_gremlin)+1]] = gr
          
          if(!is.null(gr)){ # convergence
            list_ICL[[length(list_ICL)+1]] = gr$fittedModel[[1]]$ICL
            
            clustr = gr$fittedModel[[1]]$paramEstim$Z[[1]]
            list_classif[[length(list_classif) +1 ]] = clustr
            
            list_clusters[[length(list_clusters)+1]] = length(unique(clustr))
            cat(list_clusters[[length(list_clusters)]], " - ")
            
            list_NID[[length(list_NID)+1]] = NID(clustr, classification)
            
            cat(round(list_NID[[length(list_NID)]],2)) 
            
            distrib_groups = gr$fittedModel[[1]]$paramEstim$list_theta$FeaturesFeatures
            link1 = abs(distrib_groups$mean)/max(abs(distrib_groups$mean))
            
            link = link1[clustr, clustr]
            list_link[[length(list_link)+1]] = 1-link
            
            Rho = pmax(matrix(min_pen, ncol(link), ncol(link)), lam*(1-link))*C
            gl = glassoFast(S = Sigma, rho = Rho, start = "warm",
                            wi.init = list_res_glasso[[length(list_res_glasso)]]$wi, 
                            w.init = list_res_glasso[[length(list_res_glasso)]]$w)
            
            
            if(!any(is.na(gl$wi))){
              
              list_res_glasso[[length(list_res_glasso)+1]] = gl
              list_diff_sim[[length(list_diff_sim)+1]] = sum(abs(gl$wi - Omega_sim))
              list_diff_iter_prec[[length(list_diff_iter_prec)+1]] = sum(abs(gl$wi - list_res_glasso[[length(list_res_glasso) - 1 ]]$wi))
              
              list_logLik[[length(list_logLik)+1]] = (n/2) * (sum(gl$wi*Sigma) - determinant(gl$wi, logarithm = TRUE)$modulus[1] +
                                                                 sum(Rho*abs(gl$wi)) 
              )
              list_AIC[[length(list_AIC)+1]] = list_logLik[[length(list_logLik)]] + sum(gl$wi*C!=0)
              list_BIC[[length(list_BIC)+1]] = list_logLik[[length(list_logLik)]] + sum(gl$wi*C!=0)*log(n/2)
              cat( " -- ")
              
            }else{
              cat(paste0("Glasso did not converge in interation ", i , "exiting from the loop.\n"))
              exit_status1 = paste0("Exited from loop at iter ", i, " - Glasso did not converge")
              exit_status = 0
              i = n_iter
            } # ENDIF any gl$wi == na
            
          }else{ # non convergence de toute ce qui a ete teste ?
            distrib_groups =  list("mean" = matrix(0, 3, 3))
            
            list_ICL[[length(list_ICL) +1 ]] = NA
            
            list_clusters[[length(list_clusters)+1]] = NA
            cat("DID NOT CONVERGE - ")
            
            list_NID[[length(list_NID)+1]] = NA
            
            cat(paste0("GREMLIN did not converge in interation ", i , " exiting from the loop.\n"))
            exit_status1 = paste0("Exited from loop at iter ", i, " - GREMLIN did not converge")
            exit_status = 0
            i = n_iter
            
          }
          
          if(list_diff_iter_prec[[length(list_diff_iter_prec)]]<=threshold) exit_status1 = paste0("Successfully converged at iter ", i)
          
        } # END WHILE 
        cat("\n")
      }else{
        cat("Glasso did not converge after first try of GREMLIN, did not entered the loop.\n")
        exit_status1 = paste0("Exited from algorithm before the loop - Glasso 1 did not converge")
        exit_status = 0
      }# ENDIF any gl1$wi == na
    }else{# ENDIF gr0 == NULL
      cat("GREMLIN did not converge on first try. Did not entered the loop.\n")
      exit_status1 = paste0("Exited from algorithm before the loop - GREMLIN did not converge")
      exit_status = 0
    }
    
    end_time = Sys.time()
    
    # listLambda[[length(listLambda) +1 ]] = list(
    #   "glasso" = list_res_glasso,
    #   "gremlin" = list_res_gremlin ,
    #   "clusters" = list_clusters,
    #   "classif" = list_classif,
    #   "NID" = list_NID,
    #   "diff_sim" = list_diff_sim,
    #   "diff_prev_iter" = list_diff_iter_prec,
    #   "link" = list_link,
    #   "logLik" = list_logLik,
    #   "AIC" = list_AIC,
    #   "BIC" = list_BIC,
    #   "ICL" = list_ICL,
    #   "n_iter" = n_iter,
    #   "dat" = dat,
    #   "Sigma" = Sigma,
    #   "Omega_sim" = Omega_sim,
    #   "exit_status" = exit_status,
    #   "time" = list("start" = start_time, "end" = end_time)
    # )

  list(
    "glasso" = list_res_glasso,
    "gremlin" = list_res_gremlin ,
    "clusters" = list_clusters,
    "classif" = list_classif,
    "NID" = list_NID,
    "diff_sim" = list_diff_sim,
    "diff_prev_iter" = list_diff_iter_prec,
    "link" = list_link,
    "logLik" = list_logLik,
    "AIC" = list_AIC,
    "BIC" = list_BIC,
    "ICL" = list_ICL,
    "n_iter" = n_iter,
    "dat" = dat,
    "Sigma" = Sigma,
    "Omega_sim" = Omega_sim,
    "exit_status1" = exit_status1,
    "exit_status" = exit_status,
    "time" = list("start" = start_time, "end" = end_time)
  )
  }, mc.cores = ncores) # end mclapply
  # } # END FOR lambda sequence
  names(listLambda) = lambda_seq
  return(listLambda)
}

# ---- Function res graph plot : ----

res_sim  = function(res, n_iter, lambda_seq){
  df_test = do.call("rbind",lapply(1:length(res), FUN = function(i){
    # i = 1
    liste = res[[i]]
    df = data.frame(t(do.call("cbind", list("Clusters" = c(NA, unlist(liste$clusters), rep(NA, (n_iter+1)-length(unlist(liste$clusters)))),
                                            "NID" = c(NA, unlist(liste$NID), rep(NA, (n_iter+1)-length(unlist(liste$NID)))),
                                            "diff_sim" = c(unlist(liste$diff_sim), rep(NA, 1 + (n_iter+1)-length(unlist(liste$diff_sim)))),
                                            "diff_prev_iter" = c(NA, unlist(liste$diff_prev_iter), rep(NA, (n_iter+1)-length(liste$diff_prev_iter)))))))
    colnames(df) = c("GL0", "GRGL0", paste0("GRGL", 1:(ncol(df)-2)))
    df$response = rownames(df)
    df$lambda = lambda_seq[i] 
    df
  }))
  
  df_melt = melt(df_test, id.vars = c("lambda", "response"))
  df_melt$lambda = factor(df_melt$lambda)
  df_melt = df_melt[!is.na(df_melt$value),]
  # df_melt$lambda = factor(round(as.numeric(as.character(df_melt$lambda)),2))
  
  df_melt
}

############# FONCTION ALLUVIAL LISTE #############################


res.alluv = function(list_classif, classif = NULL, title = ""){
  df_alluvial = data.frame(do.call(cbind, list_classif))
  colnames(df_alluvial) = c("GR0", paste0("GRGL", 1:(ncol(df_alluvial)-1)))
  if(!is.null(classif)) df_alluvial$classif = classif
  df_alluvial$Ind = 1:nrow(df_alluvial)
  
  
  df_melt = melt(df_alluvial, id.vars = "Ind")
  df_melt$value = factor(df_melt$value)
  colnames(df_melt) = c("Ind", "Iteration", "Groupe")
  
  ggplot(df_melt,
         aes(x = Iteration, stratum = Groupe, alluvium = Ind,
             fill = Groupe, label = Groupe)) +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    geom_flow(stat = "alluvium", lode.guidance = "frontback",
              color = "darkgray") +
    geom_stratum() +
    theme(legend.position = "bottom")  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(title)
}

############# FONCTION PLOT CRITERIA #############################

res.plot = function(res_lam, title = ""){
  df_plot = data.frame("logLik" = unlist(res_lam$logLik),  # Glasso
                       "AIC" = unlist(res_lam$AIC), # Glasso
                       "BIC" = unlist(res_lam$BIC), # Glasso
                       "ICL" = c(NA, unlist(res_lam$ICL))  # GREMLIN
  )
  df_plot$Iter = factor(c("GL0", paste0("GRGL", 1:(nrow(df_plot)-1))), 
                        levels = c("GL0", paste0("GRGL", 1:(nrow(df_plot)-1))))
  df_plot = melt(df_plot, id.vars = "Iter")
  colnames(df_plot) = c("Iter", "Criterion", "Value")
  df_plot$Method = ifelse(df_plot$Criterion=="ICL", "GREMLIN", "Glasso")
  
  df_plot = df_plot[which(!is.na(df_plot$Value)),]
  
  ggplot(df_plot,
         aes(x = Iter, y = Value, colour = Criterion, group = Criterion)) +
    geom_point() +  geom_line() + facet_grid(Method~., scales = "free") + 
    theme(legend.position = "bottom")  + theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
    ggtitle(title)
}


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

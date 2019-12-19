# ---- BM step : ----

BM_gaussian_step = function(
  membership_type, 
  adj, 
  verbosity=0,
  autosave='',
  plotting=character(0),
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
      adj = adj)
    
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
      adj = adj)
    
    BM_model$estimate()
    if(length(BM_model$memberships)>0){
      groups = apply(BM_model$memberships[[which.max(BM_model$ICL)]]$Z,1, which.max)
      
      mu_mod = BM_model$model_parameters[[which.max(BM_model$ICL)]]$mu
      # penalisation_groupe = 1/abs(mu_mod) # pour l'instant
      penalisation_groupe = 1-abs(mu_mod)/max(abs(mu_mod)) 
      mat_penalty = penalisation_groupe[groups, groups]
    }else{
      mat_penalty = matrix(1, ncol(adj), ncol(adj))
    }
  }
  
  
  
  
  return(list(BM_model= BM_model, mat_penalty = mat_penalty))
}
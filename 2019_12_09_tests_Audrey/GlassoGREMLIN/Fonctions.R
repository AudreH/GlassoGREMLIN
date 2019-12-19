### 19/11/19
### Ajout 02/12/19
# source("Fonctions_CV_glasso.R")


#### Fonction automatique GREMLIN

# qu'est-ce qui change par rapport à juste faire du gremlin ?
# ça ressort la matrice reconstituée des probas de connexions... 
# ça fait la liste des networks à partir d'une matrice entière
# est-ce bien utile ?

# param matrix_all = matrice entière
# param nb_feat = vecteur indiquant le nombre de features par block (utile pour diviser la matrice)
# blocks_names = vecteur indiquant le nom des blocks. Par defaut, pris sur le vectuer du nb_feat
# 

gremlin_list = function(A_all, nb_feat, blocks_names = names(nb_feat),
                        v_Kmin = rep(2, length(nb_feat)), v_Kmax, 
                        verbose = TRUE, plot = FALSE, nbCores = 1,
                        maxiterVE = 100){
  
  # A_all = A_all
  # nb_feat = n_f_net
  # blocks_names = paste0("Tab", 1:n_networks)
  # v_Kmin = rep(2, n_networks)
  # v_Kmax =  (rep(nb_groupes+2, n_networks))
  # verbose = TRUE
  # nbCores = 1
  # maxiterVE = 100
   
  #### Decoupage de la matrice en blocks   #### 
  colnames(A_all) = 
    rownames(A_all) = unlist(lapply(1:length(nb_feat), FUN = function(x){
      rep(blocks_names[x], times = nb_feat[x])
    }))
  
  
  # "Applatissement" de la matrice
  A_M_all = matrix(NA, ncol = 4, nrow = length(A_all))
  compteur = 1
  for(i in 1:nrow(A_all)){
    for(j in 1:ncol(A_all)){
      if(rownames(A_all)[i] == colnames(A_all)[j]) name_tab = rownames(A_all)[i]
      else name_tab = paste(rownames(A_all)[i], colnames(A_all)[j], sep = "_")
      A_M_all[compteur,] = c(name_tab, i, j, A_all[i,j])
      compteur = compteur +1
    }
  }
  
  # Creation de la liste de matrices
  listMat = by(A_M_all, INDICES = list(A_M_all[,1]), FUN = function(x){
    # x1 = x
    x = x[,-1]
    x = apply(x, 2, FUN = function(vect) as.numeric(as.character(vect)))
    # changement des indices pour se ramener a des indices de matrices neuves... 
    x[,1] = as.numeric(x[,1]) - min(as.numeric(x[,1]))+1 
    x[,2] = as.numeric(x[,2]) - min(as.numeric(x[,2]))+1
    
    mat_res = matrix(NA, nrow = length(unique(x[,1])), ncol = length(unique(x[,2])))
    for(i in 1:nrow(x)){
      mat_res[x[i,1], x[i,2]] = x[i,3]
    }
    return(mat_res)
  })
  
  #### Definition de la liste de reseaux a utiliser dans GREMLIN  #### 
  listNetworks = lapply(1:length(listMat), FUN = function(indice){
    nam = strsplit(names(listMat)[indice], "_")
    if(length(nam[[1]])==1) return(defineNetwork(listMat[[indice]], typeInter = "diradj", rowFG = nam[[1]], colFG = nam[[1]]))
    else return(defineNetwork(listMat[[indice]], typeInter = "inc", rowFG = nam[[1]][1], colFG = nam[[1]][2]))
    
  })
  names(listNetworks) = names(listMat)
  
  #### Statistiques sur les differents morceaux du réseau #### 
  list_stat = lapply(listNetworks, FUN = function(mat){
    if(ncol(mat$mat) == nrow(mat$mat)){
      gr = graph_from_adjacency_matrix(mat$mat)
      return(list("betweenness" = betweenness(gr), "density" = graph.density(gr), "degree" = degree(gr)))
    }else{
      deg = c(rowSums(mat$mat), colSums(mat$mat))
      names(deg) = c(rownames(mat), colnames(mat))
      return(list("density" = sum(mat$mat)/length(mat$mat), "degree" = deg))
    }
  })
  
  #### Application GREMLIN  #### 
  res_GREMLIN = multipartiteBM(listNetworks, 
                               v_Kmin = v_Kmin,
                               v_Kmax = v_Kmax, 
                               namesFG = blocks_names,
                               verbose = verbose,
                               nbCores = nbCores
                               # maxiterVE = maxiterVE # Fait tout planter.
                               )
  
  if(plot) plotMBM(res_GREMLIN)
  res_GREMLIN$fittedModel[[1]]$paramEstim$v_K
  class_networks = extractClustersMBM(res_GREMLIN) 
  
  classification_network = lapply(class_networks, FUN = function(x){
    classif_return = rep(NA, sum(unlist(lapply(x, length))))
    for(i in 1:length(x)){
      classif_return[x[[i]]] = i
    }
    return(classif_return)
  })
  
  #### Recuperation des probabilites de connexions  #### 
  list_theta = res_GREMLIN$fittedModel[[1]]$paramEstim$list_theta 
  names(list_theta) = names(listNetworks)
  
  groupes = res_GREMLIN$fittedModel[[1]]$paramEstim$Z
  names(groupes) = blocks_names
  
  groupes = unlist(lapply(1:length(groupes), FUN = function(x){
    paste(names(groupes)[x], groupes[[x]], sep = "_")
  }))
  
  matrix_reconstruct = matrix(NA, ncol = length(groupes), nrow = length(groupes))
  colnames(matrix_reconstruct) = rownames(matrix_reconstruct) = groupes
  
  for(i in 1:nrow(matrix_reconstruct)){ # A changer, c'est moche et long. 
    # i = 1
    matrix_reconstruct[i,] = unlist(lapply(1:ncol(matrix_reconstruct), FUN = function(j){
      
      nam1 = strsplit(rownames(matrix_reconstruct)[i], split = "_")[[1]]
      nam2 = strsplit(colnames(matrix_reconstruct)[j], split = "_")[[1]]
      
      if(nam1[1]!=nam2[1]){
        nam_tab = paste(nam1[1], nam2[1], sep = "_")
        theta_mat = list_theta[which(names(list_theta)==nam_tab)][[1]]
        matrix_reconstruct[i,j] = theta_mat[as.numeric(nam1[2]), as.numeric(nam2[2])]
      }else{
        nam_tab = nam1[1]
        theta_mat = list_theta[which(names(list_theta)==nam_tab)][[1]]
        matrix_reconstruct[i,j] = theta_mat[as.numeric(nam1[2]), as.numeric(nam2[2])]
      }
    }))
  }
  
  
  
  #### Sortie  #### 
  return(list(classification_network = classification_network,
              matrix_reconstruct = matrix_reconstruct, 
              res_GREMLIN = res_GREMLIN,
              listNetworks = listNetworks, 
              listMat = listMat
  ))
  
  
}

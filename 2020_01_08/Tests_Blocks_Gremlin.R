rm(list = ls())
setwd("/home/hulot/Documents/packages_R/GlassoGREMLIN/2020_01_08/")

# ---- Packages ----
library(Matrix)
library(GREMLIN)
# library(simone)
library(aricode)
library(glassoFast)

library(ggplot2)

library(MASS)

library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(microbenchmark)

# ---- Simulation données : ----

n_ind = 10000
n_feat = 30
nb_groupes = 3 # nombre de groupes de features

set.seed(1992)
prop_groupe = 0.20 # proportion de features dans les groupes (sauf dernier groupe qui prend le reste)
classification = base::sample(1:nb_groupes, size = n_feat, replace = TRUE, c(rep(prop_groupe, nb_groupes-1), 1-prop_groupe*(nb_groupes-1))) 
table(classification) 

Sigma_sim = matrix(NA, n_feat, n_feat) # covariance structure
for(i in 1:(nb_groupes)){
  set.seed(i*2000)
  Sigma_sim[which(classification==i), which(classification == i)] = runif(n = length(Sigma_sim[which(classification==i),which(classification == i)]), min = 0.5, max = 1)
}
Sigma_sim[is.na(Sigma_sim)] = 0
# Sigma_sim[is.na(Sigma_sim)] = runif(n = length(Sigma_sim[is.na(Sigma_sim)]), min = 0, max = 5)
Chol = matrix(0, ncol = ncol(Sigma_sim), nrow = ncol(Sigma_sim))
Chol[upper.tri(Chol, diag = TRUE)] = Sigma_sim[upper.tri(Sigma_sim, diag = TRUE)]
Chol = apply(Chol, 1, FUN = function(vect) if(sum(vect == 0) != length(vect)){sqrt(vect^2 / sum(vect^2))}else{rep(0, length(vect))})
Sigma_sim = t(Chol)%*%Chol

dat = mvrnorm(n = n_ind, mu = c(1:n_feat)*classification, Sigma = Sigma_sim)
dat_c = scale(dat, center = TRUE, scale = FALSE) # scale = FALSE ? est-ce utile de toute façon ?

# ---- Annotations heatmaps : ----

column_ha = HeatmapAnnotation(Classification = as.factor(classification),
                              col = list("Classification" = c("1" = "blue", "2" = "azure3", "3" = "darkorchid3", "4" = "darksalmon", "5" = "darkolivegreen3")),
                              annotation_legend_param = list(nrow = 1)
)

row_ha = HeatmapAnnotation(Classification = as.factor(classification),
                           col = list("Classification" = c("1" = "blue", "2" = "azure3", "3" = "darkorchid3", "4" = "darksalmon", "5" = "darkolivegreen3")), 
                           which = "row",
                           show_annotation_name = FALSE, show_legend = FALSE,  annotation_legend_param = list(nrow = 1))

# ---- Visualisation des données : ----

hc = hclust(dist(Sigma_sim), method = "complete")
ComplexHeatmap::draw(Heatmap(Sigma_sim, top_annotation = column_ha, name = "Cov", column_title = "Covariance pattern used to simulate the data",
                             col = colorRamp2(c(0,  max(Sigma_sim)), c("white",  "red")), 
                             right_annotation = row_ha, 
                             row_dend_reorder = hc$order, column_dend_reorder = hc$order,
                             cluster_rows = hc, cluster_columns = hc,
                             heatmap_legend_param = list(direction = "horizontal"))  , newpage = TRUE, heatmap_legend_side = "bottom",
                     annotation_legend_side = "bottom", merge_legend = TRUE)

ComplexHeatmap::draw(Heatmap(cov(dat), top_annotation = column_ha, name = "Cov", column_title = "Covariance pattern of the simulated the data",
                             col = colorRamp2(c(min(cov(dat)), 0, max(cov(dat))), c("blue", "white", "red")), right_annotation = row_ha,
                             heatmap_legend_param = list(direction = "horizontal"))  , newpage = TRUE, heatmap_legend_side = "bottom",
                     annotation_legend_side = "bottom", merge_legend = TRUE)


# ---- Premier Essai : ----

gl1 = glassoFast(S = cov(dat), rho = matrix(0.1, ncol = ncol(cov(dat)), nrow = nrow(cov(dat))))
ComplexHeatmap::draw(Heatmap(gl1$w, top_annotation = column_ha, name = "Cov estimation", column_title = "Covariance matrix found by glasso",
                             col = colorRamp2(c(min(gl1$w), 0, max(gl1$w)), c("blue", "white", "red")), right_annotation = row_ha,
                             heatmap_legend_param = list(direction = "horizontal"))  , newpage = TRUE, heatmap_legend_side = "bottom",
                     annotation_legend_side = "bottom", merge_legend = TRUE)


# Gremlin 1 table
Net = defineNetwork(gl1$w, typeInter = "diradj", rowFG = "Features", colFG = "Features")
# Net2 = defineNetwork(gl1$w, typeInter = "diradj", rowFG = "Feat", colFG = "Feat")
A_GLASSO <- ifelse(gl1$w!=0 & row(gl1$w)!=col(gl1$w),1,0)
NetA = defineNetwork(A_GLASSO, typeInter = "diradj", rowFG = "Features", colFG = "Features")

gr1 = multipartiteBM(list_Net = list(Net), namesFG = "Features",
                      v_distrib = "gaussian", v_Kmin = 2, v_Kmax = 3)

gr1A = multipartiteBM(list_Net = list(NetA), namesFG = "Features",
                     v_distrib = "bernoulli", v_Kmin = 2, v_Kmax = 3)

# gr1 = multipartiteBMFixedModel(list_Net = list(Net), namesFG = c("Features"),
                               # v_distrib = c("gaussian"), v_K = c(2))
# gr1 = multipartiteBMFixedModel(list_Net = list(Net, Net2), namesFG = c("Features", "Feat"),
                     # v_distrib = c("gaussian", "gaussian"), v_K = c(2, 2)) # N'a pas l'air de fonctionner quand on met qu'une seule matrice -> prévenir Sophie.

# clustr = gr1$fittedModel[[1]]$paramEstim$Z$Features
clustr = gr1$fittedModel[[1]]$paramEstim$Z[[1]]
clustrA = gr1A$fittedModel[[1]]$paramEstim$Z[[1]]
NID(clustr, classification)
NID(clustrA, classification) # au moins, le gaussien apporte de l'info... Dans ce cas là en tout cas. 

distrib_groups = gr1$fittedModel[[1]]$paramEstim$list_theta$FeaturesFeatures
set.seed(1992)
link1 = matrix(pnorm(0, as.numeric(distrib_groups$mean), sqrt(as.numeric(distrib_groups$var)), lower.tail = FALSE), ncol = ncol(distrib_groups$mean), nrow = nrow(distrib_groups$mean))
# proba que la proba de connexion tiree dans la loi soit au-dessus de 0 ? est-ce que c'est bon ? -> à demander.  (Moyenne et var peuvent yield des res negatifs ?)
link = link1[clustr, clustr] # doit servir pour la penalite. Tous les membres du meme groupe doivent avoir les memes proba de connexions/penalites

ComplexHeatmap::draw(Heatmap(link, top_annotation = column_ha, name = "Link", column_title = "Link between groups/individuals",
                             col = colorRamp2(c( 0, max(link)), c("white", "red")), right_annotation = row_ha,
                             heatmap_legend_param = list(direction = "horizontal"))  , newpage = TRUE, heatmap_legend_side = "bottom",
                     annotation_legend_side = "bottom", merge_legend = TRUE)

# glasso 2

gl2 = glassoFast(S = cov(dat), rho = 1-link, start = "warm", wi.init = gl1$wi, w.init = gl1$w)
ComplexHeatmap::draw(Heatmap(gl2$w, top_annotation = column_ha, name = "Cov estimation", column_title = "Covariance matrix found by glasso",
                             col = colorRamp2(c(min(gl2$w), 0, max(gl2$w)), c("blue", "white", "red")), right_annotation = row_ha,
                             heatmap_legend_param = list(direction = "horizontal"))  , newpage = TRUE, heatmap_legend_side = "bottom",
                     annotation_legend_side = "bottom", merge_legend = TRUE)
sum(abs(gl1$w - gl2$w)) # ok, quand meme des differences.
sum(abs(gl1$w - Sigma_sim))
sum(abs(gl2$w - Sigma_sim)) # Le besoin d'itérer est bien là. difference bien moins grande apres avoir passé le gremlin. 

# On part d'une matrice de penalisation ou tout le monde est penalise pareil.
# On ajuste par le GREMLIN/BM.
# On refait un Glasso avec une matrice de penalisation ajustee selon les groupes que l'on a trouve. 
# Mais si les groupes trouves sont "mauvais" ? c'est le cas ici. 
 
# ---- Tests sur boucles : ----

n_iter = 20
lambda_seq = seq(0.01, 0.2, length.out = 20)
# lambda_seq = seq(0.01, 0.5, length.out = 2)
listLambda = list()

for(lam in lambda_seq){
  cat("** Lambda : ", lam, "\n")
  
  list_res_glasso = list()
  list_res_gremlin = list()
  list_clusters = list()
  list_NID = list()
  list_diff_cov_sim = list()
  list_diff_iter_prec = list()
  
  gl0 = gl1 = gr0 = NULL
  
  # Iteration 0 : initialisation
  gl0 = glassoFast(S = cov(dat), rho = matrix(lam, ncol = ncol(cov(dat)), nrow = nrow(cov(dat))))
  list_res_glasso[[length(list_res_glasso)+1]] = gl0
  list_diff_cov_sim[[length(list_diff_cov_sim)+1]] = sum(abs(list_res_glasso[[length(list_res_glasso)]]$w - Sigma_sim))
  
  Net = defineNetwork(gl0$w, typeInter = "diradj", rowFG = "Features", colFG = "Features")
  gr0 = multipartiteBM(list_Net = list(Net), namesFG = "Features",
                       v_distrib = "gaussian", v_Kmin = 2, v_Kmax = 3, verbose = FALSE)
  
  list_res_gremlin[[length(list_res_gremlin)+1]] = gr0
  
  clustr = gr0$fittedModel[[1]]$paramEstim$Z[[1]]
  
  list_clusters[[length(list_clusters)+1]] = length(unique(clustr))
  list_NID[[length(list_NID)+1]] = NID(clustr, classification)
  
  distrib_groups = gr0$fittedModel[[1]]$paramEstim$list_theta$FeaturesFeatures
  set.seed(1992)
  link1 = matrix(pnorm(0, as.numeric(distrib_groups$mean), sqrt(as.numeric(distrib_groups$var)), lower.tail = FALSE), ncol = ncol(distrib_groups$mean), nrow = nrow(distrib_groups$mean))
  link = link1[clustr, clustr]
  
  gl1 = glassoFast(S = cov(dat), rho = 1-link, start = "warm", wi.init = gl0$wi, w.init = gl0$w)
  list_res_glasso[[length(list_res_glasso)+1]] = gl1
  list_diff_cov_sim[[length(list_diff_cov_sim)+1]] = sum(abs(list_res_glasso[[length(list_res_glasso)]]$w - Sigma_sim))
  list_diff_iter_prec[[length(list_diff_iter_prec)+1]] = sum(abs(list_res_glasso[[length(list_res_glasso)]]$w - list_res_glasso[[length(list_res_glasso) - 1 ]]$w))
  
  for(i in 1:n_iter){
    cat("Iter ", i, "\n")
    Net = defineNetwork(list_res_glasso[[length(list_res_glasso)]]$w , typeInter = "diradj", rowFG = "Features", colFG = "Features")
    gr = multipartiteBM(list_Net = list(Net), namesFG = "Features",
                         v_distrib = "gaussian", v_Kmin = 2, v_Kmax = 3, verbose = FALSE)
    list_res_gremlin[[length(list_res_gremlin)+1]] = gr
    
    clustr = gr$fittedModel[[1]]$paramEstim$Z[[1]]
    
    list_clusters[[length(list_clusters)+1]] = length(unique(clustr))
    list_NID[[length(list_NID)+1]] = NID(clustr, classification)
    
    distrib_groups = gr$fittedModel[[1]]$paramEstim$list_theta$FeaturesFeatures
    set.seed(1992)
    link1 = matrix(pnorm(0, as.numeric(distrib_groups$mean), sqrt(as.numeric(distrib_groups$var)), lower.tail = FALSE), ncol = ncol(distrib_groups$mean), nrow = nrow(distrib_groups$mean))
    link = link1[clustr, clustr]
    
    gl = glassoFast(S = cov(dat), rho = 1-link, start = "warm", wi.init = list_res_glasso[[length(list_res_glasso)]]$wi, w.init = list_res_glasso[[length(list_res_glasso)]]$w)
    list_res_glasso[[length(list_res_glasso)+1]] = gl
    list_diff_cov_sim[[length(list_diff_cov_sim)+1]] = sum(abs(list_res_glasso[[length(list_res_glasso)]]$w - Sigma_sim))
    list_diff_iter_prec[[length(list_diff_iter_prec)+1]] = sum(abs(list_res_glasso[[length(list_res_glasso)]]$w - list_res_glasso[[length(list_res_glasso) - 1 ]]$w))
    
  }
  
  listLambda[[length(listLambda) +1 ]] = list(
  "glasso" = list_res_glasso,
  "gremlin" = list_res_gremlin ,
  "clusters" = list_clusters,
  "NID" = list_NID,
  "diff_sim" = list_diff_cov_sim,
  "diff_prev_iter" = list_diff_iter_prec
  )
}

# length(listLambda)

df_test = do.call("rbind", lapply(listLambda, FUN = function(liste){
  # liste = listLambda[[1]]
  t(do.call("cbind", list("Clusters" = c(NA, unlist(liste$clusters)), "NID" = c(NA, unlist(liste$NID)),
                          "diff_sim" = c(unlist(liste$diff_sim)), "diff_prev_iter" = c(NA, unlist(liste$diff_prev_iter)))))
}))
colnames(df_test) = c("GL0", "GRGL0", paste0("GRGL", 1:n_iter))
df_test = data.frame(df_test)
df_test$lambda = rep(lambda_seq[1:length(listLambda)], each = 4)
df_test$response = rep(c("Clusters", "NID", "Diff_Sim", "Diff_prev_iter"), time = nrow(df_test)/4)
# listLambda[[1]]$NID

library(reshape2)
library(ggplot2)

df_melt = melt(df_test, id.vars = c("lambda", "response"))
df_melt$lambda = factor(df_melt$lambda)

ggplot(data = df_melt, aes(x = variable, y = value, color = lambda, group = lambda)) +
  geom_point() + geom_line() +
  facet_grid(response~., scales = "free")

# On voit rien, variations trop petites après la boucle 0.


df_melt2 = df_melt[-which(df_melt$variable%in%c("GL0", "GRGL0")),]

ggplot(data = df_melt2, aes(x = variable, y = value, color = lambda, group = lambda)) +
  geom_point() + geom_line() +
  facet_grid(response~., scales = "free")
# On voit toujours rien. 


rm(list = ls())
setwd("/home/hulot/Documents/packages_R/GlassoGREMLIN/2019_12_09_tests_Audrey/GlassoBM/")

# 10/12/19 : test block models

library(blockmodels)
library(Matrix)
library(aricode)
library(simone)

library(ggplot2)

library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

dyn.load("Fonctions/inference_algos.so")
sapply(dir("Fonctions/")[grep(".R", dir("Fonctions/"))], FUN = function(file) source(paste0("Fonctions/", file)))

# --- Simulation donnees : ----

n_ind = 100
n_feat = 30
nb_groupes = 5 # nombre de groupes de features

set.seed(1992)
prop_groupe = 0.20 # proportion de features dans les groupes (sauf dernier groupe qui prend le reste)
classification = base::sample(1:nb_groupes, size = n_feat, replace = TRUE, c(rep(prop_groupe, nb_groupes-1), 1-prop_groupe*(nb_groupes-1))) 
table(classification) 

Sigma_sim = matrix(NA, n_feat, n_feat) # covariance structure
for(i in 1:(nb_groupes)){
  set.seed(i*2000)
  Sigma_sim[which(classification==i), which(classification == i)] = runif(n = length(Sigma_sim[which(classification==i),which(classification == i)]), min = 30, max = 60)
}
Sigma_sim[is.na(Sigma_sim)] = runif(n = length(Sigma_sim[is.na(Sigma_sim)]), min = 0, max = 5)
diag(Sigma_sim) = apply(Sigma_sim, 2, max) + 1
Sigma_sim = as.matrix(nearPD(Sigma_sim)$mat)

library(MASS)

dat = mvrnorm(n = n_ind, mu = c(1:n_feat), Sigma = Sigma_sim)
dat_c = scale(dat, center = TRUE, scale = FALSE) # scale = FALSE ?

column_ha = HeatmapAnnotation(Classification = as.factor(classification),
                              col = list("Classification" = c("1" = "blue", "2" = "azure3", "3" = "darkorchid3", "4" = "darksalmon", "5" = "darkolivegreen3")),
                              annotation_legend_param = list(nrow = 1)
                              )

row_ha = HeatmapAnnotation(Classification = as.factor(classification),
                           col = list("Classification" = c("1" = "blue", "2" = "azure3", "3" = "darkorchid3", "4" = "darksalmon", "5" = "darkolivegreen3")), 
                           which = "row",
                           show_annotation_name = FALSE, show_legend = FALSE,  annotation_legend_param = list(nrow = 1))

ComplexHeatmap::draw(Heatmap(solve(Sigma_sim), top_annotation = column_ha, name = "Cov^-1", column_title = "partial Correlation pattern used to simulate the data",
                             col = colorRamp2(c(min(solve(Sigma_sim)),0, max(solve(Sigma_sim))), c("blue", "white", "red")), right_annotation = row_ha,
                             heatmap_legend_param = list(direction = "horizontal"))  , newpage = TRUE, heatmap_legend_side = "bottom",
                     annotation_legend_side = "bottom", merge_legend = TRUE)

ComplexHeatmap::draw(Heatmap(Sigma_sim, top_annotation = column_ha, name = "Cov", column_title = "Covariance pattern used to simulate the data",
                             col = colorRamp2(c(min(Sigma_sim), 0, max(Sigma_sim)), c("blue", "white", "red")), right_annotation = row_ha,
                             heatmap_legend_param = list(direction = "horizontal"))  , newpage = TRUE, heatmap_legend_side = "bottom",
                     annotation_legend_side = "bottom", merge_legend = TRUE)

ComplexHeatmap::draw(Heatmap(cov(dat), top_annotation = column_ha, name = "Cov", column_title = "Covariance pattern of the simulated the data",
                             col = colorRamp2(c(min(cov(dat)), 0, max(cov(dat))), c("blue", "white", "red")), right_annotation = row_ha,
                             heatmap_legend_param = list(direction = "horizontal"))  , newpage = TRUE, heatmap_legend_side = "bottom",
                     annotation_legend_side = "bottom", merge_legend = TRUE)


### ---- Block Model : ----


CV_res1 = CVgraphical.lasso(X = dat, nlam = 50, lam.min.ratio = 0.001,
                            maxit = 100,
                           diagonal = FALSE, crit.cv = "BIC",
                           cores = 4,
                           trace = "none")
plot(CV_res1)

ComplexHeatmap::draw(Heatmap(CV_res1$Sigma, top_annotation = column_ha, name = "Cov", column_title = "Covariance structure with CVglasso1",
                             right_annotation = row_ha, 
                             col = colorRamp2(c(min(CV_res1$Sigma), 0, max(CV_res1$Sigma)), c("blue", "white", "red")),
                             heatmap_legend_param = list(direction = "horizontal"))  , newpage = TRUE, heatmap_legend_side = "bottom", 
                            annotation_legend_side = "bottom", merge_legend = TRUE)


BM_model = BM_gaussian(
  membership_type = "SBM_sym", 
  verbosity = 0,
  explore_min = 1, 
  explore_max = 10,
  exploration_factor = 1.5,
  autosave = '',
  adj = CV_res1$Sigma)

# str(BM_model)
BM_model$estimate()
which.max(BM_model$ICL)

# Recuperation des groupes
groups = apply(BM_model$memberships[[which.max(BM_model$ICL)]]$Z,1, which.max)
NID(groups, classification) # ok

BM_model$model_parameters[[which.max(BM_model$ICL)]]

mu_mod = BM_model$model_parameters[[which.max(BM_model$ICL)]]$mu
sigma2_mod = BM_model$model_parameters[[which.max(BM_model$ICL)]]$sigma2

# penalisation_groupe = 1/abs(mu_mod) # pour l'instant
penalisation_groupe = 1-abs(mu_mod)/max(abs(mu_mod)) # 12/12/19

mat_penalty = penalisation_groupe[groups, groups]

CV_res = CVgraphical.lasso(X = dat, nlam = 20, lam.min.ratio = 0.0001,
                           maxit= 100,
                           diagonal = FALSE, 
                           trace = "none",
                           V = mat_penalty, crit.cv = "BIC", cores = 4)
plot(CV_res)

## Affichage des matrices de résultats

ComplexHeatmap::draw(Heatmap(CV_res$Sigma, top_annotation = column_ha, name = "Cov", column_title = "Covariance structure of the simulated data",
                             col = colorRamp2(c(min(CV_res$Sigma), 0, max(CV_res$Sigma)), c("blue", "white", "red")), 
                             right_annotation = row_ha, 
                             heatmap_legend_param = list(direction = "horizontal"))  , newpage = TRUE, heatmap_legend_side = "bottom", 
                     annotation_legend_side = "bottom", merge_legend = TRUE)


ComplexHeatmap::draw(Heatmap(mat_penalty, top_annotation = column_ha, name = "Penalty", column_title = "Structure of the penalty matrix",
                             col = colorRamp2(c(0, median(mat_penalty), max(mat_penalty)), c("white", "yellow", "red")),
                             right_annotation = row_ha,
                             heatmap_legend_param = list(direction = "horizontal"))  , newpage = TRUE, heatmap_legend_side = "bottom",
                     annotation_legend_side = "bottom", merge_legend = TRUE)


# ---- Test : alorithme 1 ----

# Algorithme 1 : cross-validation : 
# A l'intérieur d'un fold, pour un lambda donné, alternance entre glasso et BM jusqu'à "convergence" ou n.iter. 
# Recuperation du meilleur lambda, et alternance glasso et BM sur jeu de données complet jusqu'à convergence ou n.iter.
# l'étape BM donne une matrice de pénalité issue de l'estimation des moyennes des groupes, qui sert à compléter la pénalisation dans le glasso.
# Output : Sigma estimée (et Omega estimée) ainsi que la matrice de penalisation estimée par le dernier BM. 

res1 = CVglassoBM(X = dat, nlam = 50, lam.min.ratio = 0.001,
                        diagonal = FALSE,
                        maxit = 100,
                        membership_type = "SBM_sym", 
                        verbosity = 0,
                        K = 5,
                        explore_min = 4, 
                        explore_max = 6,
                        exploration_factor = 1.5,
                        autosave = '',
                        cores = 4,
                        crit.cv = "BIC",
                        trace= "none",
                        n.iter = 25,
                        thre.iter = 10^-3
                        )

plot(CVglBM_res)

## Affichage des matrices de résultats
ComplexHeatmap::draw(Heatmap(res1$Sigma, top_annotation = column_ha, name = "Cov", column_title = "Covariance structure found by 1st version of the algorithm",
                             col = colorRamp2(c(min(res1$Sigma), 0, max(res1$Sigma)), c("blue", "white", "red")),
                     right_annotation = row_ha, 
                             heatmap_legend_param = list(direction = "horizontal"))  , newpage = TRUE, heatmap_legend_side = "bottom", 
                     annotation_legend_side = "bottom", merge_legend = TRUE)

ComplexHeatmap::draw(Heatmap(res1$Omega, top_annotation = column_ha, name = "omega", column_title = "Partial Correlation structure found by 1st version of the algorithm",
                             col = colorRamp2(c(min(res1$Omega), 0, max(res1$Omega)), c("blue", "white", "red")),
                             right_annotation = row_ha, 
                             heatmap_legend_param = list(direction = "horizontal"))  , newpage = TRUE, heatmap_legend_side = "bottom", 
                     annotation_legend_side = "bottom", merge_legend = TRUE)

ComplexHeatmap::draw(Heatmap(res1$V, top_annotation = column_ha, name = "Penalty", column_title = "Structure of the penalty matrix by 1st version of the algorithm",
                             col = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red")),
                             right_annotation = row_ha,
                             heatmap_legend_param = list(direction = "horizontal"))  , newpage = TRUE, heatmap_legend_side = "bottom",
                     annotation_legend_side = "bottom", merge_legend = TRUE)

# ---- Test : alorithme 2 ----

# Algorithme 2 : 
# Alternance entre cross-validation pour choix du lambda Glasso et BM jusqu'à convergence ou n.iter.
# L'étape BM sert à estimer une matrice de pénalisation issue des moyennes des groupes trouvés, qui sert à compléter l'estimation de Sigma par la cross-validation pour le glasso.
# Output : Dernière estimation de sigma (et Omega) ainsi que dernière estimation de la matrice de penalité.

res2 = glBM_CV(X = dat, nlam = 50, lam.min.ratio = 0.001,
        diagonal = FALSE,
        maxit = 100,
        membership_type = "SBM_sym", 
        verbosity = 0,
        K = 5,
        explore_min = 4, 
        explore_max = 6,
        exploration_factor = 1.5,
        autosave = '',
        cores = 4,
        crit.cv = "BIC",
        trace= "none",
        n.iter = 25,
        thre.iter = 10^-3)

## Affichage des matrices de résultats
ComplexHeatmap::draw(Heatmap(res2$GLASSO$Sigma, top_annotation = column_ha, name = "Cov", column_title = "Covariance structure found by 2nd version of the algorithm",
                             col = colorRamp2(c(min(res2$GLASSO$Sigma), 0, max(res2$GLASSO$Sigma)), c("blue", "white", "red")),  right_annotation = row_ha, 
                             heatmap_legend_param = list(direction = "horizontal"))  , newpage = TRUE, heatmap_legend_side = "bottom", 
                     annotation_legend_side = "bottom", merge_legend = TRUE)

ComplexHeatmap::draw(Heatmap(res2$BM_step$mat_penalty, top_annotation = column_ha, name = "Penalty", column_title = "Structure of the penalty matrix found by 2nd version of the algorithm",
                             col = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red")),
                             right_annotation = row_ha,
                             heatmap_legend_param = list(direction = "horizontal"))  , newpage = TRUE, heatmap_legend_side = "bottom",
                     annotation_legend_side = "bottom", merge_legend = TRUE)


save.image("Session_final.RData")
load("Session_final.RData")


# ---- Differences matrices : -----
v = list(Sigma_sim = Sigma_sim, Sigma_dat = cov(dat), "res1" = res1$Sigma, "res2" = res2$GLASSO$Sigma)
z = outer(v,v, Vectorize(function(x,y) sum(abs(x-y)))) 
z

sum(abs(res2$BM_step$mat_penalty - res1$V)) # matrices dont les entrees sont entre 0 et 1 chacune

ComplexHeatmap::draw(Heatmap(abs(res2$GLASSO$Sigma - res1$Sigma), top_annotation = column_ha, name = "Diff",
                                 column_title = "Comparison covariance structure between algorithm",
                             col = colorRamp2(c(0, max(abs(res2$GLASSO$Sigma - res1$Sigma))), c("white", "red")),
                             right_annotation = row_ha,
                             heatmap_legend_param = list(direction = "horizontal"))  , newpage = TRUE, heatmap_legend_side = "bottom", 
                     annotation_legend_side = "bottom", merge_legend = TRUE)

# ---- Penalisation diagonale ? : ----

res1_diag = CVglassoBM(X = dat, nlam = 50, lam.min.ratio = 0.001,
                        diagonal = TRUE,
                        maxit = 100,
                        membership_type = "SBM_sym", 
                        verbosity = 0,
                        K = 5,
                        explore_min = 4, 
                        explore_max = 6,
                        exploration_factor = 1.5,
                        autosave = '',
                        cores = 4,
                        crit.cv = "BIC",
                        trace= "none",
                        n.iter = 25,
                        thre.iter = 10^-3
)

plot(res1_diag)
sum(abs(res1_diag$V  - CVglBM_res$V)) # Visiblement ne change rien de mettre diag = TRUE, la penalisation trouvee est la meme ?
# aucune option dans le grpahical.lasso puor el faire, mettre la diagonale ne sert en fait à rien.

## Affichage des matrices de résultats
ComplexHeatmap::draw(Heatmap(res1_diag$Sigma, top_annotation = column_ha, name = "Cov", column_title = "Covariance structure 1st version of the algorithm diag = TRUE",
                             col = colorRamp2(c(min(res1_diag$Sigma), 0, max(res1_diag$Sigma)), c("blue", "white", "red")),
                             right_annotation = row_ha, 
                             heatmap_legend_param = list(direction = "horizontal"))  , newpage = TRUE, heatmap_legend_side = "bottom", 
                     annotation_legend_side = "bottom", merge_legend = TRUE)

ComplexHeatmap::draw(Heatmap(res1_diag$Omega, top_annotation = column_ha, name = "omega", column_title = "Partial Correlation structure 1st version of the algorithm diag = TRUE",
                             col = colorRamp2(c(min(res1_diag$Omega), 0, max(res1_diag$Omega)), c("blue", "white", "red")),
                             right_annotation = row_ha, 
                             heatmap_legend_param = list(direction = "horizontal"))  , newpage = TRUE, heatmap_legend_side = "bottom", 
                     annotation_legend_side = "bottom", merge_legend = TRUE)

ComplexHeatmap::draw(Heatmap(res1_diag$V, top_annotation = column_ha, name = "Penalty", column_title = "Structure of the penalty matrix 1st version of the algorithm diag = TRUE",
                             col = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red")),
                             right_annotation = row_ha,
                             heatmap_legend_param = list(direction = "horizontal"))  , newpage = TRUE, heatmap_legend_side = "bottom",
                     annotation_legend_side = "bottom", merge_legend = TRUE)

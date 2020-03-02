rm(list = ls())
setwd("/home/hulot/Documents/packages_R/GlassoGREMLIN/2020_01_08/Code_Sophie/")

# ---- Packages ----
library(Matrix)
library(GREMLIN)
library(aricode)
library(glassoFast)
library(ggplot2)
library(ggpubr)
library(MASS)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
library(cluster)
source('functions_plot.R')
# ---- Simulation données : ----

n_ind = 10000
nb_groupes = 3 # nombre de groupes de features
n_feat = 30 # Nombre de variables

# Situation 1 : groupes de tailles comparables mais differentes.
#set.seed(1992)
prop_groupe = 0.33
classification = base::sample(1:nb_groupes, size = n_feat, replace = TRUE, c(rep(prop_groupe, nb_groupes-1), 1-prop_groupe*(nb_groupes-1)))
table(classification)

# classification = rep(1:nb_groupes, each = 10) ; table(classification) # Situation 2 : groupes de tailles equilibrees.

Sigma_sim = matrix(NA, n_feat, n_feat) # covariance structure
for(i in 1:(nb_groupes)){
  #set.seed(i*2000)
  Sigma_sim[which(classification==i), which(classification == i)] = runif(n = length(Sigma_sim[which(classification==i),which(classification == i)]), min =-1 , max = 1)
}
Sigma_sim[is.na(Sigma_sim)] = 0
Chol = matrix(0, ncol = ncol(Sigma_sim), nrow = ncol(Sigma_sim))
Chol[upper.tri(Chol, diag = TRUE)] = Sigma_sim[upper.tri(Sigma_sim, diag = TRUE)]
# Garder des  sur la diagonale lors de la multiplication
Chol = apply(Chol, 1, FUN = function(vect) if(sum(vect == 0) != length(vect)){sqrt(vect^2 / sum(vect^2))}else{rep(0, length(vect))})
Sigma_sim = t(Chol)%*%Chol

# dat = mvrnorm(n = n_ind, mu = c(1:n_feat)*classification, Sigma = Sigma_sim)
dat = mvrnorm(n = n_ind, mu =rep(0,n_feat), Sigma = Sigma_sim)

# ---- Visualisation matrice covariance : -----

Heatmap(Sigma_sim, col = colorRamp2(c(0, 0.5, max(Sigma_sim)), c("white", "yellow", "red")), name = "Cov_sim")
Heatmap(cov(dat), name = "Cov_dat", col = colorRamp2(c(-1, -0.5, 0, 0.5, max(cov(dat))), c("blue", "lightgray", "white", "yellow", "red")))

# ---- Recherche de groupes, nombre de groupe fixé: -----

Net = defineNetwork(cov(dat), typeInter = "adj", rowFG = "Features", colFG = "Features")
gr1 = multipartiteBMFixedModel(list_Net = list(Net), namesFG = c("Features"), v_distrib = c("gaussian"), v_K = 3)
clustr1 = gr1$fittedModel[[1]]$paramEstim$Z[[1]]

ICL1 <- gr1$fittedModel[[1]]$ICL

table(clustr1,classification)

gr1bis = multipartiteBMFixedModel(list_Net = list(Net), namesFG = "Features", classifInit = list(classification),
  v_distrib = "gaussian")
clustr1bis = gr1bis$fittedModel[[1]]$paramEstim$Z[[1]]
ICL1bis <- gr1bis$fittedModel[[1]]$ICL


table(clustr1bis,classification)
table(clustr1,clustr1bis)

# ---- Recherche de groupes et nombre de groupes : ----

gr2 = multipartiteBM(list_Net = list(Net), namesFG = "Features",
                     v_distrib = "gaussian", v_Kmin = 1, v_Kmax = 10, initBM = FALSE)

clustr2 = gr2$fittedModel[[1]]$paramEstim$Z[[1]]
ICL2 <- gr2$fittedModel[[1]]$ICL


#Heatmap(Sigma_sim[Zo,Zo], col = colorRamp2(c(0, 0.5, max(Sigma_sim)), c("white", "yellow", "red")), name = "Cov_sim")
###############################

gr3 <- BM_gaussian("SBM_sym",cov(dat) )
gr3$estimate()
Kbest <- which.max(gr3$ICL)
clustr3 <- max.col(gr3$memberships[[Kbest]]$Z)
print("=================ICL================")
print(c(ICL1,ICL1bis,ICL2))

par(mfrow = c(2,2))
graw <- plotMatrix(cov(dat),rowFG='Features',colFG='Features')
gtrue <- plotMatrix(cov(dat),rowFG='Features',colFG='Features',clustering = list(row = classification,col = classification))
g1 <- plotMatrix(cov(dat),rowFG='Features',colFG='Features',clustering = list(row=clustr1,col=clustr1))
g1bis <- plotMatrix(cov(dat),rowFG='Features',colFG='Features',clustering = list(row=clustr1bis,col=clustr1bis))
g2 <- plotMatrix(cov(dat),rowFG='Features',colFG='Features',clustering = list(row=clustr2,col=clustr2))
g3 <- plotMatrix(cov(dat),rowFG='Features',colFG='Features',clustering = list(row=clustr3,col=clustr3))


labels <- as.data.frame(c("A","B","C","D", "E","F"))
labels$metho <- c("Brut",'True','K fixe','K fixe init true classif','K estim',"blockModel Same var")
names(labels) <- c("fig","metho")
labels$ARI = c(NA,NA,ARI(clustr1,classification),ARI(clustr1bis,classification),ARI(clustr2,classification),ARI(clustr3,classification))
t.p <- ggtexttable(labels, rows = NULL,
  theme = ttheme("mOrange"))
figure <- ggarrange(graw, gtrue, g1,g1bis, g2,g3,t.p, labels = c("A","B","C","D", "E","F", "") ,nrow = 4, ncol = 2)
figure




rm(list = ls())
setwd("/home/hulot/Documents/packages_R/GlassoGREMLIN/2020_01_08/")

## ---- PACKAGES : -----

library(Matrix)
library(GREMLIN)
library(aricode)
library(glassoFast)
library(blockmodels)
library(simone)

library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggpubr)

library(MASS)

library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

library(igraph)
library(ggraph)

library(parallel)

# ---- SOURCE FONCTIONS : -----

source("Code_Sophie/functions_plot.R")
source("Function_graph.R")
source("Functions_boucles.R")

# ---- Function res graph plot : ----

res_sim  = function(res, n_iter){
  df_test = do.call("rbind",lapply(1:length(res), FUN = function(i){
    # i = 1 
    liste = res[[i]]
    df = data.frame(t(do.call("cbind", list("Clusters" = c(NA, unlist(liste$clusters), rep(NA, n_iter-length(unlist(liste$clusters)))),
                                            "NID" = c(NA, unlist(liste$NID), rep(NA, n_iter-length(unlist(liste$NID)))),
                                            "diff_sim" = c(unlist(liste$diff_sim), rep(NA, 1 + n_iter-length(unlist(liste$diff_sim)))),
                                            "diff_prev_iter" = c(NA, unlist(liste$diff_prev_iter), rep(NA, n_iter-length(liste$diff_prev_iter)))))))
    colnames(df) = c("GL0", "GRGL0", paste0("GRGL", 1:(ncol(df)-2)))
    df$response = rownames(df)
    df$lambda = lambda_seq[i] 
    df
  }))
  
  df_melt = melt(df_test, id.vars = c("lambda", "response"))
  df_melt$lambda = factor(df_melt$lambda)
  df_melt = df_melt[!is.na(df_melt$value),]
  df_melt$lambda = factor(round(as.numeric(as.character(df_melt$lambda)),2))
  
  df_melt
}


# ---- Fonctions simone: ----

# Laplacian, reprise simone
laplacian = function (M) 
{
  n <- nrow(M)
  diag(M) <- 1
  D <- colSums(M)
  L <- diag(rep(1, n)) - diag(D^(-1/2)) %*% M %*% diag(D^(-1/2))
  L[is.nan(L)] <- 0
  return(L)
}

# Changement pour signed
r_networkb =   function (p, pi, alpha = c(1), directed = FALSE, name = "a network",
                      signed = TRUE)
{
  A <- matrix(0, p, p)
  Q <- length(alpha)
  pi <- as.matrix(pi)
  fixed.edges <- FALSE
  if (max(pi) > 1) {
    if (sum(pi)%%floor(sum(pi)) == 0) {
      fixed.edges <- TRUE
    }
    else {
      stop("\nPi does not sum to one nor to an integer...")
    }
  }
  if (sum(alpha) < 1) {
    cat("\nNormalizing the vector of prior proportions which does not sum up to one\n")
    alpha <- alpha/sum(alpha)
  }
  if (!directed & !isSymmetric(pi)) {
    stop("\nPi should be symmetric for undirected graph...")
  }
  if (ncol(pi) != Q) {
    stop("\nPi and alpha do not have consistent dimensions...\n")
  }
  if (fixed.edges) {
    if (!directed) {
      pi.bounds <- floor(1/2 * (p * alpha) %*% t(p * alpha) -
                           diag(p * alpha/2, nrow = Q))
    }
    else {
      pi.bounds <- floor((p * alpha) %*% t(p * alpha))
    }
    if (any(pi > pi.bounds)) {
      cat("\nMatrix Pi is too dense: here is the adjusted version\n")
      pi[pi > pi.bounds] <- pi.bounds[pi > pi.bounds]
      print(pi)
    }
  }
  cl.levels <- factor(1:ncol(pi))
  if (is.null(rownames(pi))) {
    colnames(pi) <- rownames(pi) <- 1:ncol(pi)
  }
  names(alpha) <- cl.levels
  if (!fixed.edges) {
    clusters <- factor(sample(cl.levels, p, prob = alpha,
                              replace = TRUE), levels = cl.levels)
  }
  else {
    cl.prop <- round(p * alpha)
    clusters <- factor(sample(rep(cl.levels, cl.prop)), levels = cl.levels)
  }
  couple <- expand.grid(clusters, clusters)
  if (directed) {
    possible.edges <- rep(TRUE, nrow(couple))
  }
  else {
    up <- upper.tri(A, diag = FALSE)
    possible.edges <- rep(TRUE, nrow(couple)) & up
  }
  for (q in 1:Q) {
    for (l in 1:Q) {
      ql <- which(couple[, 1] == q & couple[, 2] == l &
                    possible.edges)
      if (!fixed.edges) {
        A[ql] <- rbinom(length(ql), 1, pi[q, l])
      }
      else {
        pi[q, l] <- min(pi[q, l], length(ql))
        if (pi[q, l] > 0) {
          A[sample(ql, pi[q, l])] <- 1
        }
      }
    }
  }
  if (!directed) {
    A <- (A | t(A)) * 1
  }
  dimnames(A) <- list(clusters, clusters)
  if (signed) {
    switchsign <- matrix(sample(c(1, -1), p * p, replace = TRUE),
                         p, p)
  }
  else {
    switchsign <- matrix(1, p, p) ## Modif here
  }
  if (!directed) {
    L <- laplacian(A)
    Theta <- -switchsign * t(switchsign) * L
    diag(Theta) <- rep(1, p)
  }
  else {
    Theta <- A
    corr <- runif(p * p)
    Theta[Theta != 0] <- corr[rank(corr) > sum(A == 0)] *
      ((switchsign * A)[A != 0])
    if (!(max(Mod(eigen(Theta)$values)) < 1)) {
      Theta <- Theta/(max(Mod(eigen(Theta)$values)) + 1e-09)
    }
  }
  nodes <- as.character(paste("g", 1:p, sep = ""))
  dimnames(Theta) <- list(nodes, nodes)
  return(structure(list(A = A, Theta = Theta, directed = directed,
                        clusters = clusters, name = name), class = "simone.network"))
}

# ---- TESTS RNETWORK : -----
p = 30
nb_groupes = 3
pi = matrix(0, nb_groupes, nb_groupes)
diag(pi) = c(1, 1, 1)
alpha = rep(1/nb_groupes, nb_groupes)

set.seed(1992)
net_1 = rNetwork(p = p,
                   pi = pi,
                   alpha    = alpha,
                   directed = FALSE, 
                   name     = "net_1")

# signed = FALSE ne fonctionne pas. -switchsign * t(switchsign) * L : tableaux de tailles inadéquates

# Theta contient les parametre du modele gaussien (variance si directed, ou Omega si undirected)
# A matrice d'adjacence, pas utile dans notre cas.
Omega = net_1$Theta
Sigma = solve(net_1$Theta)
Sigma_chol = chol2inv(net_1$Theta) # verif que ca marche, mais laplacien = ok de toute facon
classification = net_1$clusters

col2 = list("Cl" = c("1" = "blue", "2" = "azure3", "3" = "darkorchid3", "4" = "darksalmon", "5" = "darkolivegreen3"))
ha = HeatmapAnnotation(df = data.frame("Cl" = as.factor(classification)), col = col2)
ha_r = HeatmapAnnotation(df = data.frame("Cl" = as.factor(classification)), which = "row", col = col2)

draw(Heatmap(Omega - diag(diag(Omega)), top_annotation = ha, right_annotation = ha_r, name = "Omega"), merge_legend = TRUE)
draw(Heatmap(Sigma - diag(diag(Sigma)), top_annotation = ha, right_annotation = ha_r, name = "Sigma"), merge_legend = TRUE)

draw(Heatmap(Omega - diag(diag(Omega)), top_annotation = ha, right_annotation = ha_r, name = "Omega",
            row_order = order(classification), column_order = order(classification)), merge_legend = TRUE)
draw(Heatmap(Sigma - diag(diag(Sigma)), top_annotation = ha, right_annotation = ha_r, name = "Sigma",
             row_order = order(classification), column_order = order(classification)), merge_legend = TRUE)

Graphe = graph_from_adjacency_matrix(Omega, mode = "directed", weighted = TRUE, diag = FALSE)
Graphe = graph_plot(Graphe, classification)
ggraph(create_layout(Graphe, layout = "fr") ) + 
  geom_edge_link(alpha = .25, aes(width = width, colour = color)) + 
  scale_edge_width(range = c(0, 1)) +
  scale_edge_size_manual(values = 0.5) +
  geom_node_point(aes(colour = factor(classification)), size = 5) +
  theme_graph() + ggtitle("Graphe Omega simulée")


#### Est-ce que le fait que l'on ait des positifs et des négatifs gêne GREMLIN ? 
gr = multipartiteBM(list_Net = list(defineNetwork(Sigma, typeInter = "adj", rowFG = "Features", colFG = "Features")), 
                    v_distrib = "gaussian", namesFG = c("Features"), v_Kmin = 2, v_Kmax = 10)

NID(gr$fittedModel[[1]]$paramEstim$Z$Features, classification)
# Apparemment pas. Nid à 0

#### Simulation des données de transcriptomique associées à Omega/net_1
dat = rTranscriptData(1000, net_1)
# draw(Heatmap(dat$X, name = "data_sim"), merge_legend = TRUE)

draw(Heatmap(cov(dat$X) - diag(diag(cov(dat$X))), top_annotation = ha, right_annotation = ha_r, name = "Sigma_sim",
             row_order = order(classification), column_order = order(classification)), merge_legend = TRUE)

#### Essaie de boucle :
n_iter = 10
lambda_seq = seq(10^(-4), 10^(-1), length.out = 15)
res = loopGLGR(lambda_seq = lambda_seq, dat = dat$X, n_iter = n_iter, Omega_sim = Omega, classification = classification)

df_test = do.call("rbind",lapply(1:length(res), FUN = function(i){
  # i = 1 
  liste = res[[i]]
  df = data.frame(t(do.call("cbind", list("Clusters" = c(NA, unlist(liste$clusters), rep(NA, n_iter-length(unlist(liste$clusters)))),
                                          "NID" = c(NA, unlist(liste$NID), rep(NA, n_iter-length(unlist(liste$NID)))),
                                          "diff_sim" = c(unlist(liste$diff_sim), rep(NA, 1 + n_iter-length(unlist(liste$diff_sim)))),
                                          "diff_prev_iter" = c(NA, unlist(liste$diff_prev_iter), rep(NA, n_iter-length(liste$diff_prev_iter)))))))
  colnames(df) = c("GL0", "GRGL0", paste0("GRGL", 1:(ncol(df)-2)))
  df$response = rownames(df)
  df$lambda = lambda_seq[i] 
  df
}))

df_melt = melt(df_test, id.vars = c("lambda", "response"))
df_melt$lambda = factor(df_melt$lambda)
df_melt = df_melt[!is.na(df_melt$value),]
df_melt$lambda = factor(round(as.numeric(as.character(df_melt$lambda)),2))
# df_print = df_melt[as.numeric(as.character(df_melt$lambda))<10^(-4),]

ggplot(data = df_melt, aes(x = variable, y = value, color = lambda, group = lambda)) +
  geom_point(size = 3) + geom_line(size = 1) +
  facet_grid(response~lambda, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(legend.position="none", text = element_text(size = 20))

ggsave(filename = "Test_net1.jpeg", height = 20, width = 30, units = c("cm"))

# ---- Tests net_1 : données simulées avec des nombres différents d'individus : ----

n_ind_pos = c(10, 50, 100, 500, 1000)

dat_pos = lapply(n_ind_pos, rTranscriptData, graph = net_1)

n_iter = 10
lambda_seq = seq(10^(-4), 10^(-1), length.out = 10)

res = mclapply(dat_pos, FUN = function(dat){
  loopGLGR(lambda_seq = lambda_seq, dat = dat$X, n_iter = n_iter, Omega_sim = Omega, classification = classification)
}, mc.cores = 4)
# Probleme avec les convergences de GREMLIN. 

tab_res = lapply(res, res_sim, n_iter = n_iter)

lapply(1:length(tab_res), FUN = function(i){
  ggplot(data = tab_res[[1]], aes(x = variable, y = value, color = lambda, group = lambda)) +
    geom_point(size = 3) + geom_line(size = 1) +
    facet_grid(response~lambda, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    theme(legend.position="none", text = element_text(size = 20)) +
    ggtitle(paste0("N_ind" = n_ind_pos[i]))
})

# ---- Tests possibilities : ----

p = 30
nb_groupes = 3
alpha = rep(1/nb_groupes, nb_groupes)

pi_pos = seq(0.01, 1, length.out =  20)
pi_diag = expand.grid(pi_pos, pi_pos, pi_pos)
nrow(pi_diag)

net_res = mclapply(1:nrow(pi_diag), FUN = function(i){
# i = 8000
  pi = matrix(0, nb_groupes, nb_groupes)
  diag(pi) = as.numeric(pi_diag[i,])
  
  set.seed(1992)
  rNetwork(p = p, pi = pi,
                   alpha    = alpha,
                   directed = FALSE, 
                   name     = "net_1")
})

net_res_summary = mclapply(1:length(net_res), FUN = function(i){
  # i = 1
  Omega = net_res[[i]]$Theta
  Sigma = solve(Omega)
  
  Omega = Omega - diag(diag(Omega))
  Sigma = Sigma - diag(diag(Sigma))
  
  return(list("Omega_sum" = c("Min" = min(Omega), "Max" = max(Omega), "NB0" = sum(Omega == 0)),
              "Sigma_sum" = c("Min" = min(Sigma), "Max" = max(Sigma), "NB0" = sum(Sigma == 0))))
  
})

tab_omega_sum = do.call("rbind", lapply(net_res_summary, FUN = function(x) x$"Omega_sum"))
tab_sigma_sum = do.call("rbind", lapply(net_res_summary, FUN = function(x) x$"Sigma_sum"))

tab_omega_sum = data.frame("Pos" = 1:nrow(pi_diag), pi_diag, tab_omega_sum)
tab_sigma_sum = data.frame("Pos" = 1:nrow(pi_diag), pi_diag, tab_sigma_sum)

ggplot(data = tab_omega_sum, aes(x = Pos, y = Min)) + geom_point()
ggplot(data = tab_omega_sum, aes(x = Pos, y = NB0)) + geom_point()

# peu de cas où on arrive à des valeurs "correctes".
ggplot(data = tab_sigma_sum, aes(x = Pos, y = Min)) + geom_point()
# C'est pareil, c'est tout petit !

# ---- Tests breast cancer data TCGA : ----

#### Objectif : voir la répartition des corrélations partielles et des covariances.
####  Voir si on peut s'en servir pour créer une table simulée "crédible". 

load("Data_tests/tcga_brca_data.RData")
load("Data_tests/clinical.RData")

zrna2 = log2(zrna + 1)
zrna2 = zrna2[,order(apply(zrna2, 2, var), decreasing = TRUE)[1:500]]

listTab = list("mirna" = zmirna, "protein" = zprotein, "rnaseq" = zrna2)
listTab = lapply(listTab, FUN = function(tab) tab[order(rownames(tab)),])

listCov = lapply(listTab, cov)
listOmega = lapply(listCov, FUN = function(tabCov) as.matrix(solve(nearPD(tabCov)$mat)))

draw(Heatmap(listCov$mirna - diag(diag(listCov$mirna)), show_column_names = FALSE, show_row_names = FALSE, name = "mirna"))
draw(Heatmap(listCov$protein - diag(diag(listCov$protein)), show_column_names = FALSE, show_row_names = FALSE, name = "protein"))
draw(Heatmap(listCov$rnaseq - diag(diag(listCov$rnaseq)), show_column_names = FALSE, show_row_names = FALSE, name = "rnaseq"))

draw(Heatmap(listOmega$mirna - diag(diag(listOmega$mirna)), show_column_names = FALSE, show_row_names = FALSE, name = "mirna_Om"))
draw(Heatmap(listOmega$protein - diag(diag(listOmega$protein)), show_column_names = FALSE, show_row_names = FALSE, name = "protein_Om"))
draw(Heatmap(listOmega$rnaseq - diag(diag(listOmega$rnaseq)), show_column_names = FALSE, show_row_names = FALSE, name = "rnaseq_Om"))

dat_res = rbind(unlist(lapply(listCov, FUN = function(tabCov) c(quantile(tabCov, 0.10), quantile(tabCov, 0.9)))), 
      unlist(lapply(listOmega, FUN = function(tabCov) c(quantile(tabCov, 0.10), quantile(tabCov, 0.9)))))
rownames(dat_res) = c("Cov", "InvCov")

round(dat_res, 2)

par(mfrow = c(1,2))
hist(unlist(listCov$mirna), breaks = 50) ; hist(unlist(listOmega$mirna), breaks = 50)
hist(unlist(listCov$protein), breaks = 50) ; hist(unlist(listOmega$protein), breaks = 50)
hist(unlist(listCov$rnaseq), breaks = 50) ; hist(unlist(listOmega$rnaseq), breaks = 50)


# ---- Test boucle Glasso/GREmlin sur bcancer data : -----


prot_var = apply(listTab$protein, 2, var)
prot2 = listTab$protein[,order(prot_var)[1:50]]
n_iter = 25
lambda_seq = seq(10^(-5), 10^(-3), length.out = 15)
res = loopGLGR_nosim(lambda_seq = lambda_seq, dat = prot2, n_iter = n_iter)

df_test = do.call("rbind",lapply(1:length(res), FUN = function(i){
  # i = 1
  liste = res[[i]]
  df = data.frame(t(do.call("cbind", list("Clusters" = c(NA, unlist(liste$clusters), rep(NA, n_iter-length(unlist(liste$clusters)))),
                                          "diff_prev_iter" = c(NA, unlist(liste$diff_prev_iter), rep(NA, n_iter-length(liste$diff_prev_iter)))))))
  colnames(df) = c("GL0", "GRGL0", paste0("GRGL", 1:(ncol(df)-2)))
  df$response = rownames(df)
  df$lambda = lambda_seq[i] 
  df
}))

df_melt = melt(df_test, id.vars = c("lambda", "response"))
df_melt$lambda = factor(round(as.numeric(as.character(df_melt$lambda)), 6))

# df_print = df_melt[as.numeric(as.character(df_melt$lambda))<10^(-4),]

ggplot(data = df_melt, aes(x = variable, y = value, color = lambda, group = lambda)) +
  geom_point(size = 3) + geom_line(size = 1) +
  facet_grid(response~lambda, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(legend.position="none", text = element_text(size = 20))

setwd("/home/hulot/Documents/packages_R/GlassoGREMLIN/2020_03_20/")

load("res_loop.RData")

list_classif = lapply(1:26, FUN = function(i){
  res$`1e-04`$gremlin[[i]]$fittedModel[[1]]$paramEstim$Z$Features})

library(ggalluvial)
library(tidyverse)
library(reshape2)

res.alluv = function(list_classif, classif = NULL, title = ""){
  df_alluvial = data.frame(do.call(cbind, list_classif))
  colnames(df_alluvial) = c("GR0", paste0("GRGL", 1:25))
  if(!is.null(classif)) df_alluvial$classif = classif
  df_alluvial$Ind = 1:nrow(df_alluvial)
  
  
  df_melt = melt(df_alluvial, id.vars = "Ind")
  df_melt$value = factor(df_melt$value)
  colnames(df_melt) = c("Ind", "Iteration", "Groupe")
  
  ggplot(df_melt,
         aes(x = variable, stratum = value, alluvium = Ind,
             fill = value, label = value)) +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    geom_flow(stat = "alluvium", lode.guidance = "frontback",
              color = "darkgray") +
    geom_stratum() +
    theme(legend.position = "bottom")  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(title)
}




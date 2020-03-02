## Plot de graphe (igraph)

require(igraph)

graph_plot = function(Graphe, classification){
  # E(Graphe)$color[E(Graphe)$weight > 0] <- 'darkgreen'
  # E(Graphe)$color[E(Graphe)$weight < 0] <- 'red'
  E(Graphe)$color[E(Graphe)$weight > 0] <- 'positive'
  E(Graphe)$color[E(Graphe)$weight < 0] <- 'negative'
  E(Graphe)$width = abs(E(Graphe)$weight)
  E(Graphe)$width = E(Graphe)$width/max(E(Graphe)$width)
  E(Graphe)$arrow.mode = 0
  V(Graphe)$label = 1:30
  V(Graphe)$size = 3
  # V(Graphe)$label.color = classification
  V(Graphe)$color = classification
  return(Graphe)
}


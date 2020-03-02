require(ggplot2)
require(igraph)

################################ PLOT Matacency   Matrix
plotMatrix = function(Mat,rowFG,colFG, fileNameSave = NULL, clustering = NULL){

  n1 <- dim(Mat)[1]
  n2 <- dim(Mat)[2]
  u <- range(c(Mat))

  if (!is.null(clustering)){
    l <- length(clustering)
    if (l == 1){
      oRow <- oCol <- order(clustering$row)
      uRow <- cumsum(table(clustering$row)) + 0.5
      uRow <- uRow[-length(uRow)]
      sepRow <- as.data.frame(uRow)
      sepCol <- sepRow
    }
    if (l == 2){
      oRow <- order(clustering$row)
      oCol <- order(clustering$col)
      uRow <- cumsum(table(clustering$row)) + 0.5
      uRow <- uRow[-length(uRow)]
      sepRow <- as.data.frame(uRow)
      uCol <- cumsum(table(clustering$col)) + 0.5
      uCol <- uCol[-length(uCol)]
      sepCol <- as.data.frame(uCol)
    }
    Mat <- Mat[oRow,oCol]
    names(sepCol) = names(sepRow) = 'sep'


    sepRow = n1 - sepRow
  }

  index_row = rep(1:dim(Mat)[1],each = dim(Mat)[2])
  index_col = rep(1:dim(Mat)[2],dim(Mat)[1])


  melted_Mat =  data.frame(n1 - index_row , index_col)
  link = rep(-10,dim(Mat)[2]*dim(Mat)[1])
  for (k in 1:(dim(Mat)[2] * dim(Mat)[1])){ link[k] = Mat[index_row[k],index_col[k]]}
  melted_Mat$link = link
  colnames(melted_Mat) <- c('index_row', 'index_col', 'link')

  g <- ggplot(data = melted_Mat, aes(y=index_row, x=index_col, fill=link)) + geom_tile() + scale_fill_gradient(low="white", high="black", limits=u)
  g <- g + theme_bw() +  scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE)
  g <- g + theme(
    # Rotate the x-axis lables so they are legible
    axis.text.x = element_text(angle = 270, hjust = 0))
  # Force the plot into a square aspect ratio
  # Hide the legend (optional)
  # legend.position = "none")
  g <- g +  labs(x = colFG, y = rowFG)
  g <- g + theme(aspect.ratio = n1/n2)
  if (!is.null(clustering)){
    g <- g + geom_vline(data = sepCol,mapping=aes(xintercept=sep),col = 'grey')
    g <- g + geom_hline(data = sepRow,mapping=aes(yintercept=sep),col = 'grey')


  }
  if (!is.null(fileNameSave)) { ggsave(fileNameSave, width = 20, height = 20, units = "cm") }else{g}
  return(g)
}

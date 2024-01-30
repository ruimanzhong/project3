plotGraph <- function(coords, graph, title = NULL){
  edgelist = get.edgelist(graph, names = F) 
  edgedata = data.frame(coords[edgelist[, 1], ],coords[edgelist[, 2], ])
  colnames(edgedata) <- c("X1","Y1","X2","Y2")
  
  ggplot() + 
    geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2), data = edgedata, size = 0.5, colour = "grey")+
    labs(title = title, x = "lon", y = "lat")+
    theme(plot.title = element_text(hjust = 0.5))
}

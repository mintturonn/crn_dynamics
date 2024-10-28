
library(igraph)

# Convert the adjacency matrix to an igraph graph
g <- graph_from_adjacency_matrix(linkm, mode = "undirected", diag = FALSE)

# Plot the graph
plot(g, 
     vertex.label = NA, 
     vertex.size = 1, 
     edge.color = "gray",
     layout = layout_with_kk)
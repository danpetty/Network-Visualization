library(shiny)
library(igraph)
library(visNetwork)
library(pracma)
library(stringr)

graph_grouper <- function(graph) {
  for (i in 1:nrow(graph$nodes)){
  graph$nodes$group[i] <- trimws(substr(graph$nodes$label[i],1,3))
  }
  return(graph)
}
graph_color_adder <- function(graph, fire_matrix, inputVal){
  rbPal <- colorRampPalette(c('red','blue'))
  col <- rbPal(20)[as.numeric(cut(fire_matrix[,inputVal],breaks = 20))]
  for (i in 1:nrow(graph$nodes))
    {
      index <- as.numeric(str_extract(graph$nodes$id[i], pattern = "[[:digit:]]*$"))
      graph$nodes$color[i] <- col[index]
  
  }
  return(graph)
}

graph_dim_adder <- function(graph, coor) {
  for (i in 1:nrow(graph$nodes)){
    index <- as.numeric(str_extract(graph$nodes$id[i], pattern = "[[:digit:]]*$"))
    graph$nodes$x[i] <- coor[index,1]
    graph$nodes$y[i] <- coor[index,2]
  }
  return(graph)
}

betweenness_to_size <- function(graph,betweenness_list) {
  for (i in 1:nrow(graph$nodes)) {
    betweenness_list.rank <- rank(betweenness_list)/length(betweenness_list)
    graph$nodes$value[i] <- betweenness_list.rank[graph$node$id[i]]
  }
  return(graph)
}

isolated_node_remover <- function(graph){
  iso <- V(graph)[degree(graph) == 0]
  print(iso) #gets isolated nodes
  graph <- delete.vertices(graph, iso) #deletes them
  return(graph)
}

cell_type_to_node_shape <- function(graph, type_list){
  for (i in 1:nrow(graph$nodes)){
    index <- as.numeric(str_extract(graph$nodes$id[i], pattern = "[[:digit:]]*$"))
    if (type_list[index]=="+") {graph$nodes$shape[i] <- "triangle"}
    else if (type_list[index]=="-") {graph$nodes$shape[i] <- "circle"}

  }
  return(graph)
}

new_col_names <- NULL
cell_names <- read.csv("cell_names.csv", stringsAsFactors = FALSE)
new_group_names <- str_extract(cell_names,"^\\w*") #extracts the groups from the cell_names csv

for (i in 1:122){
  new_col_names[i] <- paste(new_group_names[i],i)
}
smallworld_col_names <- c(1:122)


#calculates the number of neurons in each group
numDG <- sum(new_group_names=="DG")
numEC <- sum(new_group_names=="EC")
numSUB <- sum(new_group_names=="SUB")
numCA1 <- sum(new_group_names=="CA1")
numCA2 <- sum(new_group_names=="CA2")
numCA3 <- sum(new_group_names=="CA3")

#creates some rough locations for the different neurons to appear
DG_xy_coor <- matrix(c(runif(numDG,25,35),linspace(25,35,numDG)),nrow=numDG,ncol=2)
EC_xy_coor <- matrix(c(runif(numEC,6,19),linspace(25,35,numEC)),nrow=numEC,ncol=2)
SUB_xy_coor <- matrix(c(runif(numSUB,20,23),linspace(25,28,numSUB)),nrow=numSUB,ncol=2)
CA1_xy_coor <- matrix(c(linspace(7,29,numCA1),runif(numCA1, 10,20)),nrow=numCA1,ncol=2)
CA2_xy_coor <- matrix(c(linspace(30,39,numCA2),runif(numCA2,13,24)),nrow=numCA2,ncol=2)
CA3_xy_coor <- matrix(c(runif(numCA3,36,43),linspace(25,35,numCA3)),nrow=numCA3,ncol=2)

#creates a full coordinate matrix
full_coor <- rbind(DG_xy_coor,CA3_xy_coor, CA2_xy_coor, CA1_xy_coor, SUB_xy_coor, EC_xy_coor)

firing_rate_0_125 <- read.csv(file="weight_0.125firing_of_rate.csv", header=FALSE)
firing_rate_0_3 <- read.csv(file="weight_0.3firing_of_rate.csv", header=FALSE)
firing_rate_1_0 <- read.csv(file="weight_1.0firing_of_rate.csv", header=FALSE)
firing_rate_1_5 <- read.csv(file="weight_1.5firing_of_rate.csv", header=FALSE)
firing_rate_matrix <- cbind(firing_rate_0_125,firing_rate_0_3,firing_rate_1_0,firing_rate_1_5)
colnames(firing_rate_matrix) <- c("Firing Rate 0.125","Firing Rate 0.3","Firing Rate 1.0","Firing Rate 1.5")

cell_type <- NULL
for(i in 1:length(cell_names)){
  grab_string <- str_extract(cell_names[i], pattern = "\\(.*?\\)")  #extracts inhibitory or excitatory cells
  cell_type[i] <- gsub("[()]", "", grab_string)
}

#imports ii matrix, gives it col/row names from above
ii_mat <- read.csv(file="ii.csv", header = TRUE, row.names=NULL, stringsAsFactors = FALSE, colClasses = c("logical"))
colnames(ii_mat) <- new_col_names
row.names(ii_mat) = c(colnames(ii_mat))
ii_mat <- ii_mat+0 #makes logical values into numeric

ii_graph <- graph_from_adjacency_matrix(as.matrix(ii_mat), mode = "directed")
#uses igraph function to create adj matrix

ii_graph <- isolated_node_remover(ii_graph)


ee_mat <- read.csv(file="ee.csv", header = TRUE, row.names=NULL,stringsAsFactors = FALSE, colClasses = c("logical"))
colnames(ee_mat) <- new_col_names
row.names(ee_mat) = c(colnames(ee_mat))
ee_mat <- ee_mat+0
ee_graph <- graph_from_adjacency_matrix(as.matrix(ee_mat), mode = "directed")
ee_graph <- isolated_node_remover(ee_graph)


ei_mat <- read.csv(file="ei.csv", header = TRUE, row.names=NULL,stringsAsFactors = FALSE, colClasses = c("logical"))
colnames(ei_mat) <- new_col_names
row.names(ei_mat) = c(colnames(ei_mat))
ei_mat <- ei_mat+0
ei_graph <- graph_from_adjacency_matrix(as.matrix(ei_mat), mode = "directed")
ei_graph <- isolated_node_remover(ei_graph)



ie_mat <- read.csv(file="ie.csv", header = FALSE, skip=1, row.names=NULL,stringsAsFactors = FALSE, colClasses = c("logical"))
colnames(ie_mat) <- new_col_names
row.names(ie_mat) = c(colnames(ie_mat))
ie_mat <- ie_mat+0
ie_graph <- graph_from_adjacency_matrix(as.matrix(ie_mat), mode = "directed")
ie_graph <- isolated_node_remover(ie_graph)


shinyServer(function(input, output) {
  
  graph_to_plot <- NULL
  output$network <- renderVisNetwork({
    if (input$select_matrix == "Inhibitory to Inhibitory") {
    graph_to_plot <- ii_graph
    }
    else if (input$select_matrix == "Excitatory to Excitatory") {
      graph_to_plot <- ee_graph
    }
    else if (input$select_matrix == "Excitatory to Inhibitory") {
      graph_to_plot <- ei_graph
      
    }
    else if (input$select_matrix == "Inhibitory to Excitatory") {
      graph_to_plot <- ie_graph
      
    }
    plot_betweenness <- betweenness(graph_to_plot)
    final_graph <- toVisNetworkData(graph_to_plot) #turns igraph into
    final_graph <- graph_grouper(final_graph)
    final_graph <- graph_dim_adder(final_graph, full_coor)
    final_graph <- cell_type_to_node_shape(final_graph, cell_type)
    if(input$color!='Group'){final_graph <- graph_color_adder(final_graph, firing_rate_matrix,input$color)}
    if (input$size == "Betweenness") {final_graph <- betweenness_to_size(final_graph, plot_betweenness)}
    visNetwork(nodes = final_graph$nodes, edges = final_graph$edges, height = "800px")  %>% 
      visIgraphLayout() %>% 
      visInteraction(dragNodes = FALSE, 
                     dragView = TRUE, 
                     zoomView = TRUE)%>%
      visLayout(randomSeed = 123) %>% 
      visOptions(highlightNearest = list(enabled =TRUE, degree = 1)) %>% 
      visInteraction(navigationButtons = TRUE) %>% 
      visNodes(font=list(size=25)) %>% 
      visEdges(arrows =list(to = list(enabled = TRUE, scaleFactor = 2))) %>% 
      visLegend()
    
    
      
                                                                                    
    
  })
  
  
})

library(shiny)
library(chorddiag)
library(igraph)
#library(network)
#library(sna)
#library(ggplot2)
#library(intergraph)
#library(GGally)
library(visNetwork)

shinyUI(fluidPage(
  fluidRow(
  radioButtons('select_matrix',"Select Network",inline = TRUE,
               choices = c("Excitatory to Excitatory","Inhibitory to Inhibitory","Inhibitory to Excitatory","Excitatory to Inhibitory"),
               selected = "Excitatory to Excitatory"),
  column(3,
  selectInput('size', "Select Size",
                choices = c("Uniform","Betweenness"), selected = "Uniform")
  ),
  column(5,
    selectInput('color', "Select Color",
               choices = c("Group","Firing Rate 0.125","Firing Rate 0.3","Firing Rate 1.0"), selected = "Group")
  ),
  mainPanel(visNetworkOutput("network"),width=12, column=12)
)  
)
)
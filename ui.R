library(shiny)
library(haven)
library(mice)
library(ggplot2)
library(gganimate)
library(shinycssloaders)
library(kableExtra)

# Define UI  ----
ui <- fluidPage(
  
  titlePanel("Multiple Imputation"),
  
  # Sidebarpanel
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
    
    # Input: Selector variable to plot  ----
    selectInput("Method", "Imputation Method", 
                c("Bayesian Linear Regression" = "norm")),
    
    numericInput("imp", "Imputations:", 1, min = 1, max = 1),
    numericInput("niter", "Iterations:", 20 , min = 20, max = 20),
    
    hr(),
    actionButton("gobutton","Show Imputed Values"),
    
    h4("Download original data (after convergence plot appeared)"),
    downloadLink("downloadData", "Download SPSS File"),
    
    tableOutput("Imputed_Mean_Chain")
    
    ),
  
  # Main panel for displaying outputs ----
  mainPanel(
    
    h3("Convergence plot .... may take a while"),
    withSpinner(imageOutput("plot1"), type = 5),
    
    tableOutput("table2")
    
  )
)
)

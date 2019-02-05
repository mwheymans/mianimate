library(shiny)
library(haven)
library(mice)
library(ggplot2)
library(gganimate)
library(shinycssloaders)

# Define UI  ----
ui <- fluidPage(
  
  titlePanel("Multiple Imputation"),
  
  # Sidebarpanel
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
    
    # Input: Selector variable to plot  ----
    selectInput("Method", "Imputation Method", 
                c("Mean" = "mean",
                  "Single Regression" = "norm.predict",
                  "Predictive Mean Matching" = "pmm")),
    
    numericInput("imp", "Imputations:", 1, min = 1, max = 100),
    numericInput("niter", "Iterations:", 5, min = 1, max = 100),
    
    tableOutput("table1")
    
    ),
  
  # Main panel for displaying outputs ----
  mainPanel(
    
    h3("Convergence plot"),
    withSpinner(imageOutput("plot1"), type = 5)
    
  )
)
)
library(shiny)
library(haven)
library(mice)
library(ggplot2)

# Define UI  ----
ui <- fluidPage(
  
  titlePanel("Simple Imputation methods"),
  
  # Sidebarpanel
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
    
    # Input: Selector variable to plot  ----
    selectInput("Method", "Imputation Method", 
                c("Full dataset (without missings)", "Complete Case Analysis", 
                  "Mean" = "mean",
                  "Single Regression" = "norm.predict",
                  "Single Stochastic Regression" = "norm",
                  "Predictive Mean Matching" = "pmm")),
    
    tableOutput("table1")
    
    ),
  
  # Main panel for displaying outputs ----
  mainPanel(
    
    h3("Imputed values for the Tampa scale variable"),
    plotOutput("plot2"),
    h3("Regression line in completed data"),
    plotOutput("plot1"),
    
   
    h3("Linear regression model in the Completed data"),
    verbatimTextOutput("model1")
    
  )
)
)
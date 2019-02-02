# Define server logic  ----
server <- function(input, output) {
  
  selectedData_comp <- read_sav(file="Backpain 50 missing_complete.sav")
  selectedData_comp$Group <- factor(rep(1, 50))
  levels(selectedData_comp$Group) <- c("Observed")
  
  selectedData1 <- reactive({
    data <- read_sav(file="Backpain 50 missing.sav")
    data$Group <- factor(ifelse(is.na(data$Tampascale), 1,0))
    
    levels(data$Group) <- c("Observed", "Imputed")
    imp0 <- mice(data=data, m=1, method="mean", maxit = 0)
    pred <- imp0$predictorMatrix
    pred[, c("ID", "Group", "Disability", "Radiation")] <- 0
    if(input$Method=="Full dataset (without missings)") {
      selectedData1 <- selectedData_comp 
    } else {
    if(input$Method=="Complete Case Analysis") {
      selectedData1 <- data 
    } else{
    imp <- mice(data=data, method=input$Method, predictorMatrix = pred, m=1, maxit = 1, seed=1232)
    complete(imp, 1)
      }
    }
    
    }) 
  
  Tampascale_orig <- read_sav(file="Backpain 50 missing_complete.sav")[1:20, "Tampascale"]
  
  output$plot1 <- renderPlot({ 
    
    g1 <- ggplot(data=selectedData1(), aes(x=Tampascale, y=Pain)) + 
      geom_point(size=3) + geom_smooth(method=lm) 
    plot(g1)
  })  
  
  output$plot2 <- renderPlot({ 
    
    g1 <- ggplot(data=selectedData1(), aes(x=Pain, y=Tampascale, color=factor(Group))) + 
      geom_point(size=3) + geom_smooth(method=lm, se=FALSE) 
    plot(g1)
  
    
    })  
  
  
  output$table1 <- renderTable({
    table1 <- data.frame(head(selectedData1()[, c("Pain", "Tampascale")],20), 
               Tampascale_orig)
    names(table1) <- c("Pain", "Tampascale", "Tampascale_mis")
    table1
  })
  
  output$model1 <- renderPrint({
      
      summary(lm(Pain ~ Tampascale, data=selectedData1()))
    })
  
  }

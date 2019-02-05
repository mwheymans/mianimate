# Define server logic  ----
server <- function(input, output) {
  
  selectedData1 <- reactive({
    
    dat <- read_sav(file="Backpain50 MI missing.sav")
    imp0 <- mice(data=dat, m=1, method="mean", maxit = 0)
    pred <- imp0$predictorMatrix
    pred[, c("ID", "Disability", "Radiation")] <- 0
    imp <- mice(data=dat, predictorMatrix = pred, method=input$Method, 
                m=input$imp, maxit = input$niter, seed=1232)
    
    imp_mean <- apply(imp$chainMean, 1L, c)[, "Tampascale"]
    iternr <- rep(1:input$niter, input$imp)
    impnr <- rep(1:input$imp, each=input$niter)
    dat_imp <- data.frame(impnr, iternr, imp_mean)
    dat_imp
    }) 
  
  output$table1 <- renderTable({
    table1 <- data.frame(selectedData1())
    table1
  })
  
  output$plot1 <- renderImage({ 
  outfile <- tempfile(fileext='.gif')
    
  p <- ggplot(selectedData1(), aes(iternr, imp_mean, group = impnr)) +
      geom_line(aes(colour=rep(c(1:input$imp), each=input$niter))) +
      geom_segment(aes(xend = input$niter, yend = imp_mean), linetype = 2, colour = 'grey') + 
      geom_point(size = 2) + 
      geom_text(aes(x = input$niter+0.1, label = impnr), hjust = 0) + 
      transition_reveal(iternr) + 
      coord_cartesian(clip = 'off') + 
      labs(title = "Multiple Imputation Convergence plot", 
           y = "Mean Imputed", x="Iteration number") + 
      theme(legend.position="none") #theme_minimal() #+ 
    #theme(plot.margin = margin(5.5, 40, 5.5, 5.5))
    
  anim_save("outfile.gif", animate(p)) # New
  
  # Return a list containing the filename
  list(src = "outfile.gif",
       contentType = 'image/gif',
        width = 500,
        height = 400
       # alt = "This is alternate text"
  )}, 
  
  deleteFile = TRUE
  
  )}
    #animate(g1, nframes = 100, fps=3)
  

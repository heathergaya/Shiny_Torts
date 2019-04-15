

list.of.packages <- c('rjags', 'runjags', 'mixAK', 'miscF', 'Rdistance', "lme4", "Rcpp", "shinyjs")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages,function(x){library(x,character.only=TRUE)})

js <- "
$(document).ready(function() {
$('#runmodel').on('click', function(){
var date = new Date().toLocaleString();
$('#busy').html('Model started: ' + date);
});
$('#model.ran').on('shiny:value', function(event) {
$('#busy').html('');
});
});
$(document).ready(function() {
$('#extendmodel').on('click', function(){
var date = new Date().toLocaleString();
$('#busy2').html('Model started: ' + date);
});
$('#model.extended').on('shiny:value', function(event) {
$('#busy2').html('');
});
});
"



library(shiny)

# Define server logic ----
server <- function(input, output) {
  output$min_max <- renderText({ 
    paste("You have chosen a range that goes from",
          input$slider[1], "to", input$slider[2])
  })
  output$slidermax <- renderText({ 
    paste("The model will run for a maximum of",
          input$timeslider, "hours")
  })
  output$file <- renderTable({
    head.Burrows()
  })

    
  output$file2 <- renderTable({
    head.Trans()
  })
  
  output$file3 <- renderTable({
    head.Veg()
  })
  
  observeEvent(input$VegButton, {
    if(input$VegButton == "No") {
      toggle(id = "file3")
    }})
  
  observeEvent(input$DetectionButton, {
    if(input$DetectionButton == "Perfect Detection!") {
      toggle(id = "slider")
      toggle(id = "min_max")
    }})
  
  observeEvent(input$file3, {
      inFile3 <- isolate({input$file3})
      Veg <<- read.csv(inFile3$datapath)
               })
  
  observeEvent(input$prepmodel, {
    inFile <- isolate({input$file })
    Burrows <<- read.csv(inFile$datapath)
    inFile2 <- isolate({input$file2 })
    Trans <<- read.csv(inFile2$datapath)
    Freq <<- input$ELTDSButton
    
    one <<- switch(input$DetectionButton, "Perfect Detection!" = 1, "Likely imperfect" = 2)
    two <<- switch(input$ELTDSButton, "Conventional LTDS" = 1, "Enhanced LTDS" = 2)
    three <<- switch(input$VegButton, "Yes" = 1, "No" = 2)
    
    chosenmodel <<- paste(c(one, two, three), collapse = "")
    ### choose which model is being run 
    selected.model <<- ifelse(chosenmodel == "112", "modelstring.PC",
                              ifelse(chosenmodel == "111", "modelstring.PCV",
                                     ifelse(chosenmodel == "122", "modelstring.PE",
                                            ifelse(chosenmodel == "121", "modelstring.PEV",
                                                   ifelse(chosenmodel == "212", "modelstring.IC",
                                                          ifelse(chosenmodel == "211", "modelstring.ICV",
                                                                 ifelse(chosenmodel == "222", "modelstring.IE",
                                                                        "modelstring.IEV")))))))
    
    p.min <<- input$slider[1]
    p.max <<- input$slider[2]
    maxusrtime <<- paste(input$timeslider, "h")
    source('RJMCMC.R')
    source('Bias_Adjusted_Model_412019.R')
    output$model.prepped <- renderText({
      paste(c("Model is Loaded! You have chosen this model: ", selected.model))
    })
  })
  
  observeEvent(input$runmodel, {
    Run.me(Burrows)
  })
  
  observeEvent(input$runmodel, {
    output$model.ran <- renderTable({
      summary(Foo)[,c(1:5,11)]
    }, rownames = TRUE)})
  
  observeEvent(input$extendmodel, {
    Extend.me(Burrows)
    output$model.extended <- renderTable({
      summary(Foo)[,c(1:5,11)]
    }, rownames = TRUE)})
  
  observeEvent(input$plotmodel, {
      Plot.me.shiny(Burrows)
      output$model.plotted <- renderPlot({
        plot(xs, Kall*Occ, type = "l", main = "Estimated Size Curve for Population",
             col = "black", lwd =2,  cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, xlab = "Burrow Width in cm",
             ylab = "", yaxt = "n", ylim = c(-.0002, max(Kall*Occ)+.008))
        percents <- Covmod[c("juvi1", "juvi2", "juvi3"), 4]
        legend("topright", c(paste(round(percents[1],2)*100, "%", "  ", "< 13 cm"), 
                             paste(round(percents[2],2)*100, "%", "  ", "13-22 cm"), 
                             paste(round(percents[3],2)*100, "%", "  ", "> 22 cm")), 
               text.col = c("darkgreen", "purple", "blue"), cex = 1.35)
        arrows(4, -.0002, 12.9, -.0002, length = 0.075, col = "darkgreen", code = 3, angle = 30, lwd = 2)
        arrows(13, -.0002, 21.9, -.0002, length = 0.075, col = "purple", code = 3, angle = 30, lwd = 2)
        arrows(22, -.0002, 65, -.0002, length = 0.075, col = "blue", code = 3, angle = 30, lwd = 2)
    })
    })
  
  head.Burrows <- reactive({
    infile <- input$file
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    head(read.csv(infile$datapath, header = T), n =2)
    
  })
  
  head.Trans <- reactive({
    infile2 <- input$file2
    if (is.null(infile2)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    head(read.csv(infile2$datapath, header = T), n= 2)
  })
  
  head.Veg <- reactive({
    infile3 <- input$file3
    if (is.null(infile3)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    head(read.csv(infile3$datapath, header = T), n= 2)
  })
  
  output$title <- renderText({"Bayesian LTDS Made 'Simple'"})
  
  }



# Define UI ----
ui <- fixedPage(
  useShinyjs(),
  tags$head(
    tags$script(HTML(js)),
    tags$style(type = 'text/css', '#title{font-size: 48px; font-family: calibri light; 
               background-color: rgba(255,255,255,0.40); color: blue; border-style: none;}')
  ),
  
  fixedRow(
    column(width =10, textOutput("title"), offset = 2)),
  
  fixedRow(
    column(radioButtons("DetectionButton", 
              "Is detection on the transect perfect? i.e. is g(0) = 1?", 
              c("Likely imperfect", "Perfect Detection!")), width = 4),
    
  column(width = 4, radioButtons("ELTDSButton", 
                            "Was data collected with standard LTDS (conventional survey) or Enhanced LTDS (double sampling)?", 
                            c("Conventional LTDS", "Enhanced LTDS"))),
  column(width =4, radioButtons("VegButton", 
                            "Does the data include vegetation measurements (or a similar individual covariate)?", 
                            c("Yes", "No")))),
  
  fixedRow(
    column(width = 4, fileInput("file", h3("Burrow Data (.csv)"), 
                            accept=c('text/csv', 'text/comma-separated-values,text/plain'))),
  column(tableOutput("file"), width = 8)),
  
  fixedRow(
  column(fileInput("file2", h3("Transect Data (.csv)"),
                         accept=c('text/csv', 'text/comma-separated-values,text/plain')), width = 4),
  column(tableOutput("file2"), width = 8)),
  
  fixedRow(
  column(fileInput("file3", h3("Veg Data (.csv)"),
                         accept=c('text/csv', 'text/comma-separated-values,text/plain')), width = 4),
  column(tableOutput("file3"), width = 8)),
  
  fixedRow(
    column(sliderInput("slider", label = "Estimated Detection on the Line for 5 cm Burrows",
                            min = 0, max = 1, value = c(.4, .6)), width = 4),
    column(textOutput("min_max"), width = 2),
   
    column(sliderInput("timeslider", label = "Maximum Time for Model Run",
                            min = 0, max = 10, value = 4, step = .25), width =4),
    column(textOutput("slidermax"), width = 2)),
  fixedRow(
    
    column(5, 
            actionButton("prepmodel", "Load the Model")), 
     textOutput("model.prepped")
   ),
   
  fixedRow(
       column(5, 
            actionButton("runmodel", "Run the Model")),
       tags$p(id = "busy")
   ),
   
   fixedRow(
     column(5,
            actionButton("extendmodel", "Extend the Model")),
    tags$p(id = "busy2")
    ),


  fixedRow(
    column(5, 
           actionButton("plotmodel", "Plot the Size Curve"),
           column(7, plotOutput("model.plotted"))
    )),

  fixedRow(
    column(5, align = "center", tableOutput("model.ran"))),

  fixedRow(
    column(5, align = "center", tableOutput("model.extended")))
)


# Run the app ----
shinyApp(ui = ui, server = server)

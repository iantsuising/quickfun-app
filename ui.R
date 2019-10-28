#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("quicfun: QUantIles-based Classifier for FUNctional data"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      
      # upload trainning set,
      fileInput("file1", "Trainning Set: Choose CSV File",
                # an n by (T+1) matrix without column and row names,
                # where n is the sample size, the first T columns are functional predictors observed at T time grids,
                # the last column is the lable.",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      
      fileInput("file2", "Testing Set: Choose CSV File",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      
      helpText("Both CSV files should be an n by (T+1) matrix without column and row names,
              where n is the sample size, the first T columns are functional predictors observed at T time grids,
               the last column is the lable. If data is not uploaded, the results of DTI data from refund package are demonstrated."),
      
      # Horizontal line ----
      tags$hr(),
      submitButton("run", icon("Refresh")),
      sliderInput(inputId = "FVE",
                  label = "FVE (%):",
                  min = 80,
                  max = 99.99,
                  value = 95),
      
      # numericInput("FVE", "FVE (%):", 99, min = 70, max = 99.99),
    
    sliderInput(inputId = "M0",
                label = "Number of Quantiles:",
                min = 3,
                max = 12,
                value = 10),
    

      
      # numericInput("M0",
      #              "Number of Quantiles:",
      #              min = 3,
      #              max = 12,
      #              value = 10),
      
      checkboxInput(
        "pre_smooth", 
        "Pre-smooth the functional predictors.",
        value = FALSE
      )
      
      
    ),
    
    
    
    
    
    
    
    # Show a plot of the generated distribution
    mainPanel(
      
      textOutput("text"),
      plotOutput("distPlot"),
      plotOutput("distCov"),
      plotOutput("distPdf"),
      plotOutput("distMeanVar"),
      DT::dataTableOutput(outputId = "pretable")
    )
  )
))

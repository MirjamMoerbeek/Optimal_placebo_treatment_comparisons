library(shiny)
library(shinydashboard)
library(shinythemes)

ui <- dashboardPage(skin="yellow",   
  dashboardHeader(title = "Optimal placebo-treatment comparisons in trials with an incomplete within-subject design and heterogeneous costs and variances.",disable = TRUE),
  dashboardSidebar(disable = TRUE),
  dashboardBody(
    titlePanel("Optimal placebo-treatment comparisons in trials with an incomplete within-subject design and heterogeneous costs and variances"),
       fluidRow(
         
      column(width=4,  height = 450,
        box(title = "Budgetary constraint", width = NULL, 
          numericInput("C.s", "Overhead costs per subject (C.s)", 500),
          numericInput("C.0", "Costs for control per subject (C.0)", 0),
          numericInput("C.1", "Costs for treatment 1 per subject (C.1)", 423.6),
          numericInput("C.2", "Costs for treatment 2 per subject (C.2)", 635.4),
          numericInput("B", "Budget (B)", 100000),
        )
      ),
      
      column(width=4,height = 450,
             box(title = "Correlation parameters", width = NULL, 
               box(width=6,title="Variances",
                   numericInput("var.0", "Variance of placebo (var.0)", 100),
                   numericInput("var.1", "Variance of treatment 1 (var.1)", 125),
                   numericInput("var.2", "Variance of treatment 2 (var.2)", 150)
               ),
             
               box(width=6,title="Covariances",
                   numericInput("covar.01", "Covariance of placebo and treatment 1 (covar.01)", 33.5),
                   numericInput("covar.02", "Covariance of placebo and treatment 2 (covar.02)", 36.7),
                   numericInput("covar.12", "Covariance of treatment 1 and treatment 2 (covar.12)", 41.1))
               ),
             box(title = "Calculate results",  width=NULL, background = "yellow",
                 submitButton("Submit")
             )
      ),
      
     column(width=4,height = 450,
            box(title = "Optimal design criterion", width = NULL, 
                selectInput("optcrit", "Select optimality criterion", choices = c("Compound optimal design" = 1, "D_A optimal design" = 2)),

             
                  numericInput("lambda.01", "Specify parameter lambda1 of compound optimal design", 0.75, min = 0, max = 1, step = 0.05)
             )
     )
    ),

    fluidRow(  

   column(width=12,height = 450,
         box( title = "Results for the incomplete within-subject design with scenario 1",width=4,height=540,
              tableOutput("table1"), ),
         box( title = "Results for the incomplete within-subject design with scenario 2",width=4,height=540,
              tableOutput("table2"), ),
              )
              )
    ))
  


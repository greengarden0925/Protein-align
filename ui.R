version = "1.0.0"

dashboardPage( 
  skin="yellow",
  dashboardHeader(
    title=paste0("Protein sequenceis alignment (v ",version,")"),
    titleWidth = 450
  ),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Main output", tabName = "PLOT", icon = icon("dashboard")),
      fluidRow(
        column(12, align="center",
               h5(p(strong(span(style="color:white","Target sequenceis file (.fasta):")))),
               fileInput(inputId="file1", "", multiple=FALSE),
               h5(p(strong(span(style="color:white","Reference sequenceis file (.fasta):")))),
               fileInput(inputId="file2", "", multiple=FALSE),
               actionButton("submit",strong("Start to analyze!"),icon("list-alt")),
               tags$style(type='text/css', "#submit { vertical-align: middle; height: 40px; width: 85%; font-size: 20px;}")
        )
      )
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "PLOT",
              tags$script('
                              $(document).on("shiny:connected",function(e){
                          var Width = window.screen.width;
                          Shiny.onInputChange("MonitorWidth",Width);})
                          '),
              tags$script('
                          $(document).on("shiny:connected",function(e){
                          var Height = window.screen.height;
                          Shiny.onInputChange("MonitorHeight",Height);})
                          '),
              fluidRow(
                column(4, align="center",
                       h4(p(strong(span(style="color:green","Please select a cut-point:")))),
                       uiOutput("choose_columns1")
                ),
                column(4, align="center",
                       downloadButton("download1", label = "Download text file", class = NULL),
                       tags$style(type='text/css', "#download1 { vertical-align: middle; height: 40px; width: 85%; font-size: 20px;}")
                ),
                column(4, align="center",
                       downloadButton("download2", label = "Download plot (pdf)", class = NULL),
                       tags$style(type='text/css', "#download2 { vertical-align: middle; height: 40px; width: 85%; font-size: 20px;}")
                )
              ),
              fluidRow(
                column(6, align="center",
                       htmlOutput("plot1")
                ),
                column(6, align="center",
                       htmlOutput("plot2")
                )
              )
      )
    )
  )
)

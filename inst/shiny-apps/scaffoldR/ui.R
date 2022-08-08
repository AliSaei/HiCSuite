
source("./globals/functions.R")  
source("./globals/helpers.R")

shinyUI(
  dashboardPage(
    skin = "black",
    dashboardHeader(
      disable = TRUE,
      title = "Hi-C Data Analysis"
      #class = "dropdown"),
      #tags$li(a(href = 'https://github.com/AliSaei', target="_blank", icon("github"), title = "Code on Github"),
      #class = "dropdown")
    ),
    dashboardSidebar(sidebarMenu(
      div(id = "Logo", 
          style = "margin-bottom: 5px; padding:1px; border:1px solid gray;",
          img(src='logo.png', width = '100%')
      ),
      menuItem(span("Input", style = "margin-left:0px;"), tabName = "Data", icon = icon("table")),
      menuItem("Scaffolder", tabName = "SCFLD", icon = icon("th", lib = "glyphicon"))
      #tags$li(a(href="mailto:saeibox@gmail.com", target="_top", icon("envelope"), title = "Any questions? Please email!"),
      #class = "dropdown"),
    )
    ),
    dashboardBody(
      shinyjs::useShinyjs(),
      #tags$head(tags$link(rel="shortcut icon", href="logo.png")),
      fluidRow(id = "tabItems", style = "padding-top: 0px; font-family: Tahoma; margin-top: -10px;",
               tabItems(
                 tabItem("Data",
                         source(file.path("./globals", "input_ui.R"),  local = TRUE)$value
                 ),
                 tabItem("SCFLD",
                         source(file.path("./globals", "scaf_ui.R"),  local = TRUE)$value
                 )
               )
      ),

      ## javascript
      tags$script(src="jscript.js"),
      
      ## helper function css
      tags$head(tags$style(HTML("
                            .btn-loading-container {
                            text-align: center;
                            margin-top: -20px;
                            padding: 0px;
                            font-size: 1.2em;
                            }
                            .btn-err {
                            margin-top: 10px;
                            color: red;
                            }
                                "))),
      
      # Also add some custom CSS to make the title background area the same
      # color as the rest of the header.
      tags$head(tags$style(HTML("
                           .btn-action{font-size: 12px; font-family: Tahoma; width: 100%; font-weight: 550; background-color: green; color: white;}
                           .btn-update{font-size: 12px; font-family:Tahoma; width: 100%; font-weight: 550; background-color: #ece9df;}
                           .btn-dir{font-size: 12px; font-family:Tahoma; width: 100%; background-color: #F8F8F8;}
                           .btn-picker{font-size: 12px; font-family: Tahoma; background-color: #F8F8F8; z-index: 100;}
                           .btn-accordion{font-size: 12px; font-family:Tahoma; font-weight: 550; width: 100%; height: 30px; background-color: #E8E8E8;}
                           .btn-drop{font-size: 12px; font-family:Tahoma; font-weight: 550; width: 100%; height: 30px; background-color: #E8E8E8; text-align: left;}
                           .selectize-input { font-size: 12px; font-family: Tahoma;}
                           .selectize-dropdown { font-size: 12px; font-family: Tahoma; overflow-x: auto;}
                           .selectize-dropdown-content {max-height: 300px; padding: 0; }
                           .selectize-dropdown-content .option {border-bottom: 1px solid #ccc;}
                           .selectize-dropdown-content .active {background: #2196f3 !important;}
                           .fa {font-size: 14px;}
                           .glyph-icon {font-size: 12px;}
                           .slide-btn{font-size: 12px; background-color: #ece9df; 
                              border-top-right-radius: 30px; border: 1px solid #E8E8E8; 
                              padding-top: 5px; padding-bottom: 5px; font-weight: 550; 
                              margin-top: 5px; margin-left: -5px;}

                            table.dataTable tr.selected td, table.dataTable td.selected {background-color: #C8C8C8 !important;}
                           .skin-black .main-sidebar {background-color: white; position: fixed; width: 115px;}
                           .skin-black .main-sidebar .sidebar .sidebar-menu a {
                              background-color: white;
                              color: black;
                              font-family: Tahoma;
                              font-size: 12px;
                              font-weight: bold;}
                            .skin-black .main-sidebar .sidebar .sidebar-menu .active a {
                                                    background-color: #ece9df;
                                                    border-left-color: #ff0000;
                      
                            }
                            .skin-black  .nav-tabs-custom .nav-tabs .active li a {
                                                  background-color: white;
                            }
                            .skin-black  .nav-tabs-custom .nav-tabs .active:hover a {
                                                    background-color: white;
                            }
                            /* body */
                            .content-wrapper, .right-side {
                                                    background-color: white;
                                                    margin-left: 100px;
                                                    height: 100vh;
                                                    overflow: auto;
                            }
                            .dbutt{width: 100%; background: #E8E8E8 !important; font-size: 12px; margin-top: 2px;}
                            #map2_tooltip {position: absolute; pointer-events: none; width: 300px; z-index: 100; padding: 0;}
                            #table.dataTable {font-size: 12px; position:relative; !important;}
                            #.progress-bar {color: transparent !important}

      "))),
      
      tags$script('
      var dimension = [0, 0];
      $(document).on("shiny:connected", function(e) {
        dimension[0] = window.innerWidth;
        dimension[1] = window.innerHeight;
        Shiny.onInputChange("dimension", dimension);
      });
      $(window).resize(function(e) {
        dimension[0] = window.innerWidth;
        dimension[1] = window.innerHeight;
        Shiny.onInputChange("dimension", dimension);
      });
      ')    )
  ))

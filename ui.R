source("./globals/functions.R")
source("./globals/helpers.R")
shinyUI(fluidPage(
  div(style = "position: fixed; height: 100%; background: #ece9df; padding: 10px; font-size: 12px; width: 250px; border-right: 1px solid lightgray; font-family: Tahoma;",
      fluidRow(style = "margin: 0;",
               tags$button(id = "dtInput", class = "btn btn-default btn-md btn-accordion",
                           list(NULL, label = "Data Input"),
                           htmltools::tags$span(class = "glyphicon glyphicon-triangle-bottom")
               ),
               div(style = "border: 1px solid lightgray; border-radius: 3px; padding: 2px;  background-color: white; margin-bottom: 10px;",
                   div(id =  "DataInput", style = "padding: 10px;",
                       radioGroupButtons("mapSrc", "Upload contact data:", choices = c("Local", "Remote"), selected = "Local", justified = TRUE, size = "xs",
                                         checkIcon = list(
                                           yes = tags$i(class = "fa fa-circle", 
                                                        style = "color: steelblue"),
                                           no = tags$i(class = "fa fa-circle-o", 
                                                       style = "color: steelblue"))),
                       div(style = "margin-top: -13px;",
                           conditionalPanel(
                             condition = "input.mapSrc == 'Local'",
                             #shiny::fileInput("contactMap", NULL, accept = c("csv", ".txt", ".rds"), width = '100%'),
                             shinyDirButton('directory', label='Set Working Directory', title='Please select a folder', class = "btn-dir", icon = icon("folder-open")),
                             #verbatimTextOutput("dirPath"),
                             div(style = "margin-top: 0px;",
                                 pickerInput("mapFile", NULL, choices = NULL , multiple = FALSE,
                                             options = list(style = "btn-default btn-md btn-picker", size = 20, `live-search` = TRUE,  tickIcon = "glyphicon-ok", width = "100%",showTick = TRUE)
                                 )
                             )
                           ),
                           conditionalPanel(
                             condition = "input.mapSrc == 'Remote'"
               
                           )
                       ),
                       div(
                       shinyWidgets::awesomeCheckbox("loadSeqLen", strong("Upload sequence lengths"))
                       ),
                       
                       conditionalPanel(
                         condition = "input.loadSeqLen != 0",
                         div(style = "margin-top: -15px;",
                         radioGroupButtons("lenSrc", NULL, choices = c("Local", "Remote"), selected = "Local", justified = TRUE, size = "xs",
                                           checkIcon = list(
                                             yes = tags$i(class = "fa fa-circle", 
                                                          style = "color: steelblue"),
                                             no = tags$i(class = "fa fa-circle-o", 
                                                         style = "color: steelblue"))),
                         div(style = "margin-top: -13px;",
                             conditionalPanel(
                               condition = "input.lenSrc == 'Local'",
                               #shiny::fileInput("lenFile", NULL, accept = c("csv", ".txt"), width = '100%'),
                               div(style = "margin-top: 0px;",
                                   pickerInput("lenFile", NULL, choices = NULL , multiple = FALSE,
                                               options = list(style = "btn-default btn-md btn-picker", size = 20, `live-search` = TRUE,  tickIcon = "glyphicon-ok", width = "100%",showTick = TRUE)
                                   )
                               )
                             ),
                             conditionalPanel(
                               condition = "input.lenSrc == 'Archive'"
                            
                             )
                         )
                       )
                       ),
                       
                       div(style = "margin-top: 0px;",
                           withBusyIndicatorUI(
                             actionButton("import_data", "Import Data", icon = icon('import', lib = "glyphicon"), class = "action_button", width = '100%')
                           )
                       )
                   )
               ),
               tags$button(id = "dtBin", class = "btn btn-default btn-md btn-accordion",
                           list(NULL, label = "Data Binning"),
                           htmltools::tags$span(class = "glyphicon glyphicon-triangle-bottom")
               ),
               div(style = "border: 1px solid lightgray; border-radius: 3px; padding: 2px;  background-color: white; margin-bottom: 10px;",
                   div(id = "DtBin" , style = "padding: 10px; display: none;",
                       div(style = "height:135px;",
                           fluidRow(style = "margin: 0px;",
                                    numericInput("binsize2", "Bin Size:", value = 1, min = 0, max = 500000, step= 10000),
                                    div(style = "margin-top: -13px;",
                                        sliderInput("binsize", NULL, min = 0, value = 0, max = 0, width = '100%', ticks = TRUE)
                                    ),
                                    div(style = "margin-top: -30px;",
                                        withBusyIndicatorUI(
                                          actionButton("binUpdate", "Update Binning", icon = icon('play', lib = "glyphicon"), width = '100%', class = "btn-update")
                                        )
                                    )
                           )
                       )
                   )
               )
      )
      #tags$hr()
      #p("Batch run: plot Hi-C maps of 20 most likley sequences to follow the selected sequence"),
      #actionButton("build_map", "Run", width = '100%', class = "action_button")
  ),
  div(style = "margin-left: 250px;",   
      div(
        actionButton("vwIntraMap", "Intramap View",icon = icon("table"), width = "260px"),
        tags$style(HTML("#vwIntraMap{font-size: 12px; background-color:#ece9df; border-top-right-radius: 30px; padding-top: 5px; padding-bottom: 5px; font-weight: 550; margin-top: 2px;}"))
      ),
      fluidRow(style = "margin: 0 0 0 10px; border: 1px solid #F0F0F0; border-radius: 5px; padding: 1px; font-family: Tahoma; font-size: 12px; ",
               div(id = "IntraMap", style = "height: calc(100% + 160px);",
                   div(style = "border: 1px solid lightgray; width: 100%; padding: 5px; background-color: #ece9df; border-radius: 5px; font-size: 12px; width: 250px; float: left;",
                       #tags$button(class = "btn btn-default btn-md btn-accordion",
                       #            list(NULL, label = ""),
                       #            htmltools::tags$span(class = "glyphicon glyphicon-triangle-bottom")
                       #),
                       div(style = "border: 1px solid lightgray; border-radius: 3px; padding: 2px;  background-color: white; margin-bottom: 5px;",
                           div(id =  "Scaffold", style = "padding: 10px;",
                               radioGroupButtons("dir1", "Sequence:", size = "xs",
                                                 choices = c(`<i class='fa fa-arrow-left'></i>` = "Backward", `<i class='fa fa-arrow-right'></i>` = "Forward"), 
                                                 selected = "Forward", justified = TRUE
                               ),
                               div(style = "margin-top: -15px;", 
                                   pickerInput("seq1", NULL, choices = "", width = '100%', multiple = FALSE,
                                               options = list(style = "btn-default btn-md btn-picker", size = 10, `live-search` = TRUE))
                               ),
                               conditionalPanel(
                                 condition = "input.seq1.length > 0",
                                 div(style = "margin-top: -10px;",
                                     radioGroupButtons("action", "Action:", choices = c("Cut", "Join"), selected = "Cut", justified = TRUE, size = "xs",
                                                       checkIcon = list(
                                                         yes = tags$i(class = "fa fa-circle", 
                                                                      style = "color: steelblue"),
                                                         no = tags$i(class = "fa fa-circle-o", 
                                                                     style = "color: steelblue"))
                                     ),
                                     conditionalPanel(
                                       condition = "input.action == 'Cut'",
                                       div(style = "border: 1px solid lightgray; border-radius: 3px; padding: 5px; margin-top: -15px; margin-bottom: 10px;",
                                           numericInput("cutPos", "Position:", value = 0, min = 0, max = 50000000, step= 100),
                                           withBusyIndicatorUI(
                                             actionButton("cut1", "Cut Sequence", icon = icon('scissors', lib = "glyphicon"), width = '100%', class = "action_button")
                                           )
                                       )
                                     ),
                                     conditionalPanel(
                                       condition = "input.action == 'Join'",
                                       div(style = "margin-top: -10px;",
                                           tags$button(id = "intConfig1", class = "btn btn-default btn-md btn-accordion",
                                                       list(NULL, label = "Inter-sequence Interaction"),
                                                       htmltools::tags$span(class = "glyphicon glyphicon-triangle-bottom")
                                           ),
                                           div(style = "border: 1px solid lightgray; border-radius: 3px; padding: 2px;  background-color: #ece9df; margin-bottom: 10px;",
                                               div(id ="IntConfig1", style = "padding: 5px;",
                                                   fluidRow(style = "margin: 0px;",
                                                            numericInput("edgeSize1", "Edge Size:", value = 1, min = 0, max = 500000, step= 10000),
                                                            div(style = "margin-top: -13px;",
                                                                sliderInput("edgeSize2", NULL, min = 0, value = 0, max = 0, width = '100%', ticks = TRUE)
                                                            ),
                                                            div(style = "margin-top: -30px;",
                                                                withBusyIndicatorUI(
                                                                  actionButton("calcIntraction1", "Calculate Interactions", icon = icon('refresh', lib = "glyphicon"), width = '100%', class = "btn-update")
                                                                )
                                                            )
                                                   )
                                               )
                                           ),
                                           conditionalPanel(
                                             condition = "input.calcIntraction1 != 0 | input.calcIntraction2 != 0",
                                             div(style = "border: 1px solid lightgray; border-radius: 3px; padding: 5px; margin-top: -5px;",
                                                 pickerInput("subseq1", "Subsequent sequence:", choices = NULL, width = '100%', multiple = FALSE,
                                                             options = list(style = "btn-default btn-md btn-picker", size = 10, `live-search` = TRUE)),
                                                 withBusyIndicatorUI(
                                                   actionButton("join", "Join Sequences", icon = icon('paired', lib = "glyphicon"), width = '100%', class = "action_button")
                                                 )
                                             )
                                           )
                                       )
                                     )
                                 )
                               )
                           )
                       ),
                       conditionalPanel(
                         condition = "input.seq1.length > 0",
                         downloadButton("dlChanges",label = "Download Changes" , class = "dbutt"), 
                         actionButton("svChanges", "Save Changes", width = '100%', class = "dbutt", icon = icon("floppy-disk", lib = "glyphicon")), 
                         tags$head(tags$style(".dbutt{width: 100%; background: lightgray !important; font-size: 12px; margin-top: 2px;}"))
                       )
                   ),
                   div(style = "width: calc(100% - 250px); float: left;",
                       div(style = "width: 82vh; border: 1px solid #F0F0F0; border-radius: 5px; padding: 10px; background: #ece9df; margin: 0px 5px 5px 5px; font-size: 12px; font-family: Tahoma;",
                           actionBttn("exitZoom", NULL, icon = icon('home', lib = "glyphicon"), style = "material-circle", size = "xs"),
                           actionBttn("clearBrush", NULL, icon = icon('erase', lib = "glyphicon"), style = "material-circle", size = "xs"),
                           actionBttn("cut2", NULL, icon = icon('scissors', lib = "glyphicon"), style = "material-circle", size = "xs"),
                           div(style = "border: 1px solid white; border-radius: 5px; margin: 0px; width: 80vh; margin-top: 5px;",
                               #sliderInput("range", label = NULL, min = 0, max = 100, value = c(0, 100), ticks = FALSE, width = '100%'),
                               plotOutput("intramap", width = "auto", height = "auto", dblclick = "intramap_dblclick",  brush = "intramap_brush", click = "intramap_click")
                           )
                       )
                   )
               )
      ),
      div(
        actionButton("vwScaf", "Scaffold View",icon = icon("table"), width = "260px"),
        tags$style(HTML("#vwScaf{font-size: 12px; background-color:#ece9df; border-top-right-radius: 30px; padding-top: 5px; padding-bottom: 5px; font-weight: 550; margin-top: 2px;}"))
      ),
      fluidRow(style = "margin: 0 0 0 10px; border: 1px solid #F0F0F0; border-radius: 5px; padding: 1px; font-family: Tahoma; font-size: 12px; ",
               div(id = "VwScaf", style ="height: calc(100% + 160px); display: none;",
                   div(style = "border: 1px solid lightgray; width: 100%; padding: 7px; background-color: #ece9df; border-radius: 5px; font-size: 12px; width: 250px; float: left;",
                       #div("Setting", icon("triangle-bottom", class = "glyph-icon", lib = "glyphicon")),
                       tags$button(id = "intConfig2", class = "btn btn-default btn-md btn-accordion",
                                   list(NULL, label = "Inter-sequence Interaction"),
                                   htmltools::tags$span(class = "glyphicon glyphicon-triangle-bottom")
                       ),
                       div(style = "border: 1px solid lightgray; border-radius: 3px; padding: 2px;  background-color: white; margin-bottom: 10px;",
                           div(id = "IntConfig2", style = "padding: 5px;",
                               fluidRow(style = "margin: 0px;",
                                        numericInput("edgeSize3", "Edge Size:", value = 1, min = 0, max = 500000, step= 10000),
                                        div(style = "margin-top: -13px;",
                                            sliderInput("edgeSize4", NULL, min = 0, value = 0, max = 0, width = '100%', ticks = TRUE)
                                        ),
                                        div(style = "margin-top: -30px;",
                                            withBusyIndicatorUI(
                                              actionButton("calcIntraction2", "Calculate Interactions", icon = icon('play', lib = "glyphicon"), width = '100%', class = "action_button")
                                            )
                                        )
                               )
                           )
                       ),
                       conditionalPanel(
                         condition = "input.calcIntraction2 || input.calcIntraction1",
                         tags$button(id = "scafConfig", class = "btn btn-default btn-md btn-accordion",
                                     list(NULL, label = "Scaffolding"),
                                     htmltools::tags$span(class = "glyphicon glyphicon-triangle-bottom")
                         ),
                         div(style = "border: 1px solid lightgray; border-radius: 3px; padding: 2px;  background-color: white; ",
                             div(id =  "ScafConfig", style = "padding: 5px; display: none;",
                                 div(style = "height:310px;",
                                     fluidRow(style = "margin: 0px;",
                                              radioGroupButtons("dir2", "Leading sequence:", size = "xs",
                                                                choices = c(`<i class='fa fa-arrow-left'></i>` = "Backward", `<i class='fa fa-arrow-right'></i>` = "Forward"), 
                                                                selected = "Forward", justified = TRUE,
                                                                checkIcon = list(
                                                                  yes = tags$i(class = "fa fa-circle", 
                                                                               style = "color: steelblue"),
                                                                  no = tags$i(class = "fa fa-circle-o", 
                                                                              style = "color: steelblue"))
                                              )
                                     ),
                                     fluidRow(style = "margin: -15px 0 0 0;", 
                                              pickerInput("seq2", NULL, choices = "", width = '100%', multiple = FALSE,
                                                          options = list(style = "btn-default btn-md btn-picker", size = 10, `live-search` = TRUE)),
                                              #selectInput("seq_1", "Leading sequence:", choices = NULL, width = '100%'),
                                              div(style = "margin-top: -15px;",
                                                  shinyWidgets::switchInput("orient_1", label = "Strand", offLabel = "-", onLabel = "+", value = TRUE, size = "mini", onStatus = "primary", offStatus = "primary", labelWidth = '100', handleWidth = '100')
                                                  #checkboxInput("reverse_seq.1", "Reverse initial sequence")
                                              )
                                              
                                     ),
                                     fluidRow(style = "margin: 0px;", 
                                              p("Subsequent sequence:", style = "font-weight: 550;"),
                                              div(style = "margin: 0;",
                                                  sliderInput("maxLinkDen", NULL, min = 0, value = 1, max = 1, width = '100%', ticks = FALSE, pre = "Link density:")
                                              ),
                                              div(style = "margin-top: -14px;",
                                                  div(style = "float: left; width: calc(100% - 70px);",
                                                      pickerInput("subseq2", NULL, choices = NULL, width = '100%',
                                                                  options = list(style = "btn-default btn-md btn-picker", size = 10, `live-search` = TRUE)
                                                      ),
                                                      div(style = "margin-top: -15px;",
                                                          shinyWidgets::switchInput("orient_2", label = "Strand", offLabel = "-", onLabel = "+", value = TRUE, size = "mini", onStatus = "primary", offStatus = "primary", labelWidth = '70', handleWidth = '70')
                                                      )
                                                  ),
                                                  div(style = "float: right; width: 70px; margin: 0px; border: 1px solid lightgray; padding: 1px 1px 1px 5px; border-radius: 3px; background-color: #ece9df; height: 33px;",
                                                      actionBttn("down", NULL, icon = icon('arrow-down', lib = "glyphicon"), style = "material-circle",  size = "xs"),
                                                      actionBttn("up", NULL, icon = icon('arrow-up', lib = "glyphicon"),style = "material-circle", size = "xs")
                                                      #checkboxInput("reverse_seq.1", "Reverse initial sequence")
                                                  )
                                              )
                                     ),
                                     #textInput("chr",label = "Assembled sequence:", value = "", width = '100%'),
                                     #tags$hr(), 
                                     div(style = "margin-top: 0px;",
                                         withBusyIndicatorUI(
                                           actionButton("add", "Add to scaffold", width = '100%', class = "action_button", icon = icon("paired", lib = "glyphicon"))
                                         )
                                     )
                                 ),
                                 div(style = "margin-top: 20px; margin-bottom: 10px;",
                                     fluidRow(style = "margin: 0px;",
                                              div(style = "width: 35px; float: left;",
                                                  actionBttn("clipbtn", NULL, icon = icon('copy', lib = "glyphicon"), style = "material-circle", size = "xs"),
                                                  actionBttn("edit", NULL, icon = icon('edit', lib = "glyphicon"), style = "material-circle", size = "xs"),
                                                  shinyjs::hidden(
                                                    actionBttn("check", NULL, icon = icon('check', lib = "glyphicon"), style = "material-circle", size = "xs")
                                                  ),
                                                  actionBttn("erase", NULL, icon = icon('erase', lib = "glyphicon"), style = "material-circle", size = "xs")
                                              ),
                                              div(style = "width: calc(100% - 35px); float: left;",
                                                  div(id = "CheckBox", style = "border: 1px solid #F0F0F0; padding: 5px; border-radius: 3px; min-height: 85px; max-height: 300px; overflow-y: auto; background: white;",
                                                      checkboxGroupInput("scaf_auto", "Scaffold:", choices = NULL)
                                                  ),
                                                  hidden(
                                                    div(id = "EditBox" , style = "margin-top: 0px;",
                                                        textAreaInput("scaf_edit", label = NULL, value = "", placeholder = "please enter one id per line", width = '100%', height = "300px")
                                                    )
                                                  )
                                              ),
                                              conditionalPanel(
                                                condition = "input.scaf_auto.length > 1",
                                                div("Leading sequence number:", style = "width: 160px; float: left; padding: 6px; border: 1px solid lightgray; height: 34px; border-radius: 3px; background: #F4F4F4; margin: top: 2px;"),
                                                div(style = "width: calc(100% - 160px); float: left;",
                                                    numericInput("nrSeq", NULL, value = 1, min = 1, max = 10, step= 1)
                                                )
                                              ),
                                              conditionalPanel(
                                                condition = "input.scaf_auto.length > 1",
                                                div(style ="margin-top: 10px;",
                                                    withBusyIndicatorUI(
                                                      actionButton("combineMaps", "Build scaffold map", width = '100%', class = "action_button")
                                                    )
                                                )
                                              )
                                     ),
                                     conditionalPanel(
                                       condition = "input.inputType == 'Manual'",
                                       div(style = "margin-top: -15px;",
                                           textAreaInput("scaf_man", label = NULL, value = "", placeholder = "please enter one id per line", width = '100%', height = "300px")
                                       )
                                     )
                                 )
                             )
                         )
                       )
                   ),
                   div(style = "margin-left: 5px; float: left; width: calc(100% - 305px);",  
                       div(style = "border-left: 2px solid lightgray;",
                           div(
                             actionButton("vwJointMap", "View Map",icon = icon("table"), width = "200px"),
                             tags$style(HTML("#vwJointMap{font-size: 12px; background-color:#ece9df; border-top-right-radius: 20px; padding-top: 5px; padding-bottom: 5px; font-weight: 550; margin-top: 0px; margin-left: -5px;}"))
                           ),
                           div(style = "width: 82vh; border: 1px solid #F0F0F0; border-radius: 5px; background: #ece9df;font-size: 12px; font-family: Tahoma; margin: 0px 5px 0px 5px; ",
                               div(id = "VwJointMap", style ="padding: 10px; display: none;", 
                                   actionBttn("exitZoom1", NULL, icon = icon('home', lib = "glyphicon"), style = "material-circle", size = "xs"),
                                   actionBttn("clearBrush1", NULL, icon = icon('erase', lib = "glyphicon"), style = "material-circle", size = "xs"),
                                   div(style = "border: 1px solid white; border-radius: 5px; margin: 0px; width: 80vh; margin-top: 5px;",
                                       plotOutput("map1", width = "auto", height = "auto", brush = "map1_brush",hover = hoverOpts(id = "map1_hover")) 
                                   )
                               )
                           ),
                           div(
                             actionButton("vwIntData", "View Data",icon = icon("table"), width = "200px"),
                             tags$style(HTML("#vwIntData{font-size: 12px; background-color:#ece9df; border-top-right-radius: 20px; padding-top: 5px; padding-bottom: 5px; font-weight: 550; margin-top: 0px; margin-left: -5px;}"))
                           ),
                           div(style = "width: 82vh; border: 1px solid #F0F0F0; border-radius: 5px; padding: 0px; background: #ece9df; margin: 0px 5px 0px 5px; font-size: 12px; font-family: Tahoma;",
                               div(id = "VwIntData", style = "border: 1px solid white; border-radius: 5px; margin: 0px; display: none;",
                                   DT::DTOutput("interseqData") 
                               )
                           ),
                           div(
                             actionButton("vwJointMap2", "View Scaffold Map",icon = icon("table"), width = "200px"),
                             tags$style(HTML("#vwJointMap2{font-size: 12px; background-color:#ece9df; border-top-right-radius: 20px; padding-top: 5px; padding-bottom: 5px; font-weight: 550; margin-top: 0px; margin-left: -5px;}"))
                           )
                       ),
                       div(style = "width: 82vh; border: 1px solid #F0F0F0; border-radius: 5px; background: #ece9df; margin: 0px 5px 0px 5px;  font-size: 12px; font-family: Tahoma;",
                           div(id = "VwJointMap2", style ="padding: 10px; display: none;", 
                               div(style = "border: 1px solid white; border-radius: 5px; margin: 0px; width: 80vh;",
                                   plotOutput("map2", width = "auto", height = "auto", hover = hoverOpts(id = "map2_hover")) %>%
                                     withLoader(type="html", loader="loader5")
                               ),
                               uiOutput("map2_hoverinfo")
                           )
                       )
                   )
               )
      )
  ),
  shinyjs::useShinyjs(),
  tags$style(appCSS),
  tags$script(src="jscript.js"),
  tags$head(tags$script('
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
                        ')),
  tags$head(tags$style(HTML('
                           .action_button{font-size: 12px; font-family:Tahoma; width: 100%; font-weight: 550; background-color: green;}
                           .btn-update{font-size: 12px; font-family:Tahoma; width: 100%; font-weight: 550; background-color: #ece9df;}
                           .btn-dir{font-size: 12px; font-family:Tahoma; width: 100%; background-color: #F8F8F8;}
                           .btn-picker{font-size: 12px; font-family:Tahoma; background-color: #F8F8F8;}
                           .btn-accordion{font-size: 12px; font-family:Tahoma; font-weight: 550; width: 100%; height: 30px; background-color: #F0F0F0;}
                           .selectize-input { font-size: 12px; font-family: Tahoma;}
                           .selectize-dropdown { font-size: 12px; font-family: Tahoma; overflow-x: auto;}
                           .selectize-dropdown-content {max-height: 300px; padding: 0;}
                           .selectize-dropdown-content .option {border-bottom: 1px dotted #ccc;}
                           .fa {font-size: 14px;}
                           .glyph-icon {font-size: 12px;}
                            #map2_tooltip {position: absolute; pointer-events:none; width: 300px; z-index: 100; padding: 0;}
                            '))
  )
  
)
)

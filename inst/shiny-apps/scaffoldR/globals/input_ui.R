
shinyUI(
  fluidPage(style = "font-family: Tahoma;",
            div(style = "position: fixed; height: 100%; background: #ece9df; padding: 10px; width: 250px; border-right: 1px solid #E8E8E8; font-family: Tahoma; font-size: 12px;",
                fluidRow(style = "margin: 10px 0 0 0;",
                         tags$button(id = "dtInput", class = "btn btn-default btn-md btn-accordion",
                                     list(NULL, label = "Data Input"),
                                     htmltools::tags$span(class = "glyphicon glyphicon-triangle-bottom")
                         ),
                         div(style = "border: 1px solid #E8E8E8; border-radius: 3px; padding: 2px;  background-color: white; margin-bottom: 10px;",
                             div(id =  "DataInput", style = "padding: 10px;",
                                 radioGroupButtons("mapSrc", "Contact data:", choices = c("Local", "Remote"), selected = "Local", 
                                                   justified = TRUE, size = "xs", status = "gray", 
                                                   checkIcon = list(
                                                     yes = tags$i(class = "fa fa-circle", 
                                                                  style = "color: black"),
                                                     no = tags$i(class = "fa fa-circle-o", 
                                                                 style = "color: white"))),
                                 div(style = "margin-top: -13px;",
                                     conditionalPanel(
                                       condition = "input.mapSrc == 'Local'",
                                       #shiny::fileInput("contactMap", NULL, accept = c("csv", ".txt", ".rds"), width = '100%'),
                                       shinyFiles::shinyDirButton('directory', label='Set Data Directory', 
                                                                  title='Please select a folder', 
                                                                  class = "btn-dir", icon = icon("folder-open")),
                                       #verbatimTextOutput("dirPath"),
                                       div(style = "margin-top: 0px;",
                                           pickerInput("mapFile", NULL, choices = NULL , multiple = FALSE,
                                                       options = list(style = "btn-default btn-md btn-picker", size = 20, 
                                                                      `live-search` = TRUE,  tickIcon = "glyphicon-ok", 
                                                                      width = "100%", showTick = TRUE)
                                           )
                                       )
                                     )
                                 ),
                                 div(
                                   shinyWidgets::prettySwitch("loadSeqLen", strong("Sequence lengths (optional)"), status = "success")
                                 ),
                                 
                                 conditionalPanel(
                                   condition = "input.loadSeqLen != 0",
                                   div(style = "margin-top: -10px;",
                                       radioGroupButtons("lenSrc", NULL, choices = c("Local", "Remote"), 
                                                         selected = "Local", justified = TRUE, size = "xs", status = "gray",
                                                         checkIcon = list(
                                                           yes = tags$i(class = "fa fa-circle", 
                                                                        style = "color: black"),
                                                           no = tags$i(class = "fa fa-circle-o", 
                                                                       style = "color: white"))),
                                       div(style = "margin-top: -13px;",
                                           conditionalPanel(
                                             condition = "input.lenSrc == 'Local'",
                                             #shiny::fileInput("lenFile", NULL, accept = c("csv", ".txt"), width = '100%'),
                                             div(style = "margin-top: 0px;",
                                                 pickerInput("lenFile", NULL, choices = NULL , multiple = FALSE,
                                                             options = list(style = "btn-default btn-md btn-picker", size = 20, 
                                                                            `live-search` = TRUE,  tickIcon = "glyphicon-ok", 
                                                                            width = "100%",showTick = TRUE)
                                                 )
                                             )
                                           )
                                       )
                                   )
                                 ),
                                 div(style = "margin-top: 0px;",
                                     withBusyIndicatorUI(
                                       actionButton("import_data", "Import Data", class = "btn-action",
                                                    icon = icon('import', lib = "glyphicon"))
                                     )
                                 ),
                                 div(style = "margin-top: 15px;",
                                     numericInput("binsize2", "Bin size (bp):", value = 1, min = 0, max = 500000, step= 10000),
                                     div(style = "margin-top: -13px;",
                                         sliderInput("binsize", NULL, min = 0, value = 0, max = 0, width = '100%', ticks = TRUE)
                                     ),
                                     div(style = "margin-top: -30px;",
                                         withBusyIndicatorUI(
                                           actionButton("updateBinning", "Update", class = "btn-update",
                                                        icon = icon('play', lib = "glyphicon"))
                                         )
                                     ),
                                     conditionalPanel(
                                       condition = "input.updateBinning",
                                       withBusyIndicatorUI(
                                         actionButton("saveBinData1", "Save changes to disk", width = '100%', class = "dbutt", 
                                                      icon = icon("floppy-disk", lib = "glyphicon")
                                         )
                                       )
                                     )
                                 )
                             )
                         ),
                         tags$button(id = "faInput", class = "btn btn-default btn-md btn-accordion",
                                     list(NULL, label = "Assembly Input (Optional)"),
                                     htmltools::tags$span(class = "glyphicon glyphicon-triangle-bottom")
                         ),
                         div(style = "border: 1px solid #E8E8E8; border-radius: 3px; padding: 2px;  background-color: white; margin-bottom: 10px;",
                             div(id =  "FAInput", style = "padding: 10px; display: none;",
                                 radioGroupButtons("faSrc", "FASTA file:", choices = c("Local", "Remote"), selected = "Local", 
                                                   justified = TRUE, size = "xs", status = "gray", 
                                                   checkIcon = list(
                                                     yes = tags$i(class = "fa fa-circle", 
                                                                  style = "color: black"),
                                                     no = tags$i(class = "fa fa-circle-o", 
                                                                 style = "color: white"))),
                                 div(style = "margin-top: -13px;",
                                     conditionalPanel(
                                       condition = "input.faSrc == 'Local'",
                                       #shiny::fileInput("contactMap", NULL, accept = c("csv", ".txt", ".rds"), width = '100%'),
                                       #shinyFiles::shinyDirButton('directory', label='Set Data Directory', 
                                       #                           title='Please select a folder', 
                                       #                           class = "btn-dir", icon = icon("folder-open")),
                                       #verbatimTextOutput("dirPath"),
                                       div(style = "margin-top: 0px;",
                                           pickerInput("faFile", NULL, choices = NULL , multiple = FALSE,
                                                       options = list(style = "btn-default btn-md btn-picker", size = 20, 
                                                                      `live-search` = TRUE,  tickIcon = "glyphicon-ok", 
                                                                      width = "100%",showTick = TRUE)
                                           )
                                       ),
                                     ),
                                     conditionalPanel(
                                       condition = "input.faSrc == 'Remote'",
                                       selectizeInput("faSrc", NULL, choices = "DropBox")
                                     )
                                 ),
                                 div(style = "margin-top: 0px;",
                                     withBusyIndicatorUI(
                                       actionButton("import_fasta", "Import FASTA", class = "btn-action",
                                                    icon = icon('import', lib = "glyphicon"))
                                     )
                                 )
                             )
                         )
                         
                         
                         
                )
                #tags$hr()
                #p("Batch run: plot Hi-C maps of 20 most likley sequences to follow the selected sequence"),
                #actionButton("build_map", "Run", width = '100%', class = "btn-action")
            ),
            div(style = "margin-left: 250px;",   
                div(
                  actionButton("vwData", "Data View",icon = icon("table"), width = "285px", class= "slide-btn")
                ),
                fluidRow(style = "margin: 0 0 0 10px; border: 1px solid #E8E8E8; border-radius: 5px; padding: 1px; font-family: Tahoma; font-size: 12px;",
                         div(id = "DataView", style = "margin: 0px; display: none;",  
                             div(style = "width: 450px; float: left; background: #ece9df; border: 1px solid #E8E8E8; border-radius: 3px; padding: 5px; margin: 0 2px 0 10px;",
                                 DT::DTOutput('seqLen')
                             ),
                             div(style = "width: calc(100% - 580px); float: left; margin-left: 5px;",
                                 div(style = "border-left: 3px solid #E8E8E8; border-radius: 3px;",
                                     div(
                                       actionButton("vwStats", "Stats",icon = icon("stats", lib = "glyphicon"), width = "200px"),
                                       tags$style(HTML("#vwStats{font-size: 12px; background-color: #ece9df; border-top-right-radius: 30px; padding-top: 5px; border: 1px solid #E8E8E8;
                                               padding-bottom: 5px; font-weight: 550; margin-top: 0; margin-left: -5px;}"))
                                     ),
                                     fluidRow(id = "VwStats",  style = "margin: 1px 0 0 10px; border: 1px solid #E8E8E8; border-radius: 3px; display: none;",
                                              plotOutput("lenHist", width = "500px")
                                     ),
                                     
                                     div(
                                       actionButton("vwLnkData", "Contact Data",icon = icon("table"), width = "200px"),
                                       tags$style(HTML("#vwLnkData{font-size: 12px; background-color: #ece9df; border-top-right-radius: 30px; padding-top: 5px; border: 1px solid #E8E8E8;
                                               padding-bottom: 5px; font-weight: 550; margin-top: 5px; margin-left: -5px;}"))
                                     )
                                 ),
                                 fluidRow(id = "VwLnkData",  style = "padding: 5px; margin: 0 0 0 10px; border: 1px solid #E8E8E8; border-radius: 3px; display: none;",
                                          DT::DTOutput('contactData')
                                 )
                             )
                         )
                ),
                
                div(
                  actionButton("vwIntraMap", "Single Map View",icon = icon("table"), width = "285px", class = 'slide-btn')
                ),
                fluidRow(style = "margin: 0 0 0px 10px; border: 1px solid #E8E8E8; padding: 1px;",
                         div(id = "IntraMap", style = "height: calc(100% + 160px); display: none; ",
                             fluidRow(style = "margin: 0 0 2px 0; background-color: #ece9df; height: 32px;",
                                      radioGroupButtons("action", NULL, 
                                                        choices = c("View","Cut", "Join"), 
                                                        selected = "View", size = "sm", status = "gray",
                                                        width = "250px", justified = TRUE, individual = TRUE,
                                                        #checkIcon = list(
                                                        #  yes = tags$i(class = "fa fa-circle", 
                                                        #               style = "color: black"),
                                                        #  no = tags$i(class = "fa fa-circle-o", 
                                                        #              style = "color: white"))
                                      )
                             ),
                             div(style = "border: 1px solid #E8E8E8; width: 100%; padding: 5px; background-color: #ece9df; border-radius: 3px; font-size: 12px; width: 270px; float: left;",
                                 conditionalPanel(
                                   condition = "input.action == 'Join'",
                                   tags$button(id = "intConfig1", class = "btn btn-default btn-md btn-drop",
                                               list(NULL, label = "Inter-sequence Links"),
                                               htmltools::tags$span(class = "glyphicon glyphicon-triangle-bottom")
                                   ),
                                   div(style = "border: 1px solid #E8E8E8; padding: 1px; background-color: white; margin-bottom: 10px;",
                                       div(id ="IntConfig1", style = "padding: 10px;",
                                           fluidRow(style = "margin: 0px;",
                                                    numericInput("edgeSize1", "Edge Size (bp):", value = 1, 
                                                                 min = 0, max = 500000, step= 10000),
                                                    div(style = "margin-top: -13px;",
                                                        sliderInput("edgeSize2", NULL, value = 0, 
                                                                    min = 0, max = 0, width = '100%', ticks = TRUE)
                                                    ),
                                                    div(style = "margin-top: -30px;",
                                                        withBusyIndicatorUI(
                                                          actionButton("calcIntraction1", "Calculate", class = "btn-action",
                                                                       icon = icon('calculator', lib = "font-awesome"))
                                                        )
                                                    )
                                           )
                                       )
                                   ),
                                   
                                   conditionalPanel(
                                     condition = "input.calcIntraction1 != 0 | input.calcIntraction2 != 0",
                                     div(style = "border: 1px solid #E8E8E8; padding: 7px; background-color: white; margin-bottom: 0px; margin-top: -5px; border-radius: 5px;",
                                         radioGroupButtons("dir1", "Leading Sequence:", size = "xs",
                                                           choices = c(`<i class='glyphicon glyphicon-backward'></i>` = "Backward", 
                                                                       `<i class='glyphicon glyphicon-forward'></i>` = "Forward"), 
                                                           selected = "Forward", justified = TRUE, status = "gray",
                                                           checkIcon = list(
                                                             yes = tags$i(class = "fa fa-circle", 
                                                                          style = "color: black"),
                                                             no = tags$i(class = "fa fa-circle-o", 
                                                                         style = "color: white"))
                                         ),
                                         div(style = "margin-top: -13px;",
                                             pickerInput("seq1", NULL, choices = NULL, 
                                                         width = '100%', multiple = FALSE,
                                                         options = list(style = "btn-default btn-md btn-picker", 
                                                                        size = 10, `live-search` = TRUE)
                                             ),
                                             pickerInput("subseq1", "Subsequence:", choices = NULL, width = '100%', multiple = FALSE,
                                                         options = list(style = "btn-default btn-md btn-picker", 
                                                                        size = 10, `live-search` = TRUE)
                                             ),
                                             div(style = "margin-top: -15px;",
                                                 shinyWidgets::switchInput("subseq1_strand", label = "Strand", offLabel = "-", onLabel = "+", 
                                                                           value = TRUE, size = "mini", onStatus = "primary", 
                                                                           offStatus = "primary", labelWidth = '70', handleWidth = '70')
                                             ),
                                             withBusyIndicatorUI(
                                               actionButton("join", "Join", class = "btn-action",
                                                            icon = icon('plus-sign', lib = "glyphicon"))
                                             )
                                         )
                                     )
                                   )
                                 ),
                                 conditionalPanel(
                                   condition = "input.action != 'Join'",
                                   fluidRow(style = "border: 1px solid #E8E8E8; border-radius: 3px; padding: 10px 0 10px 10px; background-color: white; margin: 0;",
                                            div(style = "float: left; width: calc(100% - 70px); margin-top: -1px;",
                                                pickerInput("seq", "Sequence:", choices = "", 
                                                            width = '100%', multiple = FALSE,
                                                            options = list(style = "btn-default btn-md btn-picker", 
                                                                           size = 10, `live-search` = TRUE)
                                                )
                                            ),
                                            div(style = "float: right; width: 65px; margin: 20px 0 0 5px; padding: 1px; height: 33px;",
                                                actionBttn("down", NULL, icon = icon('arrow-down', lib = "glyphicon"), 
                                                           style = "material-circle",  size = "xs"),
                                                actionBttn("up", NULL, icon = icon('arrow-up', lib = "glyphicon"),
                                                           style = "material-circle", size = "xs")
                                            ),
                                            
                                            conditionalPanel(
                                              condition = "input.seq.length > 0",
                                              div(style = "margin-top: -5px; margin-right: 10px;",
                                                  conditionalPanel(
                                                    condition = "input.action == 'Cut'",
                                                    numericInput("cutPos", "Break position (bp):", value = 0, 
                                                                 min = 0, max = 50000000, step= 100),
                                                    textInput("frag1Name", "First fragmnet name:", value = "Default",  width = '100%'),
                                                    textInput("frag2Name", "Second fragmnet name:", value = "Default",  width = '100%'),
                                                    withBusyIndicatorUI(
                                                      actionButton("cut1", "Break Sequence", class = "btn-action", 
                                                                   icon = icon('scissors', lib = "glyphicon"))
                                                    )
                                                  )
                                              )
                                            )
                                   )
                                 ),
                                 conditionalPanel(
                                   condition = "input.join | input.cut1 | input.cut2 ",
                                   withBusyIndicatorUI(
                                     actionButton("saveBinData2", "Save changes to disk", width = '100%', class = "dbutt", 
                                                  icon = icon("floppy-disk", lib = "glyphicon")
                                     )
                                   )
                                 )
                             ),
                             div(style = "width: calc(100% - 270px); float: left;",
                                 div(style = "width: 82vh; border: 1px solid #E8E8E8; border-radius: 5px; padding: 10px; background: #ece9df; margin: 0px 5px 5px 5px; font-size: 12px;",
                                     conditionalPanel(
                                       condition = "input.action != 'Join'",
                                       actionBttn("exitZoom1", NULL, icon = icon('zoom-out', lib = "glyphicon"), style = "material-circle", size = "xs"),
                                       actionBttn("clearBrush1", NULL, icon = icon('erase', lib = "glyphicon"), style = "material-circle", size = "xs"),
                                       actionBttn("cut2", NULL, icon = icon('scissors', lib = "glyphicon"), style = "material-circle", size = "xs"),
                                       div(style = "border: 1px solid white; border-radius: 5px; margin: 0px; width: 80vh; margin-top: 5px;",
                                           #sliderInput("range", label = NULL, min = 0, max = 100, value = c(0, 100), ticks = FALSE, width = '100%'),
                                           plotOutput("intramap", width = "auto", height = "auto", 
                                                      dblclick = "intramap_dblclick",  
                                                      brush = "intramap_brush", 
                                                      click = "intramap_click") %>%
                                             withLoader(type = "html")
                                       )
                                     ),
                                     conditionalPanel(
                                       condition = "input.action == 'Join'",
                                       div(style = "border: 1px solid white; border-radius: 5px; margin: 0px; width: 80vh; margin-top: 5px;",
                                           #sliderInput("range", label = NULL, min = 0, max = 100, value = c(0, 100), ticks = FALSE, width = '100%'),
                                           plotOutput("joinedmap", width = "auto", height = "auto", 
                                                      brush = "joinedmap_brush", 
                                                      click = "joinedmap_click") %>%
                                             withLoader(type = "html")
                                       )
                                     )
                                 )
                             )
                         )
                ),
                div(
                  actionButton("vwJointMap2", "Multi Map View",icon = icon("table"), width = "285px", class = 'slide-btn'),
                ),
                fluidRow(style = "margin: 0 0 0 10px; border: 1px solid #E8E8E8; border-radius: 5px; padding: 1px; font-family: Tahoma; font-size: 12px;",
                         div(id = "VwJointMap2", style ="height: calc(100% + 160px); display: none;",
                             div(style = "border: 1px solid #E8E8E8; width: 100%; padding: 10px; background-color: #ece9df; border-radius: 3px; font-size: 12px; width: 270px; float: left;",
                                 
                                 textAreaInput("seqList", label = NULL, value = "", rows = 6, 
                                               placeholder = "Please enter one id per line", 
                                               width = '100%', height = '300px', resize = "vertical"),
                                 textInput("titleMap2",label = "Scaffold name:", value = "", width = '100%'),
                                 div(style = "margin-bottom: 20px;",
                                     withBusyIndicatorUI(
                                       actionButton("combineMaps", "Submit", icon = icon("play", lib = "glyphicon") ,width = '100%')
                                     ),
                                     tags$style(HTML("#combineMaps{font-size: 12px; font-family: Tahoma; font-weight: 550;}"))
                                 )
                                 
                             ),
                             div(style = "width: calc(100% - 270px); float: left;",
                                 div(style = "width: 82vh; border: 1px solid #E8E8E8; border-radius: 5px; padding: 10px; background: #ece9df; margin: 5px 5px 5px 5px; font-size: 12px;",
                                     fluidRow( style = "margin: 0;",
                                               div( style = "float: left;",
                                                    dropdownButton(
                                                      
                                                      
                                                      circle = TRUE, status = "default", size = "sm",
                                                      icon = icon("cog", lib = "glyphicon"), width = "300px",
                                                      
                                                      tooltip = tooltipOptions(title = "Click to edit plot")
                                                    )
                                               )
                                     ),
                                     fluidRow(style = "border: 1px solid white; border-radius: 5px; margin: 0px; width: 80vh;",
                                              plotOutput("map2", width = "auto", height = "auto", hover = hoverOpts(id = "map2_hover")) %>%
                                                withLoader(type="html")
                                     ),
                                     uiOutput("map2_hoverinfo")
                                 )
                             )
                         )
                )
            )
  )
)

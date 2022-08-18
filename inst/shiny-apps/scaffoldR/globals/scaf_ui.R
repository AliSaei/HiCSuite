shinyUI(
  fluidPage(style = "font-family: Tahoma;",
            div(style = "position: fixed; height: 100%; background: #ece9df; padding: 10px; width: 250px; border-right: 1px solid #E8E8E8; font-family: Tahoma; font-size: 12px;",
                fluidRow(style = "margin: 10px 0 0 0;",
                         tags$button(id = "intConfig2", class = "btn btn-default btn-md btn-accordion",
                                     list(NULL, label = "Trans Interactions"),
                                     htmltools::tags$span(class = "glyphicon glyphicon-triangle-bottom")
                         ),
                         div(style = "border: 1px solid #E8E8E8; border-radius: 3px; padding: 2px;  background-color: white; margin-bottom: 10px;",
                             div(id = "IntConfig2", style = "padding: 5px;",
                                 radioGroupButtons("lnkDataSrc", "Data source:", choices = c("Calculate", "Import"), 
                                                   selected = "Calculate", 
                                                   justified = TRUE, size = "xs", status = "gray", 
                                                   checkIcon = list(
                                                     yes = tags$i(class = "fa fa-circle", 
                                                                  style = "color: black"),
                                                     no = tags$i(class = "fa fa-circle-o", 
                                                                 style = "color: white"))),
                                 conditionalPanel(
                                   condition = "input.lnkDataSrc == 'Calculate'",
                                   fluidRow(style = "margin: -12px 0 0 0; border-top: 1px solid #E8E8E8; padding-top: 5px;",
                                            numericInput("edgeSize3", "Edge Size (bp):", value = 1, 
                                                         min = 0, max = 500000, step= 10000),
                                            div(style = "margin-top: -13px;",
                                                sliderInput("edgeSize4", NULL, value = 0,
                                                            min = 0,  max = 0, width = '100%', ticks = TRUE)
                                            ),
                                            div(style = "margin-top: -30px;",
                                                withBusyIndicatorUI(
                                                  actionButton("calcIntraction2", "Calculate", class = "btn-action",
                                                               icon = icon('calculator', lib = "font-awesome"))
                                                )
                                            ),
                                            conditionalPanel(
                                              condition = "input.calcIntraction2 != 0 | input.calcIntraction1 != 0",
                                              withBusyIndicatorUI(
                                                actionButton("svInteractionCounts", "Save to disk", width = '100%', class = "dbutt", 
                                                             icon = icon("floppy-disk", lib = "glyphicon")
                                                )
                                              )
                                            )
                                   )
                                 ),
                                 conditionalPanel(
                                   condition = "input.lnkDataSrc == 'Import'",
                                   div(style = "margin-top: -12px;",
                                       pickerInput("lnkFile", NULL, choices = NULL , multiple = FALSE,
                                                   options = list(style = "btn-default btn-md btn-picker", size = 20, 
                                                                  `live-search` = TRUE,  tickIcon = "glyphicon-ok", 
                                                                  width = "100%", showTick = TRUE)
                                       ),
                                       div(style = "margin-top: 0px;",
                                           withBusyIndicatorUI(
                                             actionButton("importLnkFile", "Import", class = "btn-action",
                                                          icon = icon('import', lib = "glyphicon"))
                                           )
                                       )
                                   )
                                 )
                             )
                         ),
                         
                         
                         
                         tags$button(id = "scafConfig", class = "btn btn-default btn-md btn-accordion",
                                     list(NULL, label = "Scaffolding"),
                                     htmltools::tags$span(class = "glyphicon glyphicon-triangle-bottom")
                         ),
                         div(style = "border: 1px solid #E8E8E8; border-radius: 5px; padding: 2px;  background-color: white; margin-top: 0px; ",
                             div(id = "ScafConfig", style = "padding: 5px; display: none;",
                                 div(
                                   fluidRow(style = "margin: 0px;",
                                            radioGroupButtons("dir2", "Leading sequence:", size = "xs",
                                                              choices = c(`<i class='glyphicon glyphicon-backward'></i>` = "Backward", 
                                                                          `<i class='glyphicon glyphicon-forward'></i>` = "Forward"), 
                                                              selected = "Forward", justified = TRUE, status = "gray",
                                                              checkIcon = list(
                                                                yes = tags$i(class = "fa fa-circle", 
                                                                             style = "color: black"),
                                                                no = tags$i(class = "fa fa-circle-o", 
                                                                            style = "color: white")))
                                   ),
                                   fluidRow(style = "margin: -15px 0 0 0;", 
                                            pickerInput("seq2", NULL, choices = "", width = '100%', multiple = FALSE,
                                                        options = list(style = "btn-default btn-md btn-picker", 
                                                                       size = 10, `live-search` = TRUE)),
                                            
                                            div(style = "margin-top: -15px;",
                                                shinyWidgets::switchInput("strand_2", label = "Strand", 
                                                                          offLabel = "-", onLabel = "+", 
                                                                          value = TRUE, size = "mini", onStatus = "primary", 
                                                                          offStatus = "primary", labelWidth = '100', 
                                                                          handleWidth = '100')
                                            )
                                   ),
                                   fluidRow(style = "margin: 0px;", 
                                            div(style = "float: left; width: calc(100% - 70px); margin-top: -1px;",
                                                pickerInput("subseq2", "Subsequent Seq:", choices = NULL, width = '100%',
                                                            options = list(style = "btn-default btn-md btn-picker", size = 10, `live-search` = TRUE)
                                                ),
                                                div(style = "margin-top: -15px;",
                                                    shinyWidgets::switchInput("strand_3", label = "Strand", offLabel = "-", onLabel = "+", 
                                                                              value = TRUE, size = "mini", onStatus = "primary", 
                                                                              offStatus = "primary", labelWidth = '70', handleWidth = '70')
                                                )
                                            ),
                                            div(style = "float: right; width: 70px; margin: 20px 0 0 0; padding: 1px 1px 1px 5px; height: 33px;",
                                                actionBttn("down1", NULL, icon = icon('arrow-down', lib = "glyphicon"), 
                                                           style = "material-circle",  color = "default", size = "xs"),
                                                actionBttn("up1", NULL, icon = icon('arrow-up', lib = "glyphicon"),
                                                           style = "material-circle",  color = "default", size = "xs")
                                            )
                                   ),
                                   p("Filters:", style = "font-weight: 550;"),
                                   fluidRow(style = "margin: -5px 0 0 0; border: 1px solid #E8E8E8; padding:2px; border-radius: 3px; background: #E8E8E8;", 
                                            div(style = "margin: 0;",
                                                sliderInput("maxLinkDen", NULL, value = 1, min = 0, max = 1, 
                                                            width = '100%', ticks = FALSE, pre = "Max link density:")
                                            ),
                                            div(style = "margin-top: -10px;",
                                                sliderInput("minSeqLen", NULL, value = 0, min = 0, max = 1, 
                                                            width = '100%', ticks = FALSE, pre = "Min length: ", post = " Mb")
                                            )
                                   ),
                                   #textInput("chr",label = "Assembled sequence:", value = "", width = '100%'),
                                   #tags$hr(), 
                                   div(style = "margin-top: 10px;",
                                       withBusyIndicatorUI(
                                         actionButton("add", "Add to scaffold", class = "btn-action", 
                                                      icon = icon("plus", lib = "font-awesome"))
                                       ),
                                       withBusyIndicatorUI(
                                         actionButton("excludeSeq", "Exclude Sequence", class = "btn-action", 
                                                      icon = icon("xmark"))
                                       )
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
                  actionButton("vwData1", "Data View",icon = icon("table"), width = "285px", class= "slide-btn")
                ),
                fluidRow(style = "margin: 0 0 0 10px; border: 1px solid #E8E8E8; border-radius: 5px; padding: 1px; font-family: Tahoma; font-size: 12px;",
                         div(id = "DataView1", style = "margin: 0px; display: none;",  
                             div(style = "background: #ece9df; border: 1px solid #E8E8E8; border-radius: 3px; padding: 5px; margin: 0 2px 0 10px;",
                                 DT::DTOutput("interactionCounts") ,
                                 downloadButton("interactionCountsExp", "Download All Data")
                             )
                         )
                ),
                div(
                  actionButton("vwScaf", "Scaffolding View",icon = icon("table"), width = "285px", class= "slide-btn")
                ),
                fluidRow(style = "margin: 0 0 0 10px; border: 1px solid #E8E8E8; border-radius: 5px; padding: 1px; font-family: Tahoma; font-size: 12px; ",
                         div(id = "VwScaf", style ="height: calc(100% + 160px); display: none;",
                             div(style = "border: 1px solid #E8E8E8; padding: 5px; background-color: #ece9df; border-radius: 3px; font-size: 12px; width: 300px; float: left;",
                                 #div("Setting", icon("triangle-bottom", class = "glyph-icon", lib = "glyphicon")),
                                 
                                 fluidRow(style = "margin: 0px;",
                                          div(style = "width: 35px; float: left; padding:5px;",
                                              actionBttn("edit", NULL, icon = icon('edit', lib = "glyphicon"), 
                                                         style = "material-circle", size = "xs"),
                                              shinyjs::hidden(
                                                actionBttn("confirm", NULL, icon = icon('check', lib = "glyphicon"), 
                                                           style = "material-circle", size = "xs")
                                              ),
                                              actionBttn("erase", NULL, icon = icon('erase', lib = "glyphicon"), 
                                                         style = "material-circle", size = "xs"),
                                              actionBttn("clipbtn", NULL, icon = icon('copy', lib = "glyphicon"), 
                                                         style = "material-circle", size = "xs"),
                                              actionBttn("reverse", NULL, icon = icon('sort-by-alphabet-alt', lib = "glyphicon"), 
                                                         style = "material-circle", size = "xs"),
                                              actionBttn("export", NULL, icon = icon('save', lib = "glyphicon"), 
                                                         style = "material-circle", size = "xs")
                                          ),
                                          
                                          div(style = "width: calc(100% - 35px); float: left; padding: 2px;",
                                              div(id = "CheckBox", style = "border: 1px solid #F4F4F4; padding: 5px; border-radius: 3px; min-height: 150px; max-height: 550px; overflow: auto; background: white; margin-right: 2px; resize: vertical; overflow: auto;",
                                                  checkboxGroupInput("anchored_seqs", "Scaffold:", choices = NULL)
                                              ),
                                              hidden(
                                                div(id = "EditBox" , style = "margin-top: 0px;",
                                                    textAreaInput("scaf_edit", label = NULL, value = "", 
                                                                  placeholder = "Please enter one id per line", 
                                                                  width = '100%', height = '100%', resize = "vertical")
                                                )
                                              )
                                          )
                                 ),
                                 conditionalPanel(
                                   condition = "input.inputType == 'Manual'",
                                   div(style = "margin-top: -15px;",
                                       textAreaInput("scaf_man", label = NULL, value = "", 
                                                     placeholder = "Please enter one id per line", 
                                                     width = '100%', height = "300px", resize = "vertical")
                                   )
                                 ),
                                 div(
                                   withBusyIndicatorUI(
                                     textAreaInput("excluded_seqs", label = "Excluded sequences:", value = "", 
                                                   placeholder = "please enter one id per line", 
                                                   width = '100%', height = '200px', resize = "vertical")
                                   )
                                 )
                                 
                             ),
                             div(style = "margin-left: 5px; float: left; width: calc(100% - 320px);",  
                                 div(style = "border-left: 3px solid #E8E8E8;",
                                     div(
                                       actionButton("vwJointMap", "Scaffolding Map",icon = icon("table"), width = "200px"),
                                       tags$style(HTML("#vwJointMap{font-size: 12px; background-color:#ece9df; border-top-right-radius: 30px; padding-top: 5px; border: 1px solid #E8E8E8;
                                               padding-bottom: 5px; font-weight: 550; margin-top: 0px; margin-left: -5px;}"))
                                     ),
                                     fluidRow(id = "VwJointMap", style = "margin: 0 0 0 10px; border: 1px solid #E8E8E8; border-radius: 5px; padding: 1px; font-family: Tahoma; font-size: 12px;",
                                              div(style = "font-size: 12px; font-family: Tahoma; float: left;",
                                                  fluidRow(style = "margin: 2px 5px 0 5px; background: #ece9df; padding: 1px; border: 1px solid #F4F4F4; height: 37px; border-radius: 5px;",
                                                           div("Number of sequences to plot from the leading end:", style = "width: 300px; float: left; padding: 7px;"),
                                                           div(style = "width: calc(100% - 500px); float: left;",
                                                               numericInput("nrSeq", NULL, value = 1, min = 1, max = 10, step= 1, width = '100%')
                                                           ),
                                                           div(style = "float: right; ",
                                                               actionButton("refresh", NULL, icon = icon('refresh', lib = "glyphicon"), width = '198px')
                                                           ),
                                                  ),
                                                  div(style = "width: 82vh; border: 1px solid #E8E8E8; border-radius: 5px; padding: 10px; background: #ece9df; margin: 5px 5px 5px 5px; font-size: 12px;",
                                                      div(
                                                        actionBttn("exitZoom2", NULL, icon = icon('zoom-out', lib = "glyphicon"), 
                                                                   style = "material-circle", size = "xs"),
                                                        actionBttn("clearBrush2", NULL, icon = icon('erase', lib = "glyphicon"), 
                                                                   style = "material-circle", size = "xs")
                                                      ),
                                                      div(style = "border: 1px solid white; border-radius: 5px; margin: 0px; width: 80vh; margin-top: 5px;",
                                                          plotOutput("map1", width = "auto", height = "auto", 
                                                                     brush = "map1_brush",
                                                                     hover = hoverOpts(id = "map1_hover")) %>%
                                                            withLoader(type = "html")
                                                      )
                                                  )
                                              ),
                                              div(style = "width: calc(100% - 84vh); float: left; margin-left: 5px;",
                                                  div(style = "border-left: 3px solid #E8E8E8;",
                                                      div(
                                                        actionButton("vwStats2", "Stats",icon = icon("stats", lib = "glyphicon"), width = "200px"),
                                                        tags$style(HTML("#vwStats2{font-size: 12px; background-color: #ece9df; border-top-right-radius: 30px; padding-top: 5px; border: 1px solid #E8E8E8;
                                                                                padding-bottom: 5px; font-weight: 550; margin-top: 0; margin-left: -5px;}"))
                                                      ),
                                                      fluidRow(id = "VwStats2",  style = "margin: 1px 0 0 10px; border: 1px solid #E8E8E8; border-radius: 3px; display: none;"
                                                               
                                                      ),
                                                      div(
                                                        actionButton("vwLnkData2", "Contact Data", icon = icon("table"), width = "200px"),
                                                        tags$style(HTML("#vwLnkData2{font-size: 12px; background-color: #ece9df; border-top-right-radius: 30px; padding-top: 5px; border: 1px solid #E8E8E8;
                                                                                      padding-bottom: 5px; font-weight: 550; margin-top: 5px; margin-left: -5px;}"))
                                                      )
                                                  ),
                                                  fluidRow(id = "VwLnkData2",  style = "padding: 5px; margin: 0 0 0 10px; border: 1px solid #E8E8E8; border-radius: 3px;",
                                                           DTOutput("Subsequent")
                                                  )
                                              )
                                     ),
                                     div(
                                       actionButton("vwJointMap3", "Intra-scaffold Map",icon = icon("table"), width = "200px"),
                                       tags$style(HTML("#vwJointMap2{font-size: 12px; background-color:#ece9df; border-top-right-radius: 30px; padding-top: 5px;  border: 1px solid #E8E8E8;
                                                    padding-bottom: 5px; font-weight: 550; margin-top: 5px; margin-left: -5px;}"))
                                     )
                                 ),
                                 fluidRow(style = "margin: 0 0 0 10px; border: 1px solid #E8E8E8; border-radius: 5px; padding: 1px; font-family: Tahoma; font-size: 12px;",
                                          div(id = "VwJointMap3", style ="height: calc(100% + 160px); display: none;",
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
                                                      )
                                                     # fluidRow(style = "border: 1px solid white; border-radius: 5px; margin: 0px; width: 80vh;"
                                                        #  plotOutput("map2", width = "auto", height = "auto", hover = hoverOpts(id = "map2_hover")) %>%
                                                        #    withLoader(type="html")
                                                     # )
                                                      #uiOutput("map2_hoverinfo")
                                                  )
                                              )
                                          )
                                 )
                             )
                         )
                )
            )
  )
)

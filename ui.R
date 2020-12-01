source("./globals/functions.R")
source("./globals/helpers.R")
shinyUI(fluidPage(
  div(style = "position: fixed; height: 100%; background: #ece9df; padding-top: 10px; font-size: 12px; width: 350px; padding: 15px; border-right: 1px solid lightgray;",
      fluidRow(style = "margin: 0;",
               radioGroupButtons("fileSrc", "Upload data:", choices = c("Browse", "Path"), selected = "Path", justified = TRUE, size = "sm",
                                 checkIcon = list(
                                   yes = tags$i(class = "fa fa-circle", 
                                                style = "color: steelblue"),
                                   no = tags$i(class = "fa fa-circle-o", 
                                               style = "color: steelblue"))),
               div(style = "margin-top: -10px;",
                     conditionalPanel(
                       condition = "input.fileSrc == 'Browse'",
                       shiny::fileInput("contactMap", "Contact map:", accept = c("csv", ".txt", ".rds"), width = '100%'),
                       div(style = "margin-top: -10px;",
                       shiny::fileInput("seqLen", "Sequence length:", accept = c("csv", ".txt", ".rds"), width = '100%')
                       )
                     ),
                   conditionalPanel(
                     condition = "input.fileSrc == 'Path'",
                     pickerInput("mapPath", NULL, choices = c("./data/gillenia/Gillenia_HiC_mat_5k_scaffolds_filt.rds", 
                                                              "./data/gillenia/Gillenia_contact.map_contigs_filt.rds",
                                                              "./data/gillenia/mapped_ragoo_contact_map_1k_qualFilt.rds",
                                                              "./data/bilberry/bilbery_mat.10kb.rds",
                                                              "./data/bilberry/mapped_shasta_i1_contact_map_1k_qualFilt.rds",
                                                              "./data/bilberry/mapped_shasta_i3_contact_map_1k_qualFilt.rds",
                                                              "./data/gillenia/mapped_ragoo_ml10k_contact_map_1k_qualFilt.rds"), width = '100%', multiple = FALSE,
                                 options = list(size = 10, `selected-text-format` = "count > 1", `live-search` = TRUE)),
                     #textInput("mapPath", NULL, placeholder = "./path/to/contact_matrix_file", value = "./data/gillenia/Gillenia_HiC_mat_5k_scaffolds_filt.rds", width = "100%"),
                     div(style = "margin-top: -10px;",
                         pickerInput("lenPath", NULL, choices = c("./data/gillenia/contig_len.csv",
                                                                  "./data/gillenia/scaff_ragoo_gapfilled_noContamination_len.txt",
                                                                  "./data/bilberry/sequence_len_bilberry.csv",
                                                                  "./data/bilberry/Shasta_racon_i1_pilon-2.noBacteria.ml1000.fasta.csv",
                                                                  "./data/bilberry/Shasta_racon_i3_pilon-2.noBacteria.ml1000.fasta.csv",
                                                                  "./data/gillenia/scaff_ragoo_gapfilled_noContamination.ml10000.fasta.fai.txt"), width = '100%', multiple = FALSE,
                                     options = list(size = 10, `selected-text-format` = "count > 1", `live-search` = TRUE))
                     #textInput("lenPath", NULL, placeholder = "./path/to/sequence_length_file", value = "./data/gillenia/contig_len.csv", width = "100%")
                     )
                   )
               ),
               div(style = "margin-top: -10px;",
                   withBusyIndicatorUI(
                   actionButton("import_data", "Import data", icon = icon('import', lib = "glyphicon"), class = "action_button", width = '100%')
                   )
                   #checkboxInput("reverse_seq.1", "Reverse initial sequence")
               ),
               div(style = "margin-top: 10px;  margin-bottom: -10px;",
                   sliderInput("binsize", "Bin Size:", min = 0, value = 0, max = 0, width = '100%', ticks = TRUE)
               ),
           
               div(style = "float: left; width: 225px;",
                   pickerInput("seq_1", "Leading sequence:", choices = NULL, width = '100%', multiple = FALSE,
                               options = list(size = 10, `selected-text-format` = "count > 1", `live-search` = TRUE)),
                   #selectInput("seq_1", "Leading sequence:", choices = NULL, width = '100%'),
                   div(style = "margin-top: -15px;",
                       shinyWidgets::switchInput("orient_1", label = "Strand", offLabel = "-", onLabel = "+", value = TRUE, size = "mini", onStatus = "primary", offStatus = "primary")
                       #checkboxInput("reverse_seq.1", "Reverse initial sequence")
                   )
               ),
               div(style = "float: left; width: 90px;",
                   numericInput("nrSeq", "Number:", value = 1, min = 1, max = 10, step= 1, width = '100%')
               )
      ),

      fluidRow(style = "margin: 0;",
               p("Subsequent sequence:", style = "font-weight: 550;"),
               div(style = "float: left; width: 225px;",
                   pickerInput("seq_2", NULL, choices = NULL, width = '100%',
                               options = list(size = 10, `live-search` = TRUE)),
                   #selectInput("seq_2", NULL, choices = NULL, width = '100%'),
                   div(style = "margin-top: -15px;",
                       shinyWidgets::switchInput("orient_2", label = "Strand", offLabel = "-", onLabel = "+", value = TRUE, size = "mini", onStatus = "primary", offStatus = "primary")
                       #checkboxInput("reverse_seq.1", "Reverse initial sequence")
                   )
                   ),
               div(style = "float: right; width:90px; margin: 0px;",
                   actionButton("down", NULL, icon = icon('arrow-down', lib = "glyphicon"), class = "default"),
                   actionButton("up", NULL, icon = icon('arrow-up', lib = "glyphicon"))
                   #checkboxInput("reverse_seq.1", "Reverse initial sequence")
               )
      ),
      div(style = "margin-top: -5px;",
          radioGroupButtons("dir", "Scaffolding direction:", choices =  c("Forward", "Backward"), selected = "Forward", justified = TRUE, size = "sm",
                            checkIcon = list(
                              yes = tags$i(class = "fa fa-circle", 
                                           style = "color: steelblue"),
                              no = tags$i(class = "fa fa-circle-o", 
                                          style = "color: steelblue")))
      ),
      #textInput("chr",label = "Assembled sequence:", value = "", width = '100%'),
      #tags$hr(), 
      div(style = "margin-top: 10px;",
          withBusyIndicatorUI(
          actionButton("build_scaffold", "Combine sequences", width = '100%', class = "action_button")
          )
      ),
      div(style = "margin-top: 10px; margin-bottom: 10px;",
          fluidRow(style = "margin: 0px;",
                   div("Combined sequences:", style = "font-weight: 550; border: 1px solid #F0F0F0; border-radius: 3px; padding: 5px; width: 150px; background: #F0F0F0; float: left;"),
                   div(style = "border: 1px solid #F0F0F0; float: right; height: 29px; border-radius: 3px; padding: 3px 0 0 3px;  width: 168px;  background: #F0F0F0;",
                       radioButtons("inputType", NULL, choices = c("Auto", "Manual"), inline = TRUE)
                   )
          ),
          conditionalPanel(
            condition = "input.inputType == 'Auto'",
            div(style = "border: 1px solid #F0F0F0; padding: 5px; border-radius: 3px; margin-top: 0px; max-height: 300px; overflow-y: auto; background: white;",
                checkboxGroupInput("scaf_auto", NULL, choices = NULL)
            )
          ),
          conditionalPanel(
            condition = "input.inputType == 'Manual'",
            textAreaInput("scaf_man", label = NULL, value = "", placeholder = "please enter one id per line", width = '100%', height = "300px")
          ),
          conditionalPanel(
            condition = "input.inputType == 'Edit'",
            textAreaInput("scaf_edit", label = NULL, value = "", placeholder = "please enter one id per line", width = '100%', height = "300px")
          )
      ),
      div(
      withBusyIndicatorUI(
      actionButton("combineMaps", "Build intra chromosomal map", width = '100%', class = "action_button")
      )
      )
      #tags$hr()
      #p("Batch run: plot Hi-C maps of 20 most likley sequences to follow the selected sequence"),
      #actionButton("build_map", "Run", width = '100%', class = "action_button")
  ),
  div(style = "margin-left: 350px; position: relative;",   
      div(
        actionButton("vwIntarChr", "Intra-scaffold contact map",icon = icon("table"), width = "75%"),
        tags$style(HTML("#vwIntarChr{font-size: 12px; background-color:#ece9df; border-top-right-radius: 20px; padding-top: 5px; padding-bottom: 5px; font-weight: 550; margin-top: 2px;}"))
      ),
      
      #div(style = "background: #C0C0C0;  width: 75%; height: 30px; border-top-right-radius: 20px;margin-top: 2px;"),
      div(
        actionButton("vwMap", "Combined sequences contact map", icon = icon("table"), width = "24%"),
        tags$style(HTML("#vwMap{font-size: 12px; background-color:#ece9df; border-top-right-radius: 20px; padding-top: 5px; padding-bottom: 5px; font-weight: 550; margin-top: 2px;}"))
      ),
      fluidRow(style = "margin: 0px;",
               div(id = "Map1", class = "col-sm-6", style = "float: left;",
                   div(style = "width: 100%; border: 1px solid #F0F0F0; border-radius: 5px; padding: 10px; background: #ece9df;margin-left: 1%; position: relative;",
                   fluidRow(style = "margin: 0px; height: 30px; padding: 0px 0px 5px 5px;"
                            #checkboxInput("disable_plt", "Disable plot")
                   ),
                   fluidRow(style = "border: 1px solid white; border-radius: 5px; margin: 0px;",
                       plotOutput("map", height = "auto", dblclick = "map_dblclick", brush = "map_brush")
                   )
                   )
               ),
               div(id = "Map2", class = "col-sm-6", style = "float: right;",
                   div(style = "width: 100%; border: 1px solid #F0F0F0; border-radius: 5px; padding: 10px; background: #ece9df; margin-top: -31px; position: relative;",
                   fluidRow(style = "margin: 0px; height: 30px; padding: 0px 0px 5px 5px;"
                   ),
                   
                   # uiOutput("xslider"),
                   div(style = "border: 1px solid white; border-radius: 5px; position: relative;",
                       plotOutput("combined_maps", height = "auto")
                   )
                   )
               )
      )
  ),
  shinyjs::useShinyjs(),
  tags$style(appCSS),
  tags$script(src="jscript.js"),
  tags$head(tags$style(HTML("
         .action_button{font-size: 12px; font-family:Tahoma; width: 100%; font-weight: 550; background-color: green;}
                            .selectize-input { font-size: 11px; font-family: Tahoma;}
                            .selectize-dropdown { font-size: 12px; font-family: Tahoma; overflow-x: auto;}
                            .selectize-dropdown-content {
                            max-height: 300px;
                            padding: 0;
                            }
                            .selectize-dropdown-content .option {
                            border-bottom: 1px dotted #ccc;
                            }
                            .fa {font-size: 14px;}
                            .glyph-icon {font-size: 11px;}
                            
                            ")))
  
)
)

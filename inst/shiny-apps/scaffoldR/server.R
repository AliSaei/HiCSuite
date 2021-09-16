library(shiny)
library(data.table)
library(dplyr)
library(ggplot2)
library(shinyWidgets)
library(shinycustomloader)
library(shinyjs)
library(shinyFiles)

options(future.globals.maxSize = 300*1024^2,
        shiny.maxRequestSize=300*1024^2)

#negate %in% 
`%!in%` = Negate(`%in%`)

#source("../../../R/join_maps.R")
#source("../../../R/join_maps_plus.R")
source("./functions.R")
source("./helpers.R")

#library(readxl)
#anchored_contigs <- read_excel("/Volumes/workspace/hrpazs/bilberry_genome/chr_bilberry.xlsx", 
#                               sheet =  "Billberry (2)") %>%
#  as.data.table(.) %>%
#  data.table::melt("No.", measure = 2:ncol(.), variable.name = "chr_no", value.name = "name", na.rm = TRUE) %>%
#  .[, ':='(name = substr(name, 2, length(name)),
#           rc=ifelse(grepl("\\+", name), 0, 1),
#           q = ".", gap_size = ".")]

#groups <- as.character(unique(anchored_contigs$chr_no))
#out_dir <- "/Volumes/workspace/hrpazs/bilberry_genome/groups/"
#dir.create(out_dir)

#for (i in groups){
#  data <- anchored_contigs %>%
#    .[chr_no == i, .(name, rc, q , gap_size)]
#  file = paste0(out_dir, i,".ordering")
#  write.table(data, file, row.names = FALSE, col.names = FALSE, quote = FALSE) 
#}



server <- function(input, output, session) {
  shinyOptions(progress.style = "old")
  rv <- reactiveValues(chr = NULL, s2.1 = NULL, cut_pos = 0, 
                       intramap_range = NULL,  btn_val2 = c(0,0))
  
  ## shiny files set up
  volumes <- c(Home = fs::path_home(), getVolumes()(), `Test Data Directory` = "../../../")
  
  shinyDirChoose(input, 'directory', roots = volumes, session = session,
                 restrictions = system.file(package = "base"))
  
  observe({
    rv$projDir <- parseDirPath(volumes, input$directory)
    
    rds_list <- list.files(rv$projDir, pattern = ".rds$", recursive = TRUE)
    txt_list <- list.files(rv$projDir, pattern = ".csv$|.txt$", recursive = TRUE)
    
    updatePickerInput(session, "mapFile", choices = rds_list,                           
                      choicesOpt = list(style = paste(rep_len("font-size: 12px; line-height: 1.5; 
                                                              margin-left: -10px; border-bottom: 1px solid gray;", 
                                                              length(rds_list)), "background-color: lightgray;")))
    updatePickerInput(session, "lenFile", choices = txt_list,                           
                      choicesOpt = list(style = paste(rep_len("font-size: 12px; line-height: 1.5; 
                                                              margin-left: -10px; border-bottom: 1px solid gray;", 
                                                              length(txt_list)), "background-color: lightgray;")))
    updatePickerInput(session, "lnkFile", choices = txt_list,                           
                      choicesOpt = list(style = paste(rep_len("font-size: 12px; line-height: 1.5; 
                                                              margin-left: -10px; border-bottom: 1px solid gray;", 
                                                              length(txt_list)), "background-color: lightgray;")))
  })
  
  
  ## Import contact map genearted using the jupyternote book script from BAM file
  observeEvent(input$import_data,{
    withBusyIndicatorServer("import_data", {
      withProgress(message = 'Loading in progress',
                   detail = 'Please wait...', value = 1,{
                     #updatePickerInput(session, "seq1", selected = "")
                     
                     rv$mat <- as.data.table(readRDS(file.path(rv$projDir,input$mapFile))) 
                     
                     ## calculate or import length of sequences
                     if(input$loadSeqLen){
                       rv$sequence_length <- fread(file.path(rv$projDir,input$lenFile), 
                                                   col.names = c("rname","rname_len"), select = c(1, 2)) 
                     } else {
                       rv$sequence_length <- rv$mat[,.(rname_len = max(pos)), by = .(rname)]
                     }
                     
                     rmax = max(rv$mat$pos)
                     mmax = max(rv$mat$mpos)
                     rv$max = max(rmax, mmax)
                     binsize_ini <- sort.int(unique(rv$mat$pos), partial = 1:2)[1:2]
                     binsize_ini = as.integer(binsize_ini[2] - binsize_ini[1])
                     choices <- rv$sequence_length[rname %in% unique(rv$mat[, rname]), ][order(-rname_len), rname]
                     choices_opt <- list(style = paste(rep_len("font-size: 12px; line-height: 1.5; margin-left: -10px; 
                                                        border-bottom: 1px solid gray;", length(choices)), 
                                                       "background-color: lightgray;"))
                     
                     updateNumericInput(session, "binsize2", value = binsize_ini, min = binsize_ini,  
                                        step = binsize_ini , max = 500000)
                     
                     updateNumericInput(session, "edgeSize1",  value = 20 * binsize_ini, min = binsize_ini, 
                                        step = binsize_ini , max = 500000)
                     
                     updateNumericInput(session, "edgeSize3",  value = 20 * binsize_ini, min = binsize_ini, 
                                        step = binsize_ini , max = 500000)
                     
                     updatePickerInput(session, "seq", choices = choices, choicesOpt = choices_opt)
                     updatePickerInput(session, "seq1", choices = choices, choicesOpt = choices_opt)
                     updatePickerInput(session, "seq2", choices = c( "", choices), 
                                       selected = input$seq2, choicesOpt = choices_opt)
                     
                     rv$binsize_ini <- binsize_ini
                     rv$binsize <- binsize_ini
                     rv$mat_binned <- rv$mat
                     rv$seq_name <- choices
                   })
    })
  })
  
  observe({
    updateSliderInput(session, "binsize", value = input$binsize2, step = rv$binsize_ini , max = input$binsize2)
  })
  
  ##------------------------------------------------------------------------------------------------------------
  ##------------------------------------------------------------------------------------------------------------
  # udpate bin size
  observeEvent(input$update_bin, {
    withBusyIndicatorServer("update_bin", {
      withProgress(message = 'Binning in progress',
                   detail = 'please wait ...', value = 1, {
                     rv$binsize <- input$binsize
                     
                     if(rv$binsize_ini != rv$binsize){
                       bins =  seq(rv$binsize , rv$max, rv$binsize)
                       breaks =  c(0, bins, (length(bins) + 1) * rv$binsize)
                       
                       rv$mat_binned <- rv$mat %>%
                         .[, ':='(rbin = cut(pos,  breaks, include.lowest = TRUE),
                                  mbin = cut(mpos,  breaks, include.lowest = TRUE))] %>%
                         .[, .(n = sum(n)), by = .(rname, mrnm, rbin, mbin)] %>%
                         .[, pos :=  as.numeric(substr(rbin, 2, gregexpr(",", rbin)[[1]][1]-1)), by = "rbin"] %>%
                         .[,mpos :=  as.numeric(substr(mbin, 2, gregexpr(",", mbin)[[1]][1]-1)), by = "mbin"]
                     } else {
                       rv$mat_binned <- rv$mat
                     }
                   })
    })
  })
  
  
  ##------------------------------------------------------------------------------------------------------------
  ##------------------------------------------------------------------------------------------------------------  
  # click events on intar-map view
  observe({
    shiny::validate(need(rv$mat_binned, ""))
    if(is.null(input$intramap_dblclick)) 
      return(NULL)
    
    xy_str <- function(e) {
      if(is.null(e)) return(NULL)
      round(e$x/rv$binsize, 0) * rv$binsize
    }
    
    updateRadioGroupButtons(session, 'action', selected = "Cut")
    updateNumericInput(session, "cutPos", value = xy_str(input$intramap_dblclick))
  })
  
  
  # brush event
  observe({
    shiny::validate(need(rv$binsize > 0, ""))
    
    if(is.null(input$intramap_brush)) return(NULL)
    
    xy_range_str <- function(e) {
      if(is.null(e)) return(NULL)
      xmin <- round(e$xmin, 1)
      range <- c(ifelse(xmin < 0, 0, round(xmin/rv$binsize, 0) * rv$binsize),
                 round(e$xmax/rv$binsize, 0) * rv$binsize)
      return(range)
    }
    rv$intramap_range <- xy_range_str(input$intramap_brush)
  })
  
  
  observe(rv$cut_pos <- input$cutPos)
  
  ##--- cut sequence
  observeEvent(input$cut1 + input$cut2,{
    shiny::validate(need(input$cutPos > 0, ""))
    
    withBusyIndicatorServer("cut1", {
      withProgress(message = 'Cutting the sequence',
                   detail = 'please wait ...', value = 1,{
                     
                     tgt_contig <- input$seq
                     tgt_len <- rv$sequence_length[rname == tgt_contig, rname_len]
                     subseq_start <- rv$cut_pos + rv$binsize
                     subseq_contig <- paste0(tgt_contig, "_subseq_", subseq_start, ":", tgt_len)
                     
                     rv$mat_binned <- rv$mat_binned %>%
                       #.[(rname %in% tgt_contig | mrnm %in% tgt_contig),] %>%
                       .[,':='(rname = ifelse(rname == tgt_contig & pos >  rv$cut_pos , subseq_contig, rname),
                               mrnm = ifelse(mrnm == tgt_contig & mpos >  rv$cut_pos , subseq_contig, mrnm))] %>%
                       .[, ':='(pos = ifelse(rname == subseq_contig,  pos - subseq_start, pos ),
                                mpos = ifelse(mrnm == subseq_contig , mpos - subseq_start , mpos))]
                     
                     rv$sequence_length <- rv$sequence_length[rname == tgt_contig, rname_len := rv$cut_pos] %>%
                       rbind(.,list(subseq_contig, tgt_len - rv$cut_pos))
                     
                     updateNumericInput(session, "cutPos", value = 0)
                     
                     choices <- c(rv$seq_name, subseq_contig)
                     updatePickerInput(session, "seq", choices = choices, selected = tgt_contig,
                                       choicesOpt = list(style = paste(rep_len("font-size: 12px; line-height: 1.5; 
                                                                               margin-left: -10px; border-bottom: 1px solid gray;", 
                                                                               length(choices)), "background-color: lightgray;")))
                     #set plot size to default
                     rv$intramap_range <- NULL
                   })
    })
  })
  
  
  observe({
    updateSliderInput(session, "edgeSize2", value = input$edgeSize1, 
                      step = rv$binsize_ini , max = input$edgeSize1)
  })
  
  observe({
    updateSliderInput(session, "edgeSize4", value = input$edgeSize3, 
                      step = rv$binsize_ini , max = input$edgeSize3)
  })
  
  
  ## this allows user to calculate inter-sequence links
  observeEvent(input$calcIntraction1 + input$calcIntraction2, {
    shiny::validate(need(rv$mat_binned, ""))
    
    btn_id <- c("calcIntraction1" , "calcIntraction2")
    btn_val <- c(input$calcIntraction1, input$calcIntraction2) - rv$btn_val2 
    i <- which(btn_val == max(btn_val))
    rv$btn_val2 <- btn_val
    
    
    withBusyIndicatorServer(btn_id[i], {
      withProgress(message = 'Calculating links number',
                   detail = 'please wait ...', value = 1,{
                     
                     if(btn_id[i] == "calcIntraction1"){
                       shinyjs::hide("IntConfig2")
                       edge_slc <- input$edgeSize2
                       updateSliderInput(session, "edgeSize3", value = input$edgeSize2, step = rv$binsize_ini , max = 500000)
                     } else {
                       shinyjs::hide("IntConfig1")
                       edge_slc <- input$edgeSize4
                       updateSliderInput(session, "edgeSize1", value = input$edgeSize4, step = rv$binsize_ini , max = 500000)
                     }
                     
                     rv$interseq_links <-  rv$sequence_length[rv$mat_binned[rname != mrnm,], on = c("rname")] %>%
                       .[rv$sequence_length, on = c("mrnm" = "rname"), ':='(mrnm_len = i.rname_len)] %>%
                       .[((pos < edge_slc | pos > rname_len - edge_slc) & (mpos < edge_slc | mpos > mrnm_len - edge_slc)),] %>%
                       .[, ':='(rname_strand = ifelse(pos < edge_slc, "-", 
                                                      ifelse(pos > rname_len - edge_slc, "+", "M")),
                                mrnm_strand = ifelse(mpos < edge_slc, "+", 
                                                     ifelse(mpos > mrnm_len - edge_slc, "-", "M")),
                                edge_rname = ifelse(rname_len < edge_slc, rname_len, edge_slc),
                                edge_mrnm = ifelse(mrnm_len < edge_slc, mrnm_len, edge_slc))] %>%
                       .[order(rname_len, decreasing = TRUE), .(link_no = .N, sum = sum(n), 
                                                                edge_rname = max(edge_rname), 
                                                                edge_mrnm = max(edge_mrnm)), 
                         by = .(rname, mrnm, rname_strand, mrnm_strand, rname_len, mrnm_len)] %>%
                       .[sum > 2,.(link_density = round(link_no/(ceiling(edge_rname/rv$binsize) * 
                                                                   ceiling(edge_mrnm/rv$binsize)),2), 
                                   link_no, sum, avg = sum/link_no), 
                         by = .(rname, mrnm, rname_strand, mrnm_strand, rname_len, mrnm_len)]
                     
                     #bin_rname =(link_no/ceiling(min(edge_rname,edge_mrnm)/rv$binsize))/10
                     choices <- unique(rv$interseq_links$rname)
                     choices_opt <- list(style = paste(rep_len("font-size: 12px; line-height: 1.5; margin-left: -10px; 
                                                               border-bottom: 1px solid gray;", length(choices)), 
                                                       "background-color: lightgray;"))
                     
                     updatePickerInput(session, "seq", 
                                       choices = c( "", choices), selected = isolate(input$seq),
                                       choicesOpt = choices_opt)
                     
                     updatePickerInput(session, "seq1", 
                                       choices = c( "", choices), selected = isolate(input$seq1),
                                       choicesOpt = choices_opt)
                     
                     updatePickerInput(session, "seq2", 
                                       choices = c( "", choices), selected = isolate(input$seq2),
                                       choicesOpt = choices_opt)
                     
                     updateNumericInput(session, "edgeSize3", min = rv$binsize, value = edge_slc, 
                                        step = rv$binsize , max = 500000)
                     
                     #set plot size to default
                     rv$intramap_range <- NULL
                     rv$edge_slc <- edge_slc
                   })
    })
  })
  
  ## or users can just import pre-calculated link density file
  observeEvent(input$importLnkFile, {
    withBusyIndicatorServer("importLnkFile", {
      withProgress(message = 'Importing link density file',
                   detail = 'please wait ...', value = 1, {
                     
                     interseq_links <- fread(file.path(rv$projDir,input$lnkFile)) 
                     
                     # These columns are being used for downstream calculations
                     required_cols <- c("rname", "mrnm", "mrnm_strand", "rname_strand", "mrnm_len", "rname_len", "link_density")
                     
                     if(any(required_cols %!in% names(interseq_links)))
                       stop(paste("Columns", paste0(required_cols, collapse = ", "), 
                                  "must be present in the file"))
                     
                     rv$interseq_links <- interseq_links
                   })
    })
  })
  
  output$interseqData <- DT::renderDataTable({
    DT::datatable(rv$interseq_links ,
                  escape = FALSE, filter = 'bottom', rownames= FALSE, 
                  class = 'nowrap display compact order-column cell-border stripe', 
                  extensions = c('Buttons', 'ColReorder'), selection = "single",
                  options = list(dom = 'lBfrtip',
                                 buttons = list(list(
                                   extend  = 'collection',
                                   buttons = c('excel', 'copy'),
                                   text    = 'Download shown entries'
                                 ), 'colvis'),
                                 colReorder = list(realtime = FALSE),
                                 lengthMenu = list(c(10, 15, 25, 50, 100, 1000), c('10','15', '25', '50', '100', '1000')),
                                 pageLength = 15,
                                 fixedColumns = list(leftColumns = 1),
                                 searchHighlight = TRUE,
                                 autoWidth = FALSE,
                                 scrollX = TRUE,
                                 initComplete = htmlwidgets::JS(
                                   "function(settings, json) {",
                                   "$(this.api().table().header()).css({'background-color': '#C8C8C8', 'color': '#000'});",
                                   "}"))
    )
  })
  
  
  output$interseqLinks <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), 'interseq_Links_EdgeSize',rv$edge_slc, '.csv', sep='')
    },
    content = function(con) {
      write.csv(rv$interseq_links , con)
    }
  )
  ##------------------------------------------------------------------------------------------------------------
  ##------------------------------------------------------------------------------------------------------------
  observe({
    shiny::validate(need(input$seq1, ""))
    shiny::validate(need(rv$interseq_links, ""))
    
    seq <- input$seq1
    
    if(input$dir1 == "Backward"){
      subseq1 <- rv$interseq_links[mrnm %in% seq & (!rname %in% seq),] %>%
        .[order(link_density, link_no, avg, decreasing = TRUE), .(x = rname, x_dir = rname_strand)] 
    } else {
      subseq1 <- rv$interseq_links[rname %in% seq & (!mrnm %in% seq),] %>%
        .[order(link_density, link_no, avg, decreasing = TRUE), .(x = mrnm, x_dir = mrnm_strand)]
    }
    
    choices <- unique(subseq1$x)
    updatePickerInput(session, "subseq1", choices = choices,
                      choicesOpt = list(style = paste(rep_len("font-size: 12px; line-height: 1.5; 
                                                              margin-left: -10px; border-bottom: 1px solid gray;", 
                                                              length(choices)), "background-color: lightgray; z-index: 9999999;")))
  })
  
  ##------------------------------------------------------------------------------------------------------------
  ##------------------------------------------------------------------------------------------------------------
  
  # click event to link density
  observe({
    shiny::validate(need(rv$mat_binned, ""))
    if(is.null(input$intramap_click)) return(NULL)
    
    xy <- function(e) {
      if(is.null(e)) return(NULL)
      x <- round(e$x/rv$binsize, 0)
      y <- round(e$y/rv$binsize, 0)
      
      if(abs(x -y) > 5) return(NULL)
      
      pos <- min(x, y) * rv$binsize
      return(pos)
    }
    
    rv$tgt_pos = xy(input$intramap_click)
  })
  
  observe({
    shiny::validate(need(rv$tgt_pos, ""))
    
    tgt_contig <- input$seq1
    edge_slc <- input$edgeSize1
    
    tgt_linkDen <- rv$mat_binned[rname == tgt_contig & mrnm == tgt_contig,] %>%
      .[(pos >= rv$tgt_pos - edge_slc & pos <= rv$tgt_pos & mpos > rv$tgt_pos &  mpos <= rv$tgt_pos + edge_slc),] %>%
      .[, .(link_no = .N, sum = sum(n), 
            min_pos = min(pos), min_mpos = min(mpos), 
            max_pos = max(pos), max_mpos = max(mpos)), 
        by = .(rname)] %>%
      .[,.(link_density = round(link_no/(ceiling((max_pos - min_pos)/rv$binsize) * ceiling((max_mpos - min_mpos)/rv$binsize)),2), 
           link_no, sum, avg = sum/link_no), 
        by = .(rname)]
    print(tgt_linkDen)
  })
  
  ##------------------------------------------------------------------------------------------------------------
  ##------------------------------------------------------------------------------------------------------------
  #output first map
  
  observeEvent(input$exitZoom1, {
    session$resetBrush("intramap_brush")
    rv$intramap_range <- NULL
  })
  
  observeEvent(input$clearBrush1, {
    session$resetBrush("intramap_brush")
    
    if(rv$cut_pos > 0){
      updateNumericInput(session, "cutPos", value = 0)
    }
  })
  
  observe({
    shiny::validate(need(input$seq, ""))
    input$cut1
    input$cut2
    
    p <- rv$mat_binned[rname == input$seq & mrnm == input$seq, ]
    
    if(nrow(p) == 0) 
      return(NULL)
    
    p <- p %>% ggplot(aes(x = pos, y = mpos, fill=log10(n/2))) +
      geom_tile(size = 0.2) +
      labs(title = input$seq1 , subtitle = paste0("Size: ", 
                                                  max(p$pos/1000000), " Mb"), x = "", y = "") +
      #scale_fill_gradient(low = "#ff9999", high = "red") +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      #coord_cartesian(ylim=c(-12,12)) +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(), 
            legend.position = "none",
            panel.border = element_rect(colour = "gray", fill = NA),
            panel.background = element_rect(fill = "white", colour = "white"))
    
    rv$cut_pos <- 0
    rv$intramap_plot <- p
  })
  
  output$intramap <- renderPlot({
    shiny::validate(need(rv$intramap_plot, ""))
    
    if(!is.null(rv$intramap_range)){
      p <- rv$intramap_plot +
        coord_cartesian(ylim = rv$intramap_range, 
                        xlim = rv$intramap_range) 
      
      session$resetBrush("intramap_brush")
    } else {
      p <- rv$intramap_plot 
    }
    
    if(rv$cut_pos > 0){
      p + geom_vline(xintercept = rv$cut_pos, color = "red", size = 0.3) +
        geom_hline(yintercept = rv$cut_pos, color = "red", size = 0.3)
    } else {
      p
    }
    
  }, height = function(){
    if(is.null(input$dimension[2])){
      400
    } else {
      0.798 * input$dimension[2]
    }
  }, width = function(){
    if(is.null(input$dimension[2])){
      400
    } else {
      0.798 * input$dimension[2]
    }
  })
  
  
  ###--------------------------------
  ##----------------join-------------
  ###--------------------------------
 output$joinedmap <- renderPlot({
    shiny::validate(need(input$action == "Join", ""))
    shiny::validate(need(input$seq1, ""))
    shiny::validate(need(input$subseq1, ""))
    shiny::validate(need(rv$interseq_links, ""))
    
    seq <- input$seq1
    subseq = paste0(strand, input$subseq1)
    strand <- ifelse(input$strand_1, "+", "-")
    
    join_maps(rv$mat_binned, seq = seq, subseq = subseq, 
                        direction = isolate(input$dir1), 
                        binsize = isolate(rv$binsize), 
                        output = "graph")
    

    
  }, height = function(){
    if(is.null(input$dimension[2])){
      400
    } else {
      0.798 * input$dimension[2]
    }
  }, width = function(){
    if(is.null(input$dimension[2])){
      400
    } else {
      0.798 * input$dimension[2]
    }
  })
  
  
  
  
  
  observeEvent(input$join,{
    shiny::validate(need(input$subseq1, ""))
    
    withBusyIndicatorServer("join", {
      withProgress(message = 'Joining the sequences',
                   detail = 'please wait ...', value = 1, {
                     
                     tgt_contig <- input$seq1
                     subseq_contig <- input$subseq1
                     
                     strand <- ifelse(input$strand_1, "+", "-")
                     new_contig <- paste0(tgt_contig, "|", strand,"|", subseq_contig)
                     
                     tgt_len <- max(rv$mat_binned[rname == tgt_contig, pos])
                     subseq_len <- max(rv$mat_binned[rname == subseq_contig, pos])
                     new_contig_len <- tgt_len + subseq_len
                     subseq_start <- tgt_len + rv$binsize
                     
                     rv$mat_binned <- rv$mat_binned %>%
                       #.[(rname %in% tgt_contig | mrnm %in% tgt_contig),] %>%
                       .[, ':='(pos = ifelse(rname == subseq_contig,  pos + subseq_start, pos ),
                                mpos = ifelse(mrnm == subseq_contig , mpos + subseq_start , mpos))] %>%
                       .[,':='(rname = ifelse(rname %in% c(tgt_contig, subseq_contig), new_contig, rname),
                               mrnm = ifelse(mrnm %in% c(tgt_contig, subseq_contig), new_contig, mrnm))]
                     
                     rv$sequence_length <- rv$sequence_length %>%
                       rbind(.,list(new_contig, new_contig_len))
                     
                     choices <- rv$mat_binned[,.(rname, rname_len = max(pos)), 
                                              by = .(rname)][order(-rname_len), rname]
                     choices_style <- list(style = paste(rep_len("font-size: 12px; line-height: 1.5; margin-left: -10px; 
                                                                 border-bottom: 1px solid gray;", length(choices)), 
                                                         "background-color: lightgray;"))
                     
                     updatePickerInput(session, "seq", choices = choices, selected = new_contig,
                                       choicesOpt = choices_style)
                     updatePickerInput(session, "seq1", choices = choices, selected = new_contig,
                                       choicesOpt = choices_style)
                     
                     updatePickerInput(session, "subseq1", selected = "")
                     
                     #set plot size to default
                     rv$intramap_range <- NULL
                   })
    })
  })
  

  ###--------save changes----------- 
  observeEvent(input$svChanges, {
    withBusyIndicatorServer("svChanges", {
      saveRDS(rv$mat, paste(Sys.Date(), "_", basename(input$mapFile), sep=""))
    })
  })
  
  ##------------------------------------------------------------------------------------------------------------
  ##------------------------------------------------------------------------------------------------------------
  observe({
    shiny::validate(need(input$seq2, ""))
    shiny::validate(need(rv$interseq_links, ""))
    shiny::validate(need(input$nrSeq, ""))
    
    len <- length(isolate(input$chained_seq))
    seq2 <- input$seq2
    leading_seq <- seq2[length(seq2)]
    leadingSeq_len <-  max(rv$mat_binned[rname == leading_seq, pos])
    maxLinkDen <- input$maxLinkDen
    minSeqLen <- input$minSeqLen * 1000000
    
    
    if(len > 1 & input$nrSeq > 1){
      if(leadingSeq_len < 20 * rv$binsize){
        leading_seq <- gsub("^[-+]", "",rv$s1)
      } else {
        leading_seq <- leading_seq
      }
    }
    
    if(input$dir2 == "Backward"){
      rv$subseq2 <- rv$interseq_links[mrnm %in% leading_seq & (!rname %in% leading_seq) & link_density <= maxLinkDen & rname_len >= minSeqLen,] %>%
        .[order(link_density, link_no, avg, decreasing = TRUE), .(x = rname, x_dir = rname_strand)] 
    } else {
      rv$subseq2 <- rv$interseq_links[rname %in% leading_seq & (!mrnm %in% leading_seq) & link_density <= maxLinkDen &  mrnm_len >= minSeqLen,] %>%
        .[order(link_density, link_no, avg, decreasing = TRUE), .(x = mrnm, x_dir = mrnm_strand)]
    }
    
  })
  
  observe({
    shiny::validate(need(rv$subseq2, ""))
    
    chr <- gsub("^[-+]", "", c(input$chained_seq, input$scaf_man))
    rv$subseq2.1 <- rv$subseq2[(!(x %in% chr)), ]
    rv$choices <- c(unique(rv$subseq2.1$x), chr)
    
    updatePickerInput(session, "subseq2", choices =  rv$choices,
                      choicesOpt = list(style = paste(rep_len("font-size: 12px; line-height: 1.5; 
                                                              margin-left: -10px; border-bottom: 1px solid gray;", 
                                                              length(rv$choices)), "background-color: lightgray;")))
    shinyWidgets::updateSwitchInput(session, "strand_3", value = rv$subseq2.1$x_dir[1] == "+")
  })
  
  
  # move up and down the option list ---------------
  observeEvent(input$down,{
    shiny::validate(need(rv$choices, ""))
    
    i <- which(rv$choices == input$subseq2)
    i = i + 1
    
    if(i > length(rv$choices)) return(NULL)
    s2 <- rv$choices[i]
    strand_3 <- rv$subseq2.1[x ==  s2, x_dir][1]
    
    updatePickerInput(session, "subseq2",  selected = rv$choices[i])
    shinyWidgets::updateSwitchInput(session, "strand_3", value = strand_3 == "+")
  })
  
  observeEvent(input$up, {
    i <- which(rv$choices == input$subseq2)
    i = i - 1
    
    if(i < 1) return(NULL)
    s2 <- rv$choices[i]
    strand_3 <- rv$subseq2.1[x ==  s2, x_dir][1]
    updatePickerInput(session, "subseq2",  selected = rv$choices[i])
    shinyWidgets::updateSwitchInput(session, "strand_3", value = strand_3 == "+")
  })
  #-------------------------------------------------
  
  observe({
    shiny::validate(need(input$subseq2, ""))
    shiny::validate(need(input$seq2, ""))
    shiny::validate(need(input$nrSeq, ""))
    
    withBusyIndicatorServer("add", {
      
      strand_2 <- ifelse(input$strand_2, "+", "-")
      strand_3 <- ifelse(input$strand_3, "+", "-")
      rv$s2.1 <- paste0(strand_3, input$subseq2)
      chr <- isolate(input$chained_seq)
      
      if(input$nrSeq == 1 | is.null(chr)){
        rv$s1 <- paste0(strand_2, input$seq2)
      } else {
        len <- length(chr)
        nrSeq <- min(len, input$nrSeq)
        
        if(len == 1){
          rv$s1 <- paste0(strand_2, input$seq2)
        } else {
          if(input$dir2 == "Forward"){
            rv$s1 <- chr[((len-nrSeq) + 1):len]
          } else {
            rv$s1 <- chr[1:nrSeq]
          }
        }
      }
      rv$sel_map <- join_maps_plus(rv$mat_binned, seq = rv$s1, subseq = rv$s2.1, binsize = isolate(rv$binsize),  
                                   direction = input$dir2, output = "data")
    })
  })
  
  
  observe({
    shiny::validate(need(rv$sel_map, ""))
    shiny::validate(need(input$seq1, ""))
    
    len <- rv$sel_map[, .(len = min(pos)), by = .(rname)] %>%
      .[order(len), ]
    
    p <- ggplot(rv$sel_map, aes(x = pos, y = mpos, fill=log10(n/2))) +
      geom_tile( size = 0.2) +
      #scale_fill_gradient(low = "#ff9999", high = "red") +
      labs(title = isolate(rv$s2.1) , subtitle = paste0("Total size: ", max(len$len/1000000), " Mb"), x = "", y = "") +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      geom_vline(xintercept = len$len, color = "gray", size = 0.3) +
      geom_hline(yintercept = len$len, color = "gray", size = 0.3) +
      #coord_cartesian(ylim=c(-12,12)) +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(), 
            legend.position = "none",
            panel.border = element_rect(colour = "gray", fill = NA),
            panel.background = element_rect(fill = "white", colour = "white"))
    
    if(!is.null(rv$map1_range)){
      p <- p + coord_cartesian(ylim = rv$map1_range, xlim = rv$map1_range) 
      session$resetBrush("map1_brush")
    }
    
    rv$map1_plot <- p 
  })
  
  
  output$map1 <- renderPlot({
    shiny::validate(need(rv$map1_plot, ""))
    rv$map1_plot 
    
    
  }, height = function(){
    if(is.null(input$dimension[2])){
      400
    } else {
      0.798 * input$dimension[2]
    }
  }, width = function(){
    if(is.null(input$dimension[2])){
      400
    } else {
      0.798 * input$dimension[2]
    }
  })
  
  observe({
    shiny::validate(need(rv$sel_map, ""))
    
    if(is.null(input$map1_brush)) return(NULL)
    
    xy_range_str <- function(e) {
      if(is.null(e)) return(NULL)
      xmin <- round(e$xmin, 1)
      range <- c(ifelse(xmin < 0, 0, round(xmin/rv$binsize, 0) * rv$binsize),
                 round(e$xmax/rv$binsize, 0) * rv$binsize)
      
      return(range)
    }
    rv$map1_range <- xy_range_str(input$map1_brush)
  })
  
  observeEvent(input$exitZoom2, {
    session$resetBrush("map1_brush")
    rv$map1_range <- NULL
  })
  
  observeEvent(input$clearBrush2, {
    session$resetBrush("map1_brush")
  })
  
  output$map1_hoverinfo <- renderTable({
    shiny::validate(need(rv$combined_maps, ""))
    # For base graphics, we need to specify columns, though for ggplot2,
    # it's usually not necessary.
    res <- nearPoints(rv$sel_map, input$map_hover, "pos", "mpos", threshold = 1, maxpoints = 1)
    if (nrow(res) == 0)
      return()
    res
  })
  
  
  ## Join sequences together ---------------------------
  observeEvent(input$add,{
    withBusyIndicatorServer("add", {
      
      if(length(input$chained_seq) == 0){
        if (input$dir2 == 'Backward'){
          rv$chr <- c( rv$s2.1, rv$s1)
        } else {
          rv$chr <- c(rv$s1,  rv$s2.1)
        }
      } else {
        if (input$dir2 == 'Backward'){
          rv$chr <- c(rv$s2.1, input$chained_seq)
        } else {
          rv$chr <-  c(input$chained_seq, rv$s2.1)
        }
      }
      updateCheckboxGroupInput(session, "chained_seq", NULL, choices = unique(rv$chr), selected = unique(rv$chr))
      updateTextAreaInput(session, "scaf_edit", NULL, value = paste(unique(rv$chr), collapse = "\n"))
    })
  })
  
  ## copy to clipboard
  observeEvent(input$clipbtn, clipr::write_clip(input$chained_seq))
  
  ## erase the list 
  observeEvent(input$erase,{ 
    choices =  input$chained_seq
    updateCheckboxGroupInput(session, "chained_seq", NULL, choices =  choices[1], selected =  choices[1])
  })
  
  ## edit the list 
  observeEvent(input$edit,{ 
    shinyjs::show("EditBox")
    hide("CheckBox")
    hide("edit")
    shinyjs::show("check")
  })
  
  observeEvent(input$check,{ 
    
    choices = base::strsplit(input$scaf_edit, "\n")[[1]]
    updateCheckboxGroupInput(session, "chained_seq", NULL, choices =  choices, selected =  choices)
    
    shinyjs::hide("EditBox")
    shinyjs::show("CheckBox")
    shinyjs::show("edit")
    shinyjs::hide("check")
  })
  
  
  observeEvent(input$export,{
    shiny::validate(need(input$chained_seq, ""))
    
    out_dir <- file.path(rv$projDir,"groups")
    dir.create(out_dir, showWarnings = FALSE)
    
    super_seq <- data.table(name = input$chained_seq) %>%
      .[, ':='(name = substr(name, 2, nchar(name)),
               rc=ifelse(grepl("\\+", name), 0, 1),
               q = ".", 
               gap_size = ".")]
    
    # get next group number
    n <- length(list.files(out_dir))
    
    file = paste0(out_dir, "/group",n,".ordering")
    write.table(super_seq, file, row.names = FALSE, col.names = FALSE, quote = FALSE) 
  })
  
  
  
  observe({
    chr <- input$chained_seq
    
    if(length(chr) > 1){
      disable("seq2"); disable_switch = TRUE;
    } else {
      enable("seq2"); disable_switch = FALSE;
    }
    
    if(length(chr) == 0){
      return(NULL)
    }
    
    leading_seq <- ifelse(input$dir2 == 'Backward', chr[1], chr[length(chr)])
    
    updatePickerInput(session, "seq2", selected = gsub("^[-+]","",leading_seq))
    shinyWidgets::updateSwitchInput(session, "strand_2", value = substr(leading_seq, 1, 1) == "+", disabled = disable_switch)
    updateNumericInput(session, "n", value = 1)
  })
  
  # plot hi-c map for chr
  #------------------------
  observeEvent(input$combineMaps,{
    withBusyIndicatorServer("combineMaps",{
      withProgress(message = 'Preparing contact data', value = 1, 
                   detail = "please be patient ...", {
                     #maps <- if(input$inputType == "Manual"){base::strsplit(input$scaf_man, "\n")[[1]]} else {input$chained_seq}
                     maps <- input$chained_seq
                     
                     rv$combined_maps <- join_maps_plus(mat = rv$mat_binned, seq = maps, direction = "Forward", 
                                                        output = "data", binsize = isolate(rv$binsize))
                     
                     #fwrite(chr, paste0("./","chr",chr_no,"_",n, ".csv"), row.names = FALSE)
                     #saveRDS(chr, paste0("chr1",chr_no,".rds"))
                     #chr <- fread(paste0("chr6-", 16,".csv"), stringsAsFactors = FALSE)
                   })
    })
  })
  
  output$map2 <- renderPlot({
    shiny::validate(need(rv$combined_maps, ""))
    
    len <- rv$combined_maps[, .(len = max(pos)), by = .(rname)] %>%
      .[order(len), ]
    
    ggplot(rv$combined_maps, aes(x = pos, y = mpos, fill=log10(n))) +
      geom_tile(size = 0.2) +
      labs(title = paste("Total size:", max(len$len)/1000000, "MB"), x = "", y = "") +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      #geom_vline(xintercept = len$len	, color = "gray", size = 0.2) +
      #geom_hline(yintercept = len$len	, color = "gray", size = 0.2) +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(), 
            legend.position = "none",
            panel.background = element_rect(fill = "white", colour = "white"))
  }, height = function() {
    if(is.null(input$dimension[2])){
      400
    } else {
      0.8 * input$dimension[2]
    }
  }, width = function() {
    if(is.null(input$dimension[2])){
      400
    } else {
      0.8 * input$dimension[2]
    }
  })
  
  output$map2_hoverinfo <- renderTable({
    shiny::validate(need(rv$combined_maps, ""))
    # For base graphics, we need to specify columns, though for ggplot2,
    # it's usually not necessary.
    res <- nearPoints(rv$combined_maps, input$map2_hover, "pos", "mpos", threshold = 1, maxpoints = 1)
    
    if(nrow(res) == 0)
      return()
    res
  })
}



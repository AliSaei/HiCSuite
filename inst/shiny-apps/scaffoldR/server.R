library(shiny)
library(data.table)
library(dplyr)
library(ggplot2)
library(shinyWidgets)
library(shinycustomloader)
library(shinyjs)
library(shinyFiles)
library(HiCSuite)

options(future.globals.maxSize = 50*1024^2,
        shiny.maxRequestSize=300*1024^2)
shinyOptions(progress.style="old")
output_dir = NULL

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

#sequence_length <- fread(paste0("./data/","bilberry/sequence_len_bilberry.csv"), col.names = c("rname","len"), select = c(1, 2)) 
#a <- anchored_contigs[sequence_length, on = c(name = "rname")]

#sum(a[!is.na(contig), len])/sum(a[, len])*100

server <- function(input, output, session) {
  rv <- reactiveValues(chr = NULL, s2.1 = NULL, cut_pos = 0, 
                       intramap_range = NULL,  btn_val2 = c(0,0))
  
  volumes <- c(Home = fs::path_home(), getVolumes()(), `Data Directory` = "../../../")
  
  shinyDirChoose(input, 'directory', roots = volumes, session = session,
                 restrictions = system.file(package = "base"))
  
  observe({
    rv$projDir <- parseDirPath(volumes, input$directory)
    
    mapFiles <- list.files(rv$projDir, pattern = ".rds$", recursive = TRUE)
    lenFiles <- list.files(rv$projDir, pattern = ".csv$|.txt$", recursive = TRUE)
    
    updatePickerInput(session, "mapFile", choices = mapFiles,                           
                      choicesOpt = list(style = paste(rep_len("font-size: 12px; line-height: 1.5; 
                                                              margin-left: -10px; border-bottom: 1px solid gray;", 
                                                              length(mapFiles)), "background-color: lightgray;"))
    )
    updatePickerInput(session, "lenFile", choices = lenFiles,                           
                      choicesOpt = list(style = paste(rep_len("font-size: 12px; line-height: 1.5; 
                                                              margin-left: -10px; border-bottom: 1px solid gray;", 
                                                              length(lenFiles)), "background-color: lightgray;"))
    )
  })
  
  observeEvent(input$import_data,{
    withBusyIndicatorServer("import_data", {
      withProgress(message = 'Loading in progress',
                   detail = 'Please wait...', value = 1,{
                     updatePickerInput(session, "seq1", selected = "")
                     
                     rv$mat <- as.data.table(readRDS(file.path(rv$projDir,input$mapFile))) 
                     
                     ## calculate or import length of sequences
                     if(input$loadSeqLen){
                       rv$sequence_length <- fread(file.path(rv$projDir,input$lenFile), col.names = c("rname","len"), select = c(1, 2)) 
                     } else {
                       rv$sequence_length <- rv$mat[,.(len = max(pos)), by = .(rname)]
                     }
                     
                     rmax = max(rv$mat$pos)
                     mmax = max(rv$mat$mpos)
                     rv$max = max(rmax, mmax)
                     binsize_ini <- sort.int(unique(rv$mat$pos), partial = 1:2)[1:2]
                     binsize_ini = as.integer(binsize_ini[2] - binsize_ini[1])
                     choices <- rv$sequence_length[rname %in% unique(rv$mat[, rname]), ][order(-len), rname]
                     choices_opt <- list(style = paste(rep_len("font-size: 12px; line-height: 1.5; margin-left: -10px; 
                                                        border-bottom: 1px solid gray;", length(choices)), 
                                                       "background-color: lightgray;"))
                     
                     updateNumericInput(session, "binsize2", value = binsize_ini, min = binsize_ini,  
                                        step = binsize_ini , max = 500000)
                     
                     updateNumericInput(session, "edgeSize1",  value = 20 * binsize_ini, min = binsize_ini, 
                                        step = binsize_ini , max = 500000)
                     
                     updateNumericInput(session, "edgeSize3",  value = 20 * binsize_ini, min = binsize_ini, 
                                        step = binsize_ini , max = 500000)
                     
                     updatePickerInput(session, "seq1", choices = choices, choicesOpt = choices_opt)
                     updatePickerInput(session, "seq2", choices = c( "", choices), selected = input$seq2, choicesOpt = choices_opt)
                     
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
  
  
  # intarmap view
  observe({
    shiny::validate(need(rv$mat_binned, ""))
    if(is.null(input$intramap_dblclick)) return(NULL)
    
    xy_str <- function(e) {
      if(is.null(e)) return(NULL)
      round(e$x/rv$binsize, 0) * rv$binsize
    }
    updateRadioGroupButtons(session, 'action', selected = "Cut")
    updateNumericInput(session, "cutPos", value = xy_str(input$intramap_dblclick))
  })
  
  
  observeEvent(input$cutPos,{
    rv$cut_pos <- input$cutPos
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
  
  
  observeEvent(input$cut1 + input$cut2,{
    shiny::validate(need(rv$cut_pos > 0, ""))
    
    withBusyIndicatorServer("cut1", {
      withProgress(message = 'Cutting the sequence',
                   detail = 'please wait ...', value = 1,{
                     
                     tgt_contig <- input$seq1
                     tgt_len <- rv$sequence_length[rname == tgt_contig, len]
                     subseq_start <- rv$cut_pos + rv$binsize
                     subseq_contig <- paste0(tgt_contig, "_subseq_", subseq_start, ":", tgt_len)
                     
                     rv$mat_binned <- rv$mat_binned %>%
                       #.[(rname %in% tgt_contig | mrnm %in% tgt_contig),] %>%
                       .[,':='(rname = ifelse(rname == tgt_contig & pos >  rv$cut_pos , subseq_contig, rname),
                               mrnm = ifelse(mrnm == tgt_contig & mpos >  rv$cut_pos , subseq_contig, mrnm))] %>%
                       .[, ':='(pos = ifelse(rname == subseq_contig,  pos - subseq_start, pos ),
                                mpos = ifelse(mrnm == subseq_contig , mpos - subseq_start , mpos))]
                     
                     rv$sequence_length <- rv$sequence_length[rname == tgt_contig, len := rv$cut_pos] %>%
                       rbind(.,list(subseq_contig, tgt_len - rv$cut_pos))
                     
                     updateNumericInput(session, "cutPos", value = 0)
                     
                     choices <- c(rv$seq_name, subseq_contig)
                     updatePickerInput(session, "seq1", choices = choices, selected = tgt_contig,
                                       choicesOpt = list(style = paste(rep_len("font-size: 12px; line-height: 1.5; 
                                                                               margin-left: -10px; border-bottom: 1px solid gray;", 
                                                                               length(choices)), "background-color: lightgray;")))
                     #set plot size to default
                     rv$intramap_range <- NULL
                   })
    })
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
                     
                     choices <- rv$mat_binned[,.(rname, len = max(pos)), 
                                              by = .(rname)][order(-len), rname]
                     choices_style <- list(style = paste(rep_len("font-size: 12px; line-height: 1.5; margin-left: -10px; 
                                                                 border-bottom: 1px solid gray;", length(choices)), 
                                                         "background-color: lightgray;"))
                     
                     updatePickerInput(session, "seq1", choices = choices, selected = new_contig,
                                       choicesOpt = choices_style)
                     
                     updatePickerInput(session, "subseq1", selected = "")
                     
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
                     
                     rv$interseq_contact <-  rv$sequence_length[rv$mat_binned[rname != mrnm,], on = c("rname")] %>%
                       .[rv$sequence_length, on = c("mrnm" = "rname"), ':='(mlen = i.len)] %>%
                       .[((pos < edge_slc | pos > len - edge_slc) & (mpos < edge_slc | mpos > mlen - edge_slc)),] %>%
                       .[, ':='(rname_dir = ifelse(pos < edge_slc, "-", 
                                                   ifelse(pos > len - edge_slc, "+", "M")),
                                mrnm_dir = ifelse(mpos < edge_slc, "+", 
                                                  ifelse(mpos > mlen - edge_slc, "-", "M")),
                                edge_rname = ifelse(len < edge_slc, len, edge_slc),
                                edge_mrnm = ifelse(mlen < edge_slc, mlen, edge_slc))] %>%
                       .[order(len, decreasing = TRUE), .(link_no = .N, sum = sum(n), 
                                                          edge_rname = max(edge_rname), 
                                                          edge_mrnm = max(edge_mrnm)), 
                         by = .(rname, mrnm, rname_dir, mrnm_dir, len, mlen)] %>%
                       .[sum > 2,.(link_density = round(link_no/(ceiling(edge_rname/rv$binsize) * 
                                                                   ceiling(edge_mrnm/rv$binsize)),2), 
                                   link_no, sum, avg = sum/link_no), 
                         by = .(rname, mrnm, rname_dir, mrnm_dir, len, mlen)]
                     
                     #bin_rname =(link_no/ceiling(min(edge_rname,edge_mrnm)/rv$binsize))/10
                     choices <- unique(rv$interseq_contact$rname)
                     choices_opt <- list(style = paste(rep_len("font-size: 12px; line-height: 1.5; margin-left: -10px; 
                                                               border-bottom: 1px solid gray;", length(choices)), 
                                                       "background-color: lightgray;"))
                     
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
                   })
    })
  })
  
  output$interseqData <- DT::renderDataTable({
    DT::datatable(rv$interseq_contact ,
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
  
  ##------------------------------------------------------------------------------------------------------------
  ##------------------------------------------------------------------------------------------------------------
  observe({
    shiny::validate(need(input$seq1, ""))
    shiny::validate(need(rv$interseq_contact, ""))
    
    seq <- input$seq1
    
    if(input$dir1 == "Backward"){
      subseq1 <- rv$interseq_contact[mrnm %in% seq & (!rname %in% seq),] %>%
        .[order(link_density, link_no, avg, decreasing = TRUE), .(x = rname, x_dir = rname_dir)] 
    } else {
      subseq1 <- rv$interseq_contact[rname %in% seq & (!mrnm %in% seq),] %>%
        .[order(link_density, link_no, avg, decreasing = TRUE), .(x = mrnm, x_dir = mrnm_dir)]
    }
    
    choices <- unique(subseq1$x)
    updatePickerInput(session, "subseq1", choices = choices,
                      choicesOpt = list(style = paste(rep_len("font-size: 12px; line-height: 1.5; 
                                                              margin-left: -10px; border-bottom: 1px solid gray;", 
                                                              length(choices)), "background-color: lightgray; z-index: 9999999;")))
    #shinyWidgets::updateSwitchInput(session, "strand_3.1", value = subseq1$x_dir[1] == "+")
  })
  
  ##------------------------------------------------------------------------------------------------------------
  ##------------------------------------------------------------------------------------------------------------
  observeEvent(input$svChanges, {
    withBusyIndicatorServer("svChanges", {
      saveRDS(rv$mat, paste(Sys.Date(), "_", basename(input$mapFile), sep=""))
    })
  })
  
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
    shiny::validate(need(input$seq1, ""))
    input$cut1
    input$cut2
    
    if(isolate(input$action == "Join")){
      if(is.null(input$subseq1)) 
        return(NULL)
      
      seq <- input$seq1
      rv$seq_len <- max(rv$mat_binned[rname == seq, pos])
      
      strand <- ifelse(input$strand_1, "+", "-")
      subseq = paste0(strand,input$subseq1)
      
      p <- JoinMultiMaps(rv$mat_binned, seq = seq, subseq = subseq, 
                         direction = isolate(input$dir1), 
                         binsize = isolate(rv$binsize), 
                         output_dir = output_dir, 
                         output = "data")
    } else {
      p <- rv$mat_binned[rname == input$seq1 & mrnm == input$seq1, ]
      
      if(nrow(p) == 0) 
        return(NULL)
    }
    
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
    
    
    rv$intramap_plot <- p
    rv$tgt_pos <- NULL
  })
  
  output$intramap <- renderPlot({
    shiny::validate(need(rv$intramap_plot, ""))
    
    if(!is.null(rv$intramap_range)){
      p <- rv$intramap_plot +
        coord_cartesian(ylim = rv$intramap_range, 
                        xlim = rv$intramap_range) +
        geom_vline(xintercept = c(rv$seq_len, rv$cut_pos)	, color = "red", size = 0.3) +
        geom_hline(yintercept = c(rv$seq_len, rv$cut_pos), color = "red", size = 0.3) 
      
      session$resetBrush("intramap_brush")
    } else {
      p <- rv$intramap_plot +
        geom_vline(xintercept = c(rv$seq_len, rv$cut_pos)	, color = "red", size = 0.3) +
        geom_hline(yintercept = c(rv$seq_len, rv$cut_pos), color = "red", size = 0.3) 
    }
    
    p
    
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
  
  
  
  
  ##------------------------------------------------------------------------------------------------------------
  ##------------------------------------------------------------------------------------------------------------
  observe({
    shiny::validate(need(input$seq2, ""))
    shiny::validate(need(rv$interseq_contact, ""))
    
    len <- length(isolate(input$scaf_auto))
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
      rv$subseq2 <- rv$interseq_contact[mrnm %in% leading_seq & (!rname %in% leading_seq) & link_density < maxLinkDen & len > minSeqLen,] %>%
        .[order(link_density, link_no, avg, decreasing = TRUE), .(x = rname, x_dir = rname_dir)] 
    } else {
      rv$subseq2 <- rv$interseq_contact[rname %in% leading_seq & (!mrnm %in% leading_seq) & link_density < maxLinkDen &  mlen > minSeqLen,] %>%
        .[order(link_density, link_no, avg, decreasing = TRUE), .(x = mrnm, x_dir = mrnm_dir)]
    }
    
  })
  
  observe({
    shiny::validate(need(rv$subseq2, ""))
    
    chr <- gsub("^[-+]", "",c(input$scaf_auto, input$scaf_man))
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
      chr <- isolate(input$scaf_auto)
      
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
      rv$sel_map <- JoinMultiMaps(rv$mat_binned, seq = rv$s1, subseq = rv$s2.1, binsize = isolate(rv$binsize),  
                                  direction = input$dir2, output_dir = output_dir, output = "data")
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
  
  ## Batch run ------------------------------------
  observeEvent(input$build_map,{
    shiny::validate(need(rv$subseq2, ""))
    
    cmd <- paste0("open ", normalizePath(output_dir))
    system(cmd)
    
    direction <- input$dir2
    seq <- input$seq2
    subseq <- rv$subseq2$mrnm[1:20]
    
    #future({
    # draw the interaction map of the intial sequence with the other sequence individually and store plots on disk
    join2maps(mat, seq, subseq, direction = direction, output_dir = output_dir)
    #})
  })
  
  ## Join sequences together ---------------------------
  observeEvent(input$add,{
    withBusyIndicatorServer("add", {
      
      if(length(input$scaf_auto) == 0){
        if (input$dir2 == 'Backward'){
          rv$chr <- c( rv$s2.1, rv$s1)
        } else {
          rv$chr <- c(rv$s1,  rv$s2.1)
        }
      } else {
        if (input$dir2 == 'Backward'){
          rv$chr <- c(rv$s2.1, input$scaf_auto)
        } else {
          rv$chr <-  c(input$scaf_auto, rv$s2.1)
        }
      }
      updateCheckboxGroupInput(session, "scaf_auto", NULL, choices = unique(rv$chr), selected = unique(rv$chr))
      updateTextAreaInput(session, "scaf_edit", NULL, value = paste(unique(rv$chr), collapse = "\n"))
    })
  })
  
  ## copy to clipboard
  observeEvent(input$clipbtn, clipr::write_clip(input$scaf_auto))
  
  ## erase the list 
  observeEvent(input$erase,{ 
    choices =  input$scaf_auto
    updateCheckboxGroupInput(session, "scaf_auto", NULL, choices =  choices[1], selected =  choices[1])
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
    updateCheckboxGroupInput(session, "scaf_auto", NULL, choices =  choices, selected =  choices)
    
    shinyjs::hide("EditBox")
    shinyjs::show("CheckBox")
    shinyjs::show("edit")
    shinyjs::hide("check")
  })
  
  observe({
    chr <- input$scaf_auto
    
    if(length(chr) > 1){
      disable("seq2")
      disable_switch = TRUE
    } else {
      enable("seq2")
      disable_switch = FALSE
    }
    
    if(length(chr) == 0) return(NULL)
    
    leading_seq <- ifelse(input$dir2 == 'Backward', chr[1], chr[length(chr)])
    #leading_seq_reverse <- ifelse(grepl("^-", leading_seq), TRUE, FALSE)
    #strand_2 <- substr(leading_seq, 1, 1)
    
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
                     #maps <- if(input$inputType == "Manual"){base::strsplit(input$scaf_man, "\n")[[1]]} else {input$scaf_auto}
                     maps <- input$scaf_auto
                     
                     rv$combined_maps <- JoinMultiMaps(mat = rv$mat_binned, seq = maps, direction = "Forward", 
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



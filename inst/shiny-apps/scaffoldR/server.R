
source("./globals/functions.R")
source("./globals/helpers.R")

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
    
    rds_list <- list.files(rv$projDir, pattern = ".rds$|.bin|.txt|.csv", 
                           recursive = TRUE, ignore.case = TRUE)
    txt_list <- list.files(rv$projDir, pattern = ".csv$|.txt$|.tab$", 
                           recursive = TRUE, ignore.case = TRUE)
    fa_list <- list.files(rv$projDir, pattern = ".fa$|.fasta$", 
                           recursive = TRUE, ignore.case = TRUE)
    
    
    updatePickerInput(session, "mapFile", choices = rds_list,                           
                      choicesOpt = list(style = paste(
                        rep_len("font-size: 12px; line-height: 1.5; 
                      margin-left: -10px; border-bottom: 1px solid gray;", 
                                length(rds_list)), "background-color: lightgray;")))
    updatePickerInput(session, "lenFile", choices = txt_list,                           
                      choicesOpt = list(style = paste(
                        rep_len("font-size: 12px; line-height: 1.5; 
                      margin-left: -10px; border-bottom: 1px solid gray;", 
                                length(txt_list)), "background-color: lightgray;")))
    
    updatePickerInput(session, "faFile", choices = fa_list,                           
                      choicesOpt = list(style = paste(
                        rep_len("font-size: 12px; line-height: 1.5; 
                      margin-left: -10px; border-bottom: 1px solid gray;", 
                                length(fa_list)), "background-color: lightgray;")))
    
    updatePickerInput(session, "lnkFile", choices = rds_list,                           
                      choicesOpt = list(style = paste(
                        rep_len("font-size: 12px; line-height: 1.5; 
                      margin-left: -10px; border-bottom: 1px solid gray;", 
                                length(txt_list)), "background-color: lightgray;")))
  })
  
  
  ## Import contact map genearted using the jupyternote book script from BAM file
  observeEvent(input$import_data,{
    withBusyIndicatorServer("import_data", {
      withProgress(message = 'Loading in progress',
                   detail = 'Please wait...', value = 1,{
                     #updatePickerInput(session, "seq1", selected = "")
                     
                     if(grepl(".rds",input$mapFile, ignore.case = TRUE)){
                       rv$contact_data <- as.data.table(as.data.frame(
                         readRDS(file.path(rv$projDir, input$mapFile))
                       ))
                     } else {
                       rv$contact_data <- fread(file.path(rv$projDir, input$mapFile), 
                                                stringsAsFactors = FALSE)
                     }
                     
                     #names(rv$contact_data) <- c("rname","pos","mrnm","mpos","n")
                     
                     binsize_ini <- sort.int(unique(rv$contact_data$pos), partial = 1:2)[1:2]
                     binsize_ini = as.integer(binsize_ini[2] - binsize_ini[1])
                     
                     ## calculate or import length of sequences
                     if(input$loadSeqLen){
                       rv$seq_len <- fread(file.path(rv$projDir,input$lenFile), 
                                           col.names = c("rname","rlen"), 
                                           select = c(1, 2)) %>%
                         .[order(-rlen),]
                     } else {
                       rv$seq_len <- rv$contact_data[,.(rlen = max(pos)), by = .(rname)][order(-rlen),]
                     }
                     
                     
                     updateNumericInput(session, "binsize2", value = binsize_ini, 
                                        min = binsize_ini,  
                                        step = binsize_ini , max = 500000)
                     
                     updateNumericInput(session, "edgeSize1",  value = 20 * binsize_ini, 
                                        min = binsize_ini, 
                                        step = binsize_ini , max = 500000)
                     
                     updateNumericInput(session, "edgeSize3",  value = 20 * binsize_ini, 
                                        min = binsize_ini, 
                                        step = binsize_ini , max = 500000)
                     ##---------------------------------------------------------
                     #[rname %in% unique(rv$contact_data[, rname]), ]
                     rv$choices <- rv$seq_len[, rname]
                     choices_opt <- list(style = paste(rep_len("font-size: 12px; line-height: 1.5; margin-left: -10px; 
                                                        border-bottom: 1px solid gray;", length(rv$choices)), 
                                                       "background-color: lightgray;"))
                     
                     updatePickerInput(session, "seq", choices = rv$choices, 
                                       choicesOpt = choices_opt)
                     updatePickerInput(session, "seq1", choices = rv$choices, 
                                       choicesOpt = choices_opt)
                     updatePickerInput(session, "seq2", choices = c( "", rv$choices), 
                                       selected = input$seq2, 
                                       choicesOpt = choices_opt)
                     
                     rv$binsize_ini <- binsize_ini
                     rv$binsize <- binsize_ini
                     rv$contact_data2 <- rv$contact_data
                   })
    })
  })
  
  ## Import assembled scaffolds 
  observeEvent(input$import_fasta,{
    withBusyIndicatorServer("import_fasta", {
      withProgress(message = 'Loading in progress',
                   detail = 'Please wait...', value = 1,{
                     
                     rv$fasta <-  readDNAStringSet(file.path(rv$projDir, input$faFile))
                   })
    })
  })
  

  
  ##----------------------------------------------------------------------------
  ##----------------------------------------------------------------------------
  # Bin Size
  observe({
    updateSliderInput(session, "binsize", value = input$binsize2, 
                      step = rv$binsize_ini , max = input$binsize2)
  })
  
  observeEvent(input$update_bin, {
    withBusyIndicatorServer("update_bin", {
      withProgress(message = 'Binning in progress',
                   detail = 'please wait ...', value = 1, {
                     rv$binsize <- input$binsize
                     
                     rmax <- max(rv$contact_data$pos)
                     mmax <- max(rv$contact_data$mpos)
                     rv$max <- max(rmax, mmax)
                     
                     if(rv$binsize_ini != rv$binsize){
                       bins =  seq(rv$binsize , rv$max, rv$binsize)
                       breaks =  c(0, bins, (length(bins) + 1) * rv$binsize)
                       
                       rv$contact_data2 <- rv$contact_data %>%
                         .[, ':='(rbin = cut(pos,  breaks, include.lowest = TRUE),
                                  mbin = cut(mpos,  breaks, include.lowest = TRUE))] %>%
                         .[, .(n = sum(n)), by = .(rname, mrnm, rbin, mbin)] %>%
                         .[, pos :=  as.numeric(substr(rbin, 2, gregexpr(",", rbin)[[1]][1]-1)), by = "rbin"] %>%
                         .[,mpos :=  as.numeric(substr(mbin, 2, gregexpr(",", mbin)[[1]][1]-1)), by = "mbin"] %>%
                         .[, .(rname,pos, mrnm, mpos, n)]
                     } else {
                       rv$contact_data2 <- rv$contact_data
                     }
                   })
    })
  })
  
  output$seqLen <- DT::renderDT({
    shiny::validate(need(rv$seq_len, ""))
    
    datatable(rv$seq_len,
              rownames = FALSE, class = 'display compact row-border', 
              selection = 'single', filter = 'bottom', colnames = c("Sequence", "Length"), 
              options = list(
                pageLength = 15, dom = 'lti', autoWidth = TRUE,
                lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
                scrollY = "400px",
                initComplete = JS(
                  "function(settings, json) {",
                  "$(this.api().table().header()).css({'background-color': 'lightgray', 'color': '#000'});",
                  "}")
              )
    )
  })
  
  # output$lenHist <- renderPlot({
  #   shiny::validate(need(rv$seq_len, ""))
  #   
  #   ggplot(rv$seq_len, aes(x = rlen)) + 
  #     geom_histogram(breaks=seq(1, max(rv$seq_len$rlen), by = 10000),
  #                    fill = "blue", alpha = 0.2)
  # })
  
  output$contactData <- DT::renderDT({
    shiny::validate(need(rv$seq_len, ""))
    
    seq_index <- input$seqLen_rows_selected
    
    if(is.null(seq_index)) return(NULL)
    
    seq_selected <- rv$seq_len$rname[seq_index]

    datatable(
      rv$contact_data2[rname == seq_selected | mrnm == seq_selected,],
      rownames = FALSE, class = 'display compact cell-border', filter = 'bottom',
      options = list(
        pageLength = 15, dom = 'lti', autoWidth = TRUE,
        initComplete = JS(
          "function(settings, json) {",
          "$(this.api().table().header()).css({'background-color': '#F0F0F0', 'color': '#000'});",
          "}")
      )
    )
  })
  
  
  
  ## sequence length histogram
  output$lenHist <- renderPlot({
    shiny::validate(need(rv$seq_len, ""))
    
    rv$seq_len %>%
      ggplot(aes(x = rlen))+
      geom_histogram(bins = 50, 
                     color = "black", 
                     fill = "dodgerblue",
                     alpha = 0.5) +
      geom_vline(aes(xintercept = mean(rlen)),
                 color="red",  size=1) 
  })
  
  ##------------------------------------------------------------------------------------------------------------
  ##------------------------------------------------------------------------------------------------------------  
  # capture click events on intra-map view
  observe({
    shiny::validate(need(rv$contact_data2, ""))
    if(is.null(input$intramap_dblclick)) 
      return(NULL)
    
    xy_str <- function(e) {
      if(is.null(e)) return(NULL)
      round(e$x/rv$binsize, 0) * rv$binsize
    }
    
    updateRadioGroupButtons(session, 'action', selected = "Cut")
    updateNumericInput(session, "cutPos", value = xy_str(input$intramap_dblclick))
  })
  
  ##----------------------------------------------------------------------------
  # capture brush event
  observe({
    shiny::validate(need(rv$binsize > 0, ""))
    
    if(is.null(input$intramap_brush)) 
      return(NULL)
    
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
  
  ##----------------------------------------------------------------------------
  ##------------------- cut sequence--------------------------------------------
  observeEvent(input$cut1 + input$cut2,{
    shiny::validate(need(input$cutPos > 0, ""))
    
    withBusyIndicatorServer("cut1", {
      withProgress(message = 'Cutting the sequence',
                   detail = 'please wait ...', value = 1,{
                     
                     tgt_contig <- input$seq
                     tgt_len <- rv$seq_len[rname == tgt_contig, rlen]
                     frag_start <- rv$cut_pos + rv$binsize
                     
                     if(input$frag2Name == "Default"){
                       frag2_name <- paste0(tgt_contig, "_fragment_", frag_start, ":", tgt_len)
                     } else {
                       frag2_name <- input$frag2Name
                     }
                     
                     
                     rv$tgt_contig_data <- rv$contact_data2 %>%
                       .[(rname == tgt_contig | mrnm == tgt_contig),] %>%
                       .[,':='(rname = ifelse(rname == tgt_contig & pos >  rv$cut_pos , frag2_name, rname),
                               mrnm = ifelse(mrnm == tgt_contig & mpos >  rv$cut_pos , frag2_name, mrnm))] %>%
                       .[, ':='(pos = ifelse(rname == frag2_name,  pos - frag_start, pos ),
                                mpos = ifelse(mrnm == frag2_name , mpos - frag_start , mpos))]
                     
                     rv$contact_data2 <- rv$contact_data2 %>%
                       .[(rname != tgt_contig & mrnm != tgt_contig),] %>%
                       rbind(., rv$tgt_contig_data)
                     
                     ## Get new sequence length
                     rv$seq_len <- rv$contact_data2[,.(rlen = max(pos)), by = .(rname)] %>%
                       .[order(rname),]
                     
                     #rv$seq_len <- rv$seq_len[rname == tgt_contig, rlen := rv$cut_pos] %>%
                     #  rbind(.,list(frag2_name, tgt_len - rv$cut_pos)) %>%
                     #  .[order(-rlen),]
                     
                     rv$choices <- rv$seq_len[["rname"]]
                     updatePickerInput(session, "seq", choices = rv$choices, selected = tgt_contig,
                                       choicesOpt = list(style = paste(
                                         rep_len("font-size: 12px; line-height: 1.5; 
                                       margin-left: -10px; border-bottom: 1px solid gray;", 
                                                 length(rv$choices)), "background-color: lightgray;")))
                     ## Set plot size to default
                     rv$intramap_range <- NULL
                     ## Set break position to 0
                     updateNumericInput(session, "cutPos", value = 0)
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
  
  ##----------------------------------------------------------------------------
  ## Calculate inter-sequence interactions
  observeEvent(input$calcIntraction1 + input$calcIntraction2, {
    shiny::validate(need(rv$contact_data2, ""))
    
    btn_id <- c("calcIntraction1" , "calcIntraction2")
    btn_val <- c(input$calcIntraction1, input$calcIntraction2) - rv$btn_val2 
    i <- which(btn_val == max(btn_val))
    rv$btn_val2 <- btn_val
    
    
    withBusyIndicatorServer(btn_id[i], {
      withProgress(message = 'Calculating links number',
                   detail = 'please wait ...', value = 1, {
                     
                     if(btn_id[i] == "calcIntraction1"){
                       shinyjs::hide("IntConfig2")
                       edge_slc <- input$edgeSize2
                       updateSliderInput(session, "edgeSize3", 
                                         value = input$edgeSize2, 
                                         step = rv$binsize_ini , max = 500000)
                     } else {
                       shinyjs::hide("IntConfig1")
                       edge_slc <- input$edgeSize4
                       updateSliderInput(session, "edgeSize1", 
                                         value = input$edgeSize4, 
                                         step = rv$binsize_ini , max = 500000)
                     }
                     
                     rv$interseq_link_counts <-  rv$seq_len[rv$contact_data2[rname != mrnm,], on = c("rname")] %>%
                       .[rv$seq_len, on = c("mrnm" = "rname"), ':='(mrnm_len = i.rlen)] %>%
                       .[((pos < edge_slc | pos > rlen - edge_slc) & (mpos < edge_slc | mpos > mrnm_len - edge_slc)),] %>%
                       .[, ':='(rname_strand = ifelse(pos < edge_slc, "-", 
                                                      ifelse(pos > rlen - edge_slc, "+", "M")),
                                mrnm_strand = ifelse(mpos < edge_slc, "+", 
                                                     ifelse(mpos > mrnm_len - edge_slc, "-", "M")),
                                edge_rname = ifelse(rlen < edge_slc, rlen, edge_slc),
                                edge_mrnm = ifelse(mrnm_len < edge_slc, mrnm_len, edge_slc))] %>%
                       .[order(rlen, decreasing = TRUE), .(link_no = .N, sum = sum(n), 
                                                           edge_rname = max(edge_rname), 
                                                           edge_mrnm = max(edge_mrnm)), 
                         by = .(rname, mrnm, rname_strand, mrnm_strand, rlen, mrnm_len)] %>%
                       .[,.(link_density = round(link_no/(ceiling(edge_rname/rv$binsize) * 
                                                            ceiling(edge_mrnm/rv$binsize)),2), 
                            link_no, sum, avg = sum/link_no), 
                         by = .(rname, mrnm, rname_strand, mrnm_strand, rlen, mrnm_len)]
                     
                     #bin_rname =(link_no/ceiling(min(edge_rname,edge_mrnm)/rv$binsize))/10
                     rv$choices <- unique(rv$interseq_link_counts$rname)
                     choices_opt <- list(style = paste(rep_len("font-size: 12px; line-height: 1.5; margin-left: -10px; 
                                                               border-bottom: 1px solid gray;", length(rv$choices)), 
                                                       "background-color: lightgray;"))
                     
                     updatePickerInput(session, "seq", 
                                       choices = c( "", rv$choices), selected = isolate(input$seq),
                                       choicesOpt = choices_opt)
                     
                     updatePickerInput(session, "seq1", 
                                       choices = c( "", rv$choices), selected = isolate(input$seq1),
                                       choicesOpt = choices_opt)
                     
                     updatePickerInput(session, "seq2", 
                                       choices = c( "", rv$choices), selected = isolate(input$seq2),
                                       choicesOpt = choices_opt)
                     
                     updateNumericInput(session, "edgeSize3", min = rv$binsize, value = edge_slc, 
                                        step = rv$binsize , max = 500000)
                     
                     #set plot size to default
                     rv$edge_slc <- edge_slc
                     shinyjs::hide("IntConfig1", anim = TRUE)
                     shinyjs::hide("IntConfig2", anim = TRUE)
                   })
    })
  })
  
  ##----------------------------------------------------------------------------
  ## Or import pre-calculated interaction file
  observeEvent(input$importLnkFile, {
    withBusyIndicatorServer("importLnkFile", {
      withProgress(message = 'Importing interactions file',
                   detail = 'please wait ...', value = 1, {
                     
                     if(grepl(".rds",input$lnkFile, ignore.case = TRUE)){
                       interseq_link_counts <- readRDS(file.path(rv$projDir, input$lnkFile))
                     } else {
                       interseq_link_counts <- fread(file.path(rv$projDir,input$lnkFile))
                     }
                     
              
                     # These columns are being used for downstream calculations
                     required_cols <- c("rname", "mrnm", "mrnm_strand", "rname_strand", "mrnm_len", "rlen", "link_density")
                     
                     if(any(required_cols %!in% names(interseq_link_counts)))
                       stop(paste("Columns", paste0(required_cols, collapse = ", "), 
                                  "must be present in the file"))
                     
                     rv$interseq_link_counts <- interseq_link_counts
                   })
    })
  })
  
  output$interactionCounts <- DT::renderDataTable({
    DT::datatable(rv$interseq_link_counts[link_no > 10,],
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
  
  
  output$interactionCountsExp <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), '-interseq_link_counts-',rv$edge_slc, '.csv', sep='')
    },
    content = function(con) {
      fwrite(rv$interseq_link_counts , con, row.names = FALSE)
    }
  )
  
  observeEvent(input$svInteractionCounts, {
    withBusyIndicatorServer("svInteractionCounts", {
      file_name <- paste(Sys.Date(), '-interseq_link_counts-',rv$edge_slc, '.rds', sep='')
      
      saveRDS(rv$contact_data2, 
              file.path(rv$projDir, file_name)
      )
    })
  })
  ##------------------------------------------------------------------------------------------------------------
  ##------------------------------------------------------------------------------------------------------------
  observe({
    shiny::validate(need(input$seq1, ""))
    shiny::validate(need(rv$interseq_link_counts, ""))
    
    seq <- input$seq1
    
    if(input$dir1 == "Backward"){
      subseq1 <- rv$interseq_link_counts[mrnm %in% seq & (!rname %in% seq),] %>%
        .[order(link_density, link_no, avg, decreasing = TRUE), .(Subsequent_seq = rname, Strand = rname_strand)] 
    } else {
      subseq1 <- rv$interseq_link_counts[rname %in% seq & (!mrnm %in% seq),] %>%
        .[order(link_density, link_no, avg, decreasing = TRUE), .(Subsequent_seq = mrnm, Strand = mrnm_strand)]
    }
    
    choices <- unique(subseq1$Subsequent_seq)
    updatePickerInput(session, "subseq1", choices = choices,
                      choicesOpt = list(style = paste(
                        rep_len("font-size: 12px; line-height: 1.5; 
                                 margin-left: -10px; border-bottom: 1px solid gray;", 
                                length(choices)), 
                        "background-color: lightgray; z-index: 9999999;")))
  })
  
  ##------------------------------------------------------------------------------------------------------------
  ##------------------------------------------------------------------------------------------------------------
  
  # click event to link density
  observe({
    shiny::validate(need(rv$contact_data2, ""))
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
    
    tgt_linkDen <- rv$contact_data2[rname == tgt_contig & mrnm == tgt_contig,] %>%
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
  
  # move up and down the options list ---------------
  observeEvent(input$down,{
    shiny::validate(need(rv$choices, ""))
    
    i <- which(rv$choices == input$seq)
    i = i + 1
    
    if(i > length(rv$choices)){
      return(NULL)
    } else {
      updatePickerInput(session, "seq",  
                        selected = rv$choices[i])
    }
    
  })
  
  observeEvent(input$up, {
    shiny::validate(need(rv$choices, ""))
    
    i <- which(rv$choices == input$seq)
    i = i - 1
    
    if(i < 1){
      return(NULL)
    } else {
      updatePickerInput(session, "seq",  
                        selected = rv$choices[i])
    }
  })
  
  ##-------------------------------------------------------
  ##-------------------------------------------------------
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
    
    withProgress(message = 'Preparing plot',
                 detail = 'please wait ...', value = 1,{
                   p <- rv$contact_data2[rname == input$seq & mrnm == input$seq, ]
                   
                   if(nrow(p) == 0) 
                     return(NULL)
                   
                   p <- p %>% ggplot(aes(x = pos, y = mpos, fill=log10(n/2))) +
                     geom_tile(color = "red", size = 0.1) +
                     scale_fill_gradient(low = "#ffe6e6", high = "red") +
                     scale_x_continuous(expand=c(0,0)) +
                     scale_y_continuous(expand=c(0,0)) +
                     labs(title = input$seq , subtitle = paste0("Size: ", 
                                                                max(p$pos/1000000), " Mb"), x = "", y = "") +
                     #coord_cartesian(ylim=c(-12,12)) +
                     theme(axis.text = element_blank(),
                           axis.ticks = element_blank(), 
                           legend.position = "none",
                           panel.border = element_rect(colour = "gray", fill = NA),
                           panel.background = element_rect(fill = "white", colour = "white"))
                   
                   #updateNumericInput(session, "cutPos", value = 0)
                   # rv$cut_pos <- 0
                   rv$intramap_range <- NULL
                   rv$intramap_plot <- p
                 })
    
  })
  
  output$intramap <- renderPlot({
    shiny::validate(need(rv$intramap_plot, ""))
    
    if(!is.null(rv$intramap_range) && diff(rv$intramap_range) != 0){
      p <- rv$intramap_plot +
        coord_cartesian(ylim = rv$intramap_range, 
                        xlim = rv$intramap_range) 
      
      session$resetBrush("intramap_brush")
    } else {
      p <- rv$intramap_plot 
    }
    
    if(rv$cut_pos > 0){
      p + geom_vline(xintercept = rv$cut_pos, color = "blue", size = 0.3) +
        geom_hline(yintercept = rv$cut_pos, color = "blue", size = 0.3)
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
  
  joinedmap_data <- reactive({
    shiny::validate(need(input$action == "Join", ""))
    shiny::validate(need(input$seq1, ""))
    shiny::validate(need(input$subseq1, ""))
    shiny::validate(need(rv$interseq_link_counts, ""))
    
    strand <- ifelse(input$subseq1_strand, "+", "-")
    subseq = paste0(strand, input$subseq1)
    
    join_maps(rv$contact_data2, seq = input$seq1, subseq = subseq, 
              direction = input$dir1, 
              binsize = isolate(rv$binsize))
  })
  
  
  output$joinedmap <- renderPlot({
    shiny::validate(need(joinedmap_data(), ""))
    
    len <- joinedmap_data()[, .(len= max(pos)), 
                            by = .(rname)] %>%
      .[len == min(len), ]
    
    ggplot(joinedmap_data(),aes(x = pos, y = mpos, fill=log10(n))) +
      geom_tile(color = "red", size = 0.1) +
      scale_fill_gradient(low = "#ffe6e6", high = "red") +
      labs(title = paste("Total size:", max(len$len)/1000000, "MB"), 
           subtitle = len$rname[1], x = "", y = "") +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      geom_vline(xintercept = len$len	, color = "blue", size = 0.1) +
      geom_hline(yintercept = len$len	, color = "blue", size = 0.1) +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(), 
            legend.position = "none",
            panel.border = element_rect(colour = "gray", fill = NA),
            panel.background = element_rect(fill = "white", colour = "white"))
    
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
                     
                     strand <- ifelse(input$subseq1_strand, "+", "-")
                     
                     if(input$dir1 == "Forward"){
                       new_scaffold <- paste0("+", tgt_contig," & ",strand, subseq_contig)
                     } else {
                       new_scaffold <- paste0(strand, subseq_contig, " & +",tgt_contig)
                     }
                     
                     tgt_len <- max(rv$contact_data2[rname == tgt_contig, pos])
                     subseq_len <- max(rv$contact_data2[rname == subseq_contig, pos])
                     new_contig_len <- tgt_len + subseq_len
                     subseq_start <- tgt_len + rv$binsize
                     
                     rv$contact_data2 <- rv$contact_data2 %>%
                       #.[(rname %in% tgt_contig | mrnm %in% tgt_contig),] %>%
                       .[, ':='(pos = ifelse(rname == subseq_contig,  pos + subseq_start, pos ),
                                mpos = ifelse(mrnm == subseq_contig , mpos + subseq_start , mpos))] %>%
                       .[,':='(rname = ifelse(rname %in% c(tgt_contig, subseq_contig), new_scaffold, rname),
                               mrnm = ifelse(mrnm %in% c(tgt_contig, subseq_contig), new_scaffold, mrnm))]
                     
                     #rv$seq_len <- rv$seq_len %>%
                     #  rbind(.,list(new_scaffold, new_contig_len))
                     
                     rv$seq_len <- rv$contact_data2[,.(rlen = max(pos)), by = .(rname)] %>%
                       .[order(-rlen),]
                     
                     rv$choices <- rv$contact_data2[,.(rname, rlen = max(pos)), 
                                                    by = .(rname)][order(-rlen), rname]
                     choices_style <- list(style = 
                                             paste(rep_len("font-size: 12px; line-height: 1.5; margin-left: -10px; 
                                                                 border-bottom: 1px solid gray;", 
                                                           length(rv$choices)
                                             ), 
                                             "background-color: lightgray;")
                     )
                     
                     updatePickerInput(session, "seq", choices = rv$choices, 
                                       selected = new_scaffold,
                                       choicesOpt = choices_style)
                     updatePickerInput(session, "seq1", choices = rv$choices,
                                       choicesOpt = choices_style)
                     updatePickerInput(session, "subseq1", selected = "")
                     
                     #set plot size to default
                     rv$intramap_range <- NULL
                   })
    })
  })
  
  
  ###--------save changes to disk------------ 
  observeEvent(input$svChanges, {
    withBusyIndicatorServer("svChanges", {
      base_name <- sub("\\d+-\\d+-\\d+_(\\w+).rds|.bin|.txt*", "\\1",input$mapFile, 
                       ignore.case = TRUE)
      
      saveRDS(rv$contact_data2, 
              file.path(rv$projDir, paste0(Sys.Date(), "_", base_name, ".rds"))
      )
    })
  })
  
  ##------------------------------------------------------------------------------------------------------------
  ##------------------------------------------------------------------------------------------------------------
  observe({
    shiny::validate(need(input$seq2, ""))
    shiny::validate(need(rv$interseq_link_counts, ""))
    shiny::validate(need(input$nrSeq, ""))
    
    len <- length(isolate(input$joined_seqs))
    seq2 <- input$seq2
    leading_seq <- seq2[length(seq2)]
    leadingSeq_len <-  max(rv$contact_data2[rname == leading_seq, pos])
    maxLinkDen <- input$maxLinkDen
    minSeqLen <- input$minSeqLen * 1000000
    
    
    if(len > 1 & input$nrSeq > 1){
      if(leadingSeq_len < 20 * rv$binsize){
        leading_seq <- gsub("^[-+]", "",rv$s1)
      } else {
        leading_seq <- leading_seq
      }
    }
    
    if(input$dir2 == "Forward"){
      rv$subseq2 <- rv$interseq_link_counts[rname %in% leading_seq & 
                                        (!mrnm %in% leading_seq) & 
                                        link_density <= maxLinkDen &  mrnm_len >= minSeqLen,] %>%
        .[order(link_density, link_no, avg, decreasing = TRUE), 
          .(Subsequent_seq = mrnm, Length = mrnm_len, Strand = mrnm_strand, link_no, link_density)]
    } else {
      rv$subseq2 <- rv$interseq_link_counts[mrnm %in% leading_seq & 
                                        (!rname %in% leading_seq) & 
                                        link_density <= maxLinkDen & 
                                        rlen >= minSeqLen,] %>%
        .[order(link_density, link_no, avg, decreasing = TRUE), 
          .(Subsequent_seq = rname, Length = mrnm_len, Strand = rname_strand, link_no, link_density)] 
    }
    
  })
  
  observe({
    shiny::validate(need(rv$subseq2, ""))
    
    chr <- gsub("^[-+]", "", c(input$joined_seqs, input$scaf_man))
    rv$subseq2.1 <- rv$subseq2[(!(Subsequent_seq %in% chr)), ]
    rv$choices2 <- c(unique(rv$subseq2.1$Subsequent_seq), chr)
    
    updatePickerInput(session, "subseq2", choices =  rv$choices2,
                      choicesOpt = list(
                        style = paste(
                          rep_len("font-size: 12px; line-height: 1.5; 
                                  margin-left: -10px; border-bottom: 1px solid gray;", 
                                  length(rv$choices2)), "background-color: lightgray;")))
    shinyWidgets::updateSwitchInput(session, "strand_3", value = rv$subseq2.1$Strand[1] == "+")
  })
  
  ##----------------------------------------------------------------------------
  ##----------------------------------------------------------------------------
  output$Subsequent <- DT::renderDT({
    datatable(rv$subseq2,
              rownames = FALSE, class = 'display compact row-border', 
              selection = 'single', filter = 'bottom',
              options = list(
                pageLength = 15, dom = 'lti', 
                lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
                ScrollY = "500px",
                scrollX = TRUE,
                initComplete = JS(
                  "function(settings, json) {",
                  "$(this.api().table().header()).css({'background-color': '#F0F0F0', 'color': '#000'});",
                  "}")
              )
    )
  })
  ##----------------------------------------------------------------------------  
  ##---------------------------------------------------------------------------- 
  
  # move up and down the options list ------------------------------------------
  observeEvent(input$down1,{
    shiny::validate(need(rv$choices2, ""))
    
    i <- which(rv$choices2 == input$subseq2)
    i = i + 1
    
    if(i > length(rv$choices2)) return(NULL)
    s2 <- rv$choices2[i]
    strand_3 <- rv$subseq2.1[Subsequent_seq ==  s2, Strand][1]
    
    updatePickerInput(session, "subseq2",  selected = rv$choices2[i])
    shinyWidgets::updateSwitchInput(session, "strand_3", value = strand_3 == "+")
  })
  
  observeEvent(input$up1, {
    i <- which(rv$choices2 == input$subseq2)
    i = i - 1
    
    if(i < 1) return(NULL)
    s2 <- rv$choices2[i]
    strand_3 <- rv$subseq2.1[Subsequent_seq ==  s2, Strand][1]
    updatePickerInput(session, "subseq2",  selected = rv$choices2[i])
    shinyWidgets::updateSwitchInput(session, "strand_3", value = strand_3 == "+")
  })
  #-----------------------------------------------------------------------------
  
  scaffold_map <- reactive({
    shiny::validate(need(input$subseq2, ""))
    shiny::validate(need(input$seq2, ""))
    shiny::validate(need(input$nrSeq, ""))
    
    withBusyIndicatorServer("add", {
      strand_2 <- ifelse(input$strand_2, "+", "-")
      strand_3 <- ifelse(input$strand_3, "+", "-")
      rv$s2.1 <- paste0(strand_3, input$subseq2)
      chr <- isolate(input$joined_seqs)
      direction <- isolate(input$dir2)
      
      if(input$nrSeq == 1 | is.null(chr)){
        rv$s1 <- paste0(strand_2, input$seq2)
      } else {
        len <- length(chr)
        nrSeq <- min(len, input$nrSeq)
        
        if(len == 1){
          rv$s1 <- paste0(strand_2, input$seq2)
        } else {
          if(direction == "Forward"){
            rv$s1 <- chr[((len-nrSeq) + 1):len]
          } else {
            rv$s1 <- chr[1:nrSeq]
          }
        }
      }
      
      join_maps_plus(rv$contact_data2, seq = rv$s1, subseq = rv$s2.1, 
                     binsize = isolate(rv$binsize),  
                     direction = direction)
    })
  })
  
  
  
  ##----------------------------------------------------------------------------  
  ## Hi-C map-------------------------------------------- 
  observe({
    shiny::validate(need(scaffold_map(), ""))
    
    withBusyIndicatorServer("add", {
      
      # get total length
      total_len <- max(scaffold_map()[["pos"]])
      
      # get joining points
      len <- scaffold_map()[,.(len = min(pos)), by = .(rname)][len != 0, ]
      
      join_col <- if(nrow(len) == 1){"blue"} else {c(rep_len("gray", nrow(len)-1),"blue")}
      
      if(isolate(input$dir2) == "Backward"){
        join_col <- sort.default(join_col, decreasing = FALSE)
      }
      
      
      # plot
      rv$map1_plot <- ggplot(scaffold_map(), 
                             aes(x = pos, y = mpos, fill=log10(n/2))) +
        geom_tile(color = "red", size = 0.1) +
        scale_fill_gradient(low = "#ffe6e6", high = "red") +
        labs(title = isolate(rv$s2.1) , subtitle = 
               paste0("Total size: ", total_len/1000000, " Mb"), 
             x = "", y = "") +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        geom_vline(xintercept = len$len, color = join_col, size = 0.3) +
        geom_hline(yintercept = len$len, color = join_col, size = 0.3) +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(), 
              legend.position = "none",
              panel.border = element_rect(colour = "gray", fill = NA),
              panel.background = element_rect(fill = "white", colour = "white"))
    })
  })
  
  
  output$map1 <- renderPlot({
    shiny::validate(need(rv$map1_plot, ""))
    
    if(!is.null(rv$map1_range)){
      p <- rv$map1_plot + coord_cartesian(ylim = rv$map1_range, xlim = rv$map1_range) 
      session$resetBrush("map1_brush")
    } else {
      p <- rv$map1_plot
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
  
  observe({
    shiny::validate(need(scaffold_map(), ""))
    
    if(is.null(input$map1_brush)) 
      return(NULL)
    
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
    res <- nearPoints(scaffold_map(), input$map_hover, "pos", "mpos", threshold = 1, maxpoints = 1)
    if (nrow(res) == 0)
      return()
    res
  })
  
  ##----------------------------------------------------------------------------
  ## Join sequences together ---------------------------------------------------
  observeEvent(input$add,{
    withBusyIndicatorServer("add", {
      
      direction <- input$dir2
      
      if(length(input$joined_seqs) == 0){
        if (direction == 'Backward'){
          rv$chr <- c( rv$s2.1, rv$s1)
        } else {
          rv$chr <- c(rv$s1,  rv$s2.1)
        }
      } else {
        if (direction == 'Backward'){
          rv$chr <- c(rv$s2.1, input$joined_seqs)
        } else {
          rv$chr <-  c(input$joined_seqs, rv$s2.1)
        }
      }
      
      updateCheckboxGroupInput(session, "joined_seqs", NULL, 
                               choices = unique(rv$chr), selected = unique(rv$chr))
      updateTextAreaInput(session, "scaf_edit", NULL, 
                          value = paste(unique(rv$chr), collapse = "\n"))
    })
  })
  
  observe({
    chr <- input$joined_seqs
    
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
    shinyWidgets::updateSwitchInput(session, "strand_2", 
                                    value = substr(leading_seq, 1, 1) == "+", 
                                    disabled = disable_switch)
    updateNumericInput(session, "n", value = 1)
  })
  
  ##----------------------------------------------------------------------------
  ##----------------------------------------------------------------------------
  ## Copy to clipboard
  observeEvent(input$clipbtn, {
    clipr::write_clip(input$joined_seqs, 
                      allow_non_interactive = TRUE)
  })
  
  ## Erase the list 
  observeEvent(input$erase,{ 
    choices =  input$joined_seqs
    updateCheckboxGroupInput(session, "joined_seqs", NULL, 
                             choices =  choices[1], 
                             selected =  choices[1])
  })
  
  ## Edit the list manually 
  observeEvent(input$edit,{ 
    shinyjs::show("EditBox")
    hide("CheckBox")
    hide("edit")
    shinyjs::show("check")
  })
  
  observeEvent(input$check,{ 
    
    choices = base::strsplit(input$scaf_edit, "\n")[[1]]
    updateCheckboxGroupInput(session, "joined_seqs", NULL, 
                             choices =  choices, 
                             selected =  choices)
    
    shinyjs::hide("EditBox")
    shinyjs::show("CheckBox")
    shinyjs::show("edit")
    shinyjs::hide("check")
  })
  
  ## output groups to build fasta file using CreateScaffoldFasta.pl script
  observeEvent(input$export,{
    shiny::validate(need(input$joined_seqs, ""))
    
    out_dir <- file.path(rv$projDir,"groups")
    dir.create(out_dir, showWarnings = FALSE)

    rv$new_scaffold <- data.table(rname = sub("^[-+]","", input$joined_seqs),
                               rc = ifelse(grepl("^\\+", input$joined_seqs), 0, 1)) %>%
      .[rv$seq_len, on = "rname", nomatch=0] %>%
      .[, ':='(start = ifelse(grepl("fragment", rname), as.numeric(gsub(".*_(\\d+):.*", "\\1", rname)), 0),
          end = ifelse(grepl("fragment", rname), as.numeric(gsub(".*:(\\d+)", "\\1", rname)), rlen),
          contig = sub("_fragment.*", "", rname),
          id = seq_len(.N),
          name = "")] %>%
      .[, .(name, contig, fragment = rname, id, rc, start, end)]
    
    # Get next group number
    n <- length(list.files(out_dir))
    
    file = paste0(out_dir, "/group",n,".ordering")
    fwrite(rv$new_scaffold, file, quote = FALSE) 
  })
  
  observeEvent(input$export, {
    # Create fasta file
    hic_scaffolds <- list()
    
    for (s in unique(groups$name))
    {
      # init vars
      scaffold_seq <- NULL
      gap <- NULL
      gap_size <- 100
      
      # filter data
      data <- groups[groups$name == s,]
      
      for (i in 1:nrow(data))
      {
        contig_seq <- subseq(fasta[[data$contig[i]]], start = data$start[i], end = data$end[i])
        
        if(data$rc[i]){contig_seq <- reverseComplement(contig_seq)}
        
        scaffold_seq <- paste0(scaffold_seq, gap ,as.character(contig_seq))
        gap <- paste(rep("N", gap_size), collapse = '')
      }
      hic_scaffolds[[s]] <- scaffold_seq
    }
    
    hic_scaffolds <- DNAStringSet(sapply(hic_scaffolds, `[[`, 1))
    
    Biostrings::writeXStringSet(hic_scaffolds, "C:/Users/saeia/OneDrive - AgResearch/shared_with_andrew/T.repens_genome/hic_assembly.fasta")
  })
  
  
  
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # plot hi-c map for whole sequence

  observeEvent(input$combineMaps,{
    withBusyIndicatorServer("combineMaps",{
      withProgress(message = 'Preparing contact data', value = 1, 
                   detail = "please be patient ...", {
                     #maps <- if(input$inputType == "Manual"){base::strsplit(input$scaf_man, "\n")[[1]]} else {input$joined_seqs}
                     maps <- input$joined_seqs
                     
                     rv$combined_maps <- join_maps_plus(mat = rv$contact_data2, 
                                                        seq = maps, direction = "Forward", 
                                                        binsize = isolate(rv$binsize))
                     
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
    
    ggplot(rv$combined_maps, 
           aes(x = pos, y = mpos, fill=log10(n/2))) +
      geom_tile(color = "red", size = 0.1) +
      labs(title = paste("Total size:", max(len$len)/1000000, "MB"), x = "", y = "") +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      #geom_vline(xintercept = len$len	, color = "gray", size = 0.2) +
      #geom_hline(yintercept = len$len	, color = "gray", size = 0.2) +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(), 
            legend.position = "none",
            panel.border = element_rect(colour = "gray", fill = NA),
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
    
    res <- nearPoints(rv$combined_maps, input$map2_hover, "pos", "mpos", threshold = 1, maxpoints = 1)
    
    if(nrow(res) == 0)
      return()
    res
  })
}



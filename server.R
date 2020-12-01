
shinyOptions(progress.style="old")

source("./globals/functions.R")
source("./globals/helpers.R")
options(future.globals.maxSize = 50*1024^2,
        shiny.maxRequestSize=50*1024^2)
output_dir = NULL
#a <- readRDS("./data/bilberry/sequence_len_bilberry.rds")
#write.csv(a, "/Volumes/workspace/hrpazs/petunia_genome/HiC/HicTool/data/bilberry/sequence_len_bilberry.csv", row.names = FALSE)
#mat <- readRDS("./data/gillenia/Gillenia_contact.map_contigs_filt.rds")
#mat <- data.table::as.data.table(Gillenia_HiC_mat_5k_contigs_filt)
#mat <- Gillenia_HiC_mat_5k_contigs_filt
#mat <- as.data.table(readRDS("./data/gillenia/Gillenia_HiC_mat_5k_scaffolds_filt.rds"))
#mat <- contact.mat.5kb.noFilt
#sequence_length <- readRDS("./data/bilberry/sequence_len_bilberry.rds")
#mat <- as.data.table(Gillenia_contact.map_contigs_filt)
#sequence_length <- fread("./data/gillenia/contig_len.csv") 

#mat <- mat[!(mrnm %in% c("S4.1.26486", "S4.1.26547", "S4.1.26592", "S4.1.47370", "S4.1.262", "S4.1.47371", "S4.1.24488", "S4.1.225","S4.1.152" , "S4.1.27149", "S4.1.26973", "S4.1.27169", "S4.1.2973","S4.1.25215", "S4.1.139", "S4.1.26547", "S4.1.26592", "S4.1.47370", "S4.1.262")),]
#rmax = max(mat$pos)
#mmax = max(mat$mpos)
#max = max(rmax, mmax)
#binsize <- sort(unique(mat$pos))[1:2]
#binsize = as.integer(binsize[2] - binsize[1])


server <- function(input, output, session) {
  rv <- reactiveValues(chr = NULL, s2.1 = NULL)
  
  observeEvent(input$import_data,{
    withBusyIndicatorServer("import_data", {
      withProgress(message = 'Loading in progress',
                   detail = 'Please wait...', value = 1,{
                     rv$mat <- as.data.table(readRDS(input$mapPath)) 
                      # .[!(mrnm %in% c("S4.1.16096", "S4.1.25910","S4.1.27131","S4.1.27159","S4.1.27116","S4.1.26615","S4.1.153","S4.1.26088","S4.1.26615", "S4.1.26486", "S4.1.26547", "S4.1.26592", "S4.1.47370", "S4.1.262", "S4.1.47371", "S4.1.24488", "S4.1.225","S4.1.152" , "S4.1.27149", "S4.1.26973", "S4.1.27169", "S4.1.2973","S4.1.25215", "S4.1.139", "S4.1.26547", "S4.1.26592", "S4.1.47370", "S4.1.262")),]
                     
                     #sequence_length <- mat[,.(len = max(pos)), by = "rname"]
                     rv$sequence_length <- fread(input$lenPath, col.names = c("rname","len")) 
                     
                     rmax = max(rv$mat$pos)
                     mmax = max(rv$mat$mpos)
                     rv$max = max(rmax, mmax)
                     
                     binsize_ini <- sort.int(unique(rv$mat$pos), partial = 1:2)[1:2]
                     binsize_ini = as.integer(binsize_ini[2] - binsize_ini[1])
                     updateSliderInput(session, "binsize", min = binsize_ini, value = binsize_ini, step = binsize_ini , max = 20 * binsize_ini)
                     rv$binsize_ini <- binsize_ini
                   })
    })
  })
  
  observe({
    shiny::validate(need(input$binsize != 0, ""))
    
    withProgress(message = 'Binning in progress',
                 detail = 'This may take a while...', value = 1,{
                   rv$binsize <- input$binsize 
                   
                   if(rv$binsize_ini != input$binsize){
                     bins =  seq(rv$binsize , rv$max, rv$binsize)
                     #binsize = ifelse(max < binsize, max, binsize)
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
                   
                   edge_slc <- 10 * rv$binsize 
                   #system.time(
                   rv$edge_contact <- rv$sequence_length[rv$mat_binned[rname != mrnm,], on = c("rname")] %>%
                     .[rv$sequence_length, on = c("mrnm" = "rname"), ':='(len2 = i.len)] %>%
                     .[((pos < edge_slc | pos > len - edge_slc) & (mpos < edge_slc | mpos > len2 - edge_slc)),] %>%
                     .[, ':='(rname_dir = ifelse(pos < edge_slc, "-", 
                                                 ifelse(pos > len - edge_slc, "+", "M")),
                              mrnm_dir = ifelse(mpos < edge_slc, "+", 
                                                ifelse(mpos > len2 - edge_slc, "-", "M")),
                              edge_rname = ifelse(len < edge_slc, len, edge_slc),
                               edge_mrnm = ifelse(len2 < edge_slc, len2, edge_slc))] %>%
                     .[, .(no = .N, sum = sum(n), edge_rname = max(edge_rname), edge_mrnm = max(edge_mrnm)), .(rname, mrnm, rname_dir, mrnm_dir)] %>%
                     .[,.(edge_no = sum(edge_rname,edge_mrnm)/edge_slc), .(rname, mrnm, rname_dir, mrnm_dir, no)]
                 #  )
                   
                   updatePickerInput(session, "seq_1", choices = unique(rv$edge_contact$rname), selected = isolate(input$seq_1))
                 })
  })
  
  
  
  observe({
    shiny::validate(need(input$seq_1, ""))
    
    seq_1 <- input$seq_1
    leading_seq <- seq_1[length(seq_1)]
    
    if(input$dir == "Backward"){
      rv$sub_seq <- rv$edge_contact[mrnm == leading_seq,] %>%
        .[order(no, decreasing = TRUE), .(x = rname, x_dir = rname_dir)] 
    } else {
      rv$sub_seq <- rv$edge_contact[rname == leading_seq,] %>%
        .[order(no, decreasing = TRUE), .(x = mrnm, x_dir = mrnm_dir)] 
    }
  })
  
  observe({
    shiny::validate(need(rv$sub_seq, ""))
    
    chr <- gsub("^[-+]", "",c(input$scaf_auto, input$scaf_man))
    rv$sub_seq2 <- rv$sub_seq[(!(x %in% chr)), ]
    rv$choices <- c(unique(rv$sub_seq2$x), chr)
  
    updatePickerInput(session, "seq_2", choices =  rv$choices)
    shinyWidgets::updateSwitchInput(session, "orient_2", value = rv$sub_seq2$x_dir[1] == "+")
  })
  
  
  # move up and down the option list ---------------
  observeEvent(input$down,{
    i <- which(rv$choices == input$seq_2)
    i = i + 1
    
    if(i > length(rv$choices)) return(NULL)
    s2 <- rv$choices[i]
    orient_2 <- rv$sub_seq2[x ==  s2, x_dir][1]
    
    updatePickerInput(session, "seq_2",  selected = rv$choices[i])
    shinyWidgets::updateSwitchInput(session, "orient_2", value = orient_2 == "+")
  })
  
  observeEvent(input$up, {
    i <- which(rv$choices == input$seq_2)
    i = i - 1
    
    if(i < 1) return(NULL)
    s2 <- rv$choices[i]
    orient_2 <- rv$sub_seq2[x ==  s2, x_dir][1]
    updatePickerInput(session, "seq_2",  selected = rv$choices[i])
    shinyWidgets::updateSwitchInput(session, "orient_2", value = orient_2 == "+")
  })
  
  observe({
    shiny::validate(need(input$seq_2, ""))
    shiny::validate(need(input$seq_1, ""))
    shiny::validate(need(input$nrSeq, ""))
    withBusyIndicatorServer("build_scaffold", {
      
      orient_1 <- ifelse(input$orient_1, "+", "-")
      orient_2 <- ifelse(input$orient_2, "+", "-")
      rv$s2.1 <- paste0(orient_2, input$seq_2)
      
      
      chr <- isolate(input$scaf_auto)
      
      
      if(input$nrSeq == 1 | is.null(chr)){
        rv$s1 <- paste0(orient_1,input$seq_1)
      } else {
        len <- length(chr)
        nrSeq <- min(len, input$nrSeq)
        
        if(len == 1){
          rv$s1 <- paste0(orient_1,input$seq_1)
        } else {
          if(input$dir == "Forward"){
            rv$s1 <- chr[((len-nrSeq) + 1):len]
          } else {
            rv$s1 <- chr[1:nrSeq]
          }
        }
      }
      
      rv$sel_map <- JoinMultiMaps(rv$mat_binned, map1 = rv$s1, map2 = rv$s2.1, binsize = isolate(rv$binsize),  direction = isolate(input$dir), output_dir = output_dir, output = "data")
    })
  })
  
  output$map <- renderPlot({
    shiny::validate(need(rv$sel_map, ""))
    shiny::validate(need(input$seq_1, ""))
    
    len <- rv$sel_map[, .(len= max(pos)), by = .(rname)] %>%
      .[order(len), ]
    
    ggplot(rv$sel_map, aes(x = pos, y = mpos, fill=log10(n/2))) +
      geom_tile(color = "red", size = 0.2) +
      #scale_fill_gradient(low = "#ff9999", high = "red") +
      labs(title = isolate(rv$s2.1) , subtitle = paste0("Total size: ", max(len$len/1000000), " Mb"), x = "", y = "") +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      geom_vline(xintercept = len$len[1:(nrow(len)-1)]	, color = "gray", size = 0.3) +
      geom_hline(yintercept = len$len[1:(nrow(len)-1)], color = "gray", size = 0.3) +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(), 
            legend.position = "none",
            panel.border = element_rect(colour = "gray", fill = NA),
            panel.background = element_rect(fill = "white", colour = "white"))
  }, height = function() {
    session$clientData$output_map_width
  })
  
  
  ## Batch run ------------------------------------
  observeEvent(input$build_map,{
    shiny::validate(need(rv$sub_seq, ""))
    
    cmd <- paste0("open ", normalizePath(output_dir))
    system(cmd)
    
    direction <- input$dir
    s1 <- input$seq_1
    s2 <- rv$sub_seq$mrnm[1:20]
    
    #future({
    # draw the interaction map of the intial sequence with the other sequence individually and store plots on disk
    join2maps(mat, map1 = s1, map2 = s2, direction = direction, output_dir = output_dir)
    #})
  })
  
  
  ## Join sequences together ---------------------------
  observeEvent(input$build_scaffold,{
    withBusyIndicatorServer("build_scaffold", {
      
      if(length(input$scaf_auto) == 0){
        if (input$dir == 'Backward'){
          rv$chr <- c( rv$s2.1, rv$s1)
        } else {
          rv$chr <- c(rv$s1,  rv$s2.1)
        }
      } else {
        if (input$dir == 'Backward'){
          rv$chr <- c(rv$s2.1, input$scaf_auto)
        } else {
          rv$chr <-  c(input$scaf_auto, rv$s2.1)
        }
      }
      updateCheckboxGroupInput(session, "scaf_auto", NULL, choices = unique(rv$chr), selected = unique(rv$chr))
    })
    #updateTextInput(session, "chr", value = rv$chr)
  })
  
  
  observe({
    shiny::validate(need(input$scaf_auto,""))
    chr <- input$scaf_auto
    
    leading_seq <- ifelse(input$dir == 'Backward', chr[1], chr[length(chr)])
    #leading_seq_reverse <- ifelse(grepl("^-", leading_seq), TRUE, FALSE)
    #orient_1 <- substr(leading_seq, 1, 1)
    
    updatePickerInput(session, "seq_1", selected = gsub("^[-+]","",leading_seq))
    shinyWidgets::updateSwitchInput(session, "orient_1", value = substr(leading_seq, 1, 1) == "+")
    updateNumericInput(session, "n", value = 1)
  })
  
  
  
  # plot hi-c map for chr
  #------------------------
  observeEvent(input$combineMaps,{
    withBusyIndicatorServer("combineMaps",{
      withProgress(message = 'Preparing contact data', value = 1, 
                   detail = "please be patient ...", {
                     maps <- if(input$inputType == "Manual"){base::strsplit(input$scaf_man, "\n")[[1]]} else {input$scaf_auto}
                     
                     rv$combined_maps <- JoinMultiMaps(mat = rv$mat_binned, map1 = maps, direction = "Forward", output = "data", binsize = isolate(rv$binsize))
                     
                     #fwrite(chr, paste0("./","chr",chr_no,"_",n, ".csv"), row.names = FALSE)
                     #saveRDS(chr, paste0("chr1",chr_no,".rds"))
                     #chr <- fread(paste0("chr6-", 16,".csv"), stringsAsFactors = FALSE)
                   })
    })
  })
  
  output$combined_maps <- renderPlot({
    shiny::validate(need(rv$combined_maps, ""))
    
    len <- rv$combined_maps[, .(len= max(pos)), by = .(rname)] %>%
      .[order(len), ]
    
    ggplot(rv$combined_maps,aes(x = pos, y = mpos, fill=log10(n))) +
      geom_tile(color = "red", size = 0.2) +
      labs(title = paste("Total size:", max(len$len)/1000000, "MB"), x = "", y = "") +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      geom_vline(xintercept = len$len	, color = "gray", size = 0.2) +
      geom_hline(yintercept = len$len	, color = "gray", size = 0.2) +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(), 
            legend.position = "none",
            panel.background = element_rect(fill = "white", colour = "white"))
  }, height = function() {
    session$clientData$output_combined_maps_width
  })
}



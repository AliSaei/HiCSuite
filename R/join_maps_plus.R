#' @export

join_maps_plus <- function(mat, seq, subseq = NULL,  binsize, direction = "Forward", output_dir = NULL, output = list("data", "graph")){
  # name the chromosme being assembled
  
  
  maps <- gsub("^\\s+|\\s$", "", c(seq, subseq))
  n = length(maps)
  
  
  d <- setdiff(gsub("^-|^[+]","",maps),unique(mat$rname))
  if(length(d) > 0 ) stop(paste0('"', paste0(d, collapse = ","),'"', "does not exist"))
  
  
  if(n < 3) {
    combined_maps <- join_maps(mat, maps[1], maps[2], direction = direction , output_dir = output_dir, binsize = binsize, output = output)
  } else {
    if(direction == "Backward") maps <- gsub("^\\s+|\\s$", "", c(subseq, seq))
    
    combined_maps <- join_maps(mat, maps[1], maps[2], binsize = binsize,  output = "data")
    for(i in 3:n) {
      
      # remove '+-" signs
      scaffold <- maps[1:i]
      seq_name <- gsub('^-|^[+]',"", scaffold) 
      
      message(sprintf("Plotting %s",seq_name[i]))
      
      len <- combined_maps[, .(len = as.integer(min(pos)), len.m = as.integer(max(pos))), by = .(rname)] %>%
        .[order(len), ] 
      
      max_len <- as.integer(max(len$len.m) + binsize)
      
      sel_mat <-   mat[rname %in% seq_name & mrnm %in% seq_name,] %>% 
        .[(rname == seq_name[i] & mrnm %in% seq_name) | (mrnm == seq_name[i] & rname %in% seq_name),] %>%
        .[, ':='(pos = ifelse(rname == seq_name[i],  pos + max_len,
                              ifelse(rname %in% len$rname, pos + len$len[len$rname == rname], pos)),
                 mpos = ifelse(mrnm == seq_name[i],  mpos +  max_len ,
                               ifelse(mrnm %in% len$rname, mpos +  len$len[len$rname == mrnm], mpos))),  
          by = .(rname, pos, mrnm, mpos)]
      
      # check if scaffolds are rotated
      flipped.i = substring(scaffold[i],1,1) == "-"
      
      if(flipped.i){
        max.i <- as.integer(max(sel_mat[rname %in% seq_name[i], pos]))
        max.i0 <- max_len 
        
        sel_mat <- sel_mat[, ':='(pos = ifelse(rname == seq_name[i], (max.i-pos) + max.i0, pos),
                                  mpos = ifelse(mrnm == seq_name[i], (max.i-mpos) + max.i0, mpos)), by = .(rname, pos, mrnm, mpos)]
      }
      
      for(j in 1:(i-1)){
        flipped.j = substring(scaffold[j],1,1) == "-"
        if(flipped.j){
          
          sel_mat <- sel_mat[, ':='(pos = ifelse(rname == len$rname[j], (len$len.m[j]-pos) + len$len[j], pos),
                                    mpos = ifelse(mrnm ==  len$rname[j], (len$len.m[j]-mpos) + len$len[j], mpos)), by = .(rname, pos, mrnm, mpos)]
        }
      }
      combined_maps <- rbindlist(list(combined_maps, sel_mat))
    }
    
    if(output == "data"){
      return(combined_maps)
    } else {
      len <- combined_maps[, .(len= max(pos)), by = .(rname)] %>%
        .[order(len), ]
      
      p <- ggplot(combined_maps,aes(x = pos, y = mpos, fill=log10(n))) +
        geom_tile(color = "red", size = 0.2) +
        labs(title = paste("Total size:", max(len$len)/1000000, "MB"), x = "", y = "") +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        geom_vline(xintercept = len$len	, color = "blue", size = 0.1) +
        geom_hline(yintercept = len$len	, color = "blue", size = 0.1) +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(), 
              legend.position = "none",
              panel.background = element_rect(fill = "white", colour = "white"))
      
      if(!is.null(output_dir)){
        filename = paste0(n,i, ".tiff")
        ggsave(file.path(output_dir, seq, filename), scale = 1, width = 10, height = 10, units = "cm")
        return(NULL)
      } else {
        p
      }
    }
  }
}
join_maps <- function(mat, seq, subseq, binsize, direction = "Forward", output = "data", output_dir = NULL){
  
  if(is.null(seq)) 
    return(NULL)
  
  if(!is.null(output_dir)){
    
    # cretae a dir to store tiff files
    if(!dir.exists(file.path(output_dir, seq))){
      dir.create(file.path(output_dir, seq), showWarnings = FALSE, recursive = TRUE)
    }
  }
  
  seq <- seq[length(seq)]
  
  for(i in subseq){
    if(direction == "Backward"){
      maps <- gsub("^\\s+|\\s$", "", c(i, seq))
      ref_name = gsub('^-|^[+]',"",maps)
    } else {
      maps <- gsub("^\\s+|\\s$", "", c(seq, i))
      ref_name = gsub('^-|^[+]',"",maps)
    }
    
    max_ref1 <- max(mat[rname == ref_name[1], pos]) + binsize
    
    sel_mat <- mat[rname %in% ref_name & mrnm %in% ref_name,] %>%
      .[, ':='(pos = ifelse(rname == ref_name[2], (pos + max_ref1 ), pos),
               mpos = ifelse(mrnm == ref_name[2], (mpos + max_ref1), mpos))]
    
    
    if(substring(maps[1],1,1) == "-"){
      max_rf <- max(sel_mat[rname %in% ref_name[1], pos])
      sel_mat <- sel_mat[, ':='(pos = ifelse(rname == ref_name[1], max_rf-pos, pos),
                                mpos = ifelse(mrnm == ref_name[1], max_rf-mpos, mpos))]
    }
    
    if(substring(maps[2],1,1) == "-"){
      max_rf2 <- max(sel_mat[rname %in% ref_name[2], pos])
      sel_mat <- sel_mat[, ':='(pos = ifelse(rname == ref_name[2], max_rf2-pos, pos),
                                mpos = ifelse(mrnm == ref_name[2], max_rf2-mpos, mpos))] %>%
        .[, ':='(pos = ifelse(rname == ref_name[2], (pos + max_ref1), pos),
                 mpos = ifelse(mrnm == ref_name[2], (mpos + max_ref1), mpos))]
      
    } 
    
    if(output == "data"){
      return(sel_mat)
    } else {
      p <- ggplot(sel_mat,aes(x = pos, y = mpos, fill=log10(n/2))) +
        geom_tile(color = "red", size = 0.2) +
        #scale_fill_gradient(low = "#ff9999", high = "red") +
        labs(title = i, x = "", y = "") +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        geom_vline(xintercept = max(sel_mat[rname == ref_name[1], pos])	, color = "blue", size = 0.2) +
        geom_hline(yintercept = max(sel_mat[rname == ref_name[1], pos])	, color = "blue", size = 0.2) +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(), 
              legend.position = "none",
              panel.border = element_rect(colour = "gray", fill = NA),
              panel.background = element_rect(fill = "white", colour = "white"))
      
      if(!is.null(output_dir)){
        filename = paste0(n,i, ".tiff")
        ggsave(file.path(output_dir, seq, filename), scale = 1, width = 10, height = 10, units = "cm")
      } else {
        p
      }
    }
  }
}


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
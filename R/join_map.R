#' @export
join_map <- function(mat, seq, subseq, binsize, direction = "Forward", output = "data", output_dir = NULL){

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

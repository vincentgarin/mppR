# detect chunk
##############

# function to detect sequence of consistent genotype and set as missing the one
# length below a certain threshold.

detect_chunk <- function(d, thre = 5){
  
  for(j in 1:nrow(d)){
    
    gj <- d[j, ]
    gj <- gj[!is.na(gj)]
    
    if(length(gj) > 1){
      
      dif <- diff(gj)
      
      ref_id <- 1
      chunk <- c(ref_id, rep(NA, length(dif)))
      
      for(i in 1:length(dif)){
        
        if(dif[i] == 0){chunk[i+1] <- ref_id} else {chunk[i+1] <- ref_id + 1; ref_id <- ref_id + 1}
        
      }
      
      names(chunk) <- names(gj)
      
      # detect the chunk with less
      chunk_size <- table(chunk)
      rem_chunk <- names(chunk_size)[chunk_size <= thre]
      
      rem_mk <- names(chunk[chunk %in% as.numeric(rem_chunk)])
      
      d[j, colnames(d) %in% rem_mk] <- NA
      
    }
    
  }
  
  d
  
}
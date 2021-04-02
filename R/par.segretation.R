###################
# par.segregation #
###################

par.segretation <- function(x){
  
  # sub-function to test homozygosity
  
  al.sc.homo <- function(x) {
    
    strsplit(x, split = "")[[1]][1] == strsplit(x, split = "")[[1]][2]
    
  }
  
  miss <- is.na(x)
  
  if(sum(miss) == 2){return(1)
    
  } else if (sum(miss) == 1) {
    
    # test if the non-missing value is hetero or homozygous
    
    if(al.sc.homo(x = x[!miss])){return(2)} else {return(1)}
    
  }
  
  all.sc <- unique(x)
  
  if(length(all.sc) == 1) { return(1)
    
  } else {
    
    # test if homo/homo (AA/TT) or homo/hetero (AA/AT)
    
    test <- unlist(lapply(X = x, FUN = al.sc.homo))
    
    if (sum(!test) == 1) {return(2)} else{return(3)}
    
  }
  
  
}
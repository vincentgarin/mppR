#############
# allele.sc #
#############

# Function to determine the (homozygous) minor and major alleles and the
# heterozygous using SNP bi-allelic markers.

allele.sc <- function(x) {
  
  x <- x[!is.na(x)] # remove missing values
  
  if (length(x) == 0){ # Only missing values
    
    return(c(NA, NA, NA, NA))
    
  } else {
    
    # split all marker score into single allele scores
    
    alleles <- c(vapply(X = x, FUN = function(x) unlist(strsplit(x, split = "")),
                        FUN.VALUE = character(2)))
    
    alleles.sum <- table(alleles)
    
    if (length(alleles.sum) == 1) { # monomorphic case
      
      return(c(NA, NA, NA, NA))
      
    } else {
      
      all.sc <- unlist(attr(sort(alleles.sum, decreasing = FALSE), "dimnames"))
      all.sc <- c(paste0(all.sc[1], all.sc[1]), paste0(all.sc[2], all.sc[2]),
                  paste0(all.sc[1], all.sc[2]), paste0(all.sc[2], all.sc[1]))
      
      return(all.sc)
      
    }
    
  } 
  
}
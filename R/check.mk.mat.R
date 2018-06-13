################
# check.mk.mat #
################

# function to check if the format of the marker matrix is correct. It checks if:
# a) it is a matrix ; b) it is a character matrix ; c) it all marker scores are
# written with double letters.

check.mk.mat <- function(mk.mat){
  
  ### 1.1 check if it is a matrix or a data frame
  
  if (!(is.matrix(mk.mat))){
  
    stop("'mk.mat' is not a matrix")  
    
  }
  
  ### 1.2 check all columns are character
  
  if (!(is.character(mk.mat))){
    
    stop("'mk.mat' is not character")  
    
  }
  
  
  ### 1.3 Check that the format of marker scores is correct
  
  
    mk.scores <- attr(table(mk.mat), "dimnames")[[1]]
    test.format <- vapply(X = mk.scores, FUN = function(x) nchar(x) != 2,
                          FUN.VALUE = logical(1))
    
    if (sum(test.format) > 0) {
      
      pbmk <- paste(mk.scores[which(test.format)], collapse = ", ")
      
      message <- sprintf(ngettext(length(test.format),
                                  "The following marker score %s is not allowed",
                                  "The following marker scores %s are not allowed"),
                         pbmk)
      
      stop(message) 
      
    }    
  
}
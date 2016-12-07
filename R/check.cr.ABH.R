################
# check.cr.ABH #
################

# function to check the data that are introduced in cross.ABH function

check.cr.ABH <- function(par.sc, off.sc, cross.ind, par.per.cross){
  
  # check data
  ############
  
  # check if the parent and offspring matrix are characters
  
  if(!is.character(par.sc)){
    
    stop("The matrix of parent score par.sc argument is not
         character matrix.")
    
  }
  
  if(!is.character(off.sc)){
    
    stop("The matrix of offspring score off.sc argument is not
         character matrix.")
    
  }
  
  # check if the cross.ind is character
  
  if(!is.character(cross.ind)){
    
    stop("The cross.ind vector is not character vector.")
    
  }
  
  # check if the cross.indicator as the same length as the list of genotype
  
  if(length(cross.ind) != dim(off.sc)[1]){
    
    stop("The length of the cross.ind vector does not correspond
         to the number of genotype present in the genotype matrix.")
    
  }
  
  # check if the cross.ind is character
  
  if(!is.character(par.per.cross)){
    
    stop("The par.per.cross matrix is not a character matrix.")
    
  }
  
  # remove the eventual rownames of par.per.cross
  
  if(!is.null(rownames(par.per.cross))){
    
    rownames(par.per.cross) <- NULL
    
  }
  
  
  
  if (!identical(unique(cross.ind), par.per.cross[, 1])){
    
    stop("The cross identifiers used in cross.ind and in par.per.cross differ")
    
  }
  
  
  # check if the parent list is the same in the par.per.cross
  # and in the parents score rownames
  
  parents <- union(par.per.cross[,2], par.per.cross[,3])
  
  if(sum(!(parents %in% rownames(par.sc)))>0){
    
    list.par <- paste(parents[!(parents %in% rownames(par.sc))])
    
    stop("The following parents indicators: ", list.par, " (is) are present in
         par.per.cross object but not in the rownames of the par.sc matrix.")
    
  }
  
  if(sum(!(rownames(par.sc) %in% parents))>0){
    
    list.par <- paste(rownames(par.sc)[!(rownames(par.sc) %in% parents)])
    
    stop("The following parents indicators: ", list.par, " (is) are present in
          the rownames of the par.sc matrix but not in par.per.cross object.")
    
  }
  
}
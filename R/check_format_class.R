#################################
# check object format and class #
#################################

is_mppData <- function(x){
  
  inherits(x = x, what = "mppData")
  
}

# check trait

check_trait <- function(trait, mppData){
  
  if(!(is.character(trait) || is.numeric(trait))){
    
    stop("'trait' must be numeric or character")
    
  }
  
  if(is.numeric(trait)){
    
    nb.trait <- dim(mppData$pheno)[2]
    
    if(!((0 < trait) && (trait <=nb.trait))){
      
      stop("trait must be between 1 and ", nb.trait,
           " the total number of traits in 'mppData'")
      
    }
    
  }
  
  if(is.character(trait)){
    
    trait.names <- colnames(mppData$pheno)
    
    if (!(trait %in% trait.names)){
      
      t_nms <- paste(trait.names, collapse = ', ')
      
      stop("'trait' must be one of: ", t_nms)
      
    }
    
  }
  
}

# peak the trait
################

sel_trait <- function(mppData, trait){
  
  if(is.numeric(trait)){
    
    t_val <- mppData$pheno[, trait]
    
  } else {
    
    trait.names <- colnames(mppData$pheno)
    t_val <- mppData$pheno[, which(trait.names == trait)]
    
  }
  
  return(t_val)
  
}

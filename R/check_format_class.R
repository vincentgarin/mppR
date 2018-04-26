#################################
# check object format and class #
#################################

is_mppData <- function(x){
  
  inherits(x = x, what = "mppData")
  
}

# check trait

check_trait <- function(trait, gp){
  
  if(!(is.character(trait) || is.numeric(trait))){
    
    stop("The trait object must be either numeric or character")
    
  }
  
  if(is.numeric(trait)){
    
    nb.trait <- dim(gp$pheno)[3]
    
    if(!((0 < trait) && (trait <=nb.trait))){
      
      stop(paste("The trait indicator should be between 1 and", nb.trait,
                 "the total number of traits in the gp object"))
      
    }
    
  }
  
  if(is.character(trait)){
    
    trait.names <- attr(gp$pheno, "dimnames")[[2]]
    
    if (!(trait %in% trait.names)){
      
      stop(paste("trait must be one of:", trait.names))
      
    }
    
  }
  
}
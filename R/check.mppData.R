#################
# check.mppData #
#################

# function to check data format before the formation of mppData object


check.mppData <- function(geno, geno.par = NULL, type, type.mating, nb.gen,
                          biall, IBS.format, trait, map, cross.ind,
                          par.per.cross, dir){
  
  
  # 1. control format of genotype matrix
  ######################################
  
  if (!(is.matrix(geno))){
    
    stop("The marker matrix (geno) is not a matrix.")  
    
  }
  
  if((biall) & (IBS.format == "012")){
    
    if(!is.numeric(geno)) {
      
      stop("The marker matrix (geno) is not a numeric matrix.")
      
    }
    
  } else {
    
    if(!is.character(geno)) {
      
      stop("The marker matrix (geno) is not a character matrix.")
      
    }
    
  }
  
  
  
  mk.scores <- attr(table(geno), "dimnames")[[1]]
  
  if(!biall) { # cross, parental and ancestral model
    
    ### 1.1 Control the validity of the specified path
    
    if(!file.exists(dir)){
      
      stop("The path specified in the argument dir is not valid.")
      
    }
  
    if(sum(!(mk.scores %in% c("A", "B", "H", "-", NA))) > 0){
      
      prob.scores <- mk.scores[!(mk.scores %in% c("A", "B", "H", "-", NA))]
      
      info <- paste("The following marker score(s):",
                    paste(prob.scores, collapse = ", "),
                    "are not allowed")
      
      stop(info)
      
    }  
    
  } else { # bi-allelic model
   
     
    
  }
  
  # 2. control the type of population
  ###################################
  
  if (!(type %in% c("F", "bc", "RIL", "dh"))) {
    
    info <- paste("The type of population specified in the argument type:",
                   type, "is not allowed.", "Please use 'F', 'bc','RIL' or 'dh'")
    
    stop(info)
    
  }
  
  if ((type == "F") || (type == "bc")){
    
    if(is.null(nb.gen)){
      
      stop("The number of generation (nb.gen) is not specified.") }
    
  }
  
  if((type == "RIL") && is.null(type.mating)){
    
    stop("The type of mating (argument type.mating) must be provided.")
    
  }
  
  
  # 3. control format of the map
  ##############################
  
  if(!is.character(map[, 1])) {
    
    stop(paste("the marker list of the map argument is not a character vector.",
         "If if is a factor vector you can transform it using as.character()."))
    
  }
  
  if(sum(map[, 2] %% 1) > 0){
    
    stop(paste("chromosome indicators of the map argument are not all integers",
         "(1, 2, 3, ...)"))
    
  }
  
  if (!is.numeric(map[, 3])){
    
    stop("The marker position in centi-Morgan are not numeric.")
    
  }
  
  
  if( sum(diff(map[, 3]) == 0) > 0 ){
    
    warning("Some markers have the same position in the map.")
    
  }
  
  # 4. control format of the phenotypic data
  ##########################################
  
  if(!is.character(trait[, 1])) {
    
    stop(paste("The genotype identifiers of the trait argument are not",
              "a character vector. If if is a factor vector you can transform",
              "it using as.character()."))
    
  }
  
  if (!is.numeric(trait[, 2])){
    
    stop("The phenotypic values are not numeric.")
    
  }
  
  # 5. Control the similarity of marker list
  ##########################################
  
  
  if (length(colnames(geno)) != length(map[, 1])){
    
    stop("Genotype marker list and map data differ in length.")
  }
  
  # Test if the list are identical (same content, same order)
  
  if(!identical(colnames(geno), map[, 1])){
    
    stop(paste("Marker list in geno (colnames(geno)) and map (map[, 1]) are",
                "not identical"))
    
  }
  
  # 6. Control the similarity of genotype list
  ############################################
  
  
  if (length(rownames(geno)) != length(trait[, 1])){
    
    stop("Genotype data and trait data differ in length.")
  }
  
  if(!identical(rownames(geno), trait[, 1])){
   
    stop(paste("The genotypes list in geno (rownames(geno)) and trait",
               "(trait[, 1]) are not identical"))  
     
  }
  
  # 7. Test on cross.ind
  ######################
  
  if (length(cross.ind) != length(trait[, 1])){
    
    stop(paste("The length of the cross indicator is not the same as the",
                "one of the genotype list"))
    
  }
  
  # test par.per.cross
  
  if(!is.matrix(par.per.cross)){
    
    stop("The par.per.cross argument is not a matrix.")
    
  }
  
  if(!is.character(par.per.cross)){
    
    stop("The par.per.cross argument is not a character matrix.")
    
  }
  
  # remove the eventual rownames of par.per.cross
  
  if(!is.null(rownames(par.per.cross))){
    
    rownames(par.per.cross) <- NULL
    
  }
  
  
  if (!identical(unique(cross.ind), par.per.cross[, 1])){
    
    stop("The cross identifiers used in cross.ind and in par.per.cross differ")
    
  }
  
  if (!is.null(geno.par)){
    
    # 8. Test on geno.par
    #####################
    
    if(!identical(colnames(geno.par), map[, 1])){
      
      stop(paste("The marker list in geno.par (colnames(geno.par)) and map",
           "(map[, 1]) are not identical"))
      
    }
    
    # test similarity of parent indicators
    
    parents <- union(par.per.cross[,2], par.per.cross[,3])
    
    if(sum(!(parents %in% rownames(geno.par)))>0){
      
      list.par <- paste(parents[!(parents %in% rownames(geno.par))])
      
      stop(paste("The following parents indicators:", list.par,
                 "(is) are present in par.per.cross object but not in the",
                 "rownames of the geno.par matrix"))
      
    }
    
    if(sum(!(rownames(geno.par) %in% parents))>0){
      
      list.par <- paste(rownames(geno.par)[!(rownames(geno.par) %in% parents)])
      
      stop(paste("The following parents indicators:", list.par,
                 "(is) are present in the rownames of the geno.par matrix",
                 "but not in par.per.cross object"))
      
    }
    
  }
  
}
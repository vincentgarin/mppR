#############
# check_IBD #
#############

# Function to check the argument provided to the function IBS.mppData which
# convert the genotype data into 0, 1, 2 format and do an optional
# genotype imputation.

check_IBD <- function(mppData, het.miss.par, subcross.ind, par.per.subcross,
                      type, F.gen, BC.gen, type.mating, map.function){
  
  # 1. check mppData format
  #########################
  
  if(!is_mppData(mppData)){
    
    stop('the mppData provided provided is not a mppData object.')
    
  }
  
  # test if correct step in the mppData processing
  
  if(mppData$status != 'IBS'){
    
    stop(paste('You have to process the mppData objects in a strict order:',
               'create.mppData(), QC.mppData(), IBS.mppData(), IBD.mppData(),',
               'parent_cluster.mppData(). You can only use IBD.mppData()',
               'after performing create.mppData(), QC.mppData(), and',
               'IBS.mppData()'))
    
  }
  
  geno.par <- mppData$geno.par
  
  # 2. check ABH coding arguments
  ###############################
  
  # Test if there are some heterozygous parents and the user want to make the
  # ABH assignement without specifying het_miss_par = TRUE
  
  if(!het.miss.par){
    
    # test if there are some heterozygous parents
    
    par.het <- QC_hetero(mk.mat = geno.par)
    
    het.mk <- which(par.het != 0)
    
    if(length(het.mk) > 0){
      
      message <- paste("The following markers:",
                       paste(colnames(geno.par)[het.mk], collapse = ", "),
                       "are heterozygous for at least one parent. In order",
                       "to proceed to the ABH assignement either remove them",
                       "or use option het.miss.par = TRUE.")
      
      stop(message)
      
    }
    
  }
  
  # test par.subcross
  
  if((!is.null(subcross.ind)) || (!is.null(par.per.subcross))){
    
    if((!is.null(subcross.ind)) && (is.null(par.per.subcross))){
      
      stop("You must also provide the par.per.subcross argument.")
      
    }
    
    if((!is.null(par.per.subcross)) && (is.null(subcross.ind))){
      
      stop("You must also provide the subcross.ind argument.")
      
    }
    
    
    if(!is.matrix(par.per.subcross)){
      
      stop("The par.per.subcross argument is not a matrix.")
      
    }
    
    if(!is.character(par.per.subcross)){
      
      stop("The par.per.subcross argument is not a character matrix.")
      
    }
    
    # remove the eventual rownames of par.per.subcross
    
    if(!is.null(rownames(par.per.subcross))){
      
      rownames(par.per.subcross) <- NULL
      
    }
    
    if (!identical(unique(subcross.ind), par.per.subcross[, 1])){
      
      stop(paste("The subcross identifiers used in subcross.ind and",
                 "in par.per.subcross differ"))
      
    }
    
    # test the similarity of parents list between par.per.subcross and
    # rownames(geno.par)
    
    parents <- union(par.per.subcross[,2], par.per.subcross[,3])
    
    if(sum(!(parents %in% rownames(geno.par)))>0){
      
      list.par <- paste(parents[!(parents %in% rownames(geno.par))])
      
      stop(paste("The following parents indicators:", list.par,
                 "(is) are present in par.per.subcross object but not in the",
                 "rownames of the geno.par matrix"))
      
    }
    
    
  }
  
  
  # 3. check argument for IBD computation
  #######################################
  
  
  if (!(type %in% c("F", "BC", "RIL", "DH", "BCsFt"))) {
    
    info <- paste("The type of population specified in the argument type:",
                  type, "is not allowed.",
                  "Please use 'F', 'BC','RIL','DH' or 'BCsFt'.")
    
    stop(info)
    
  }
  
  ### Check specification of generation number for BC F and BCsFt populations
  
  if (type == "F"){
    
    if(is.null(F.gen)){
      
      stop("The number of generation (F.gen) is not specified.") }
    
  }
  
  if (type == "BC"){
    
    if(is.null(BC.gen)){
      
      stop("The number of generation (BC.gen) is not specified.") }
    
  }
  
  if (type == "BCsFt"){
    
    if(is.null(BC.gen)|| is.null(F.gen)){
      
      stop("The number of generation (BC.gen or F.gen) is/are not specified.") }
    
  }
  
  
  
  if((type == "RIL") && is.null(type.mating)){
    
    stop("The type of mating (argument type.mating) must be provided.")
    
  }
  
  # type.mating
  
  if(!is.null(type.mating)){
    
    if(!is.character(type.mating)){stop('type.mating must be a character string.')}
    
    if(!(type.mating %in% c('selfing', 'sib.mat'))){
      
      stop("type.mating must be either 'selfing' or 'sib.mat'.")
      
    }
    
  }
  
  # map.function
  
  if(!is.character(map.function)){
    
    stop('map.function must be a character string.')
    
  }
  
  if(!(map.function %in% c('haldane', 'kosambi', 'c-f', 'morgan'))){
    
    stop("type.mating must be one of: 'haldane', 'kosambi', 'c-f', 'morgan'.")
    
  }
  
}
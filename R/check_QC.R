############
# check_QC #
############

# check the data format for the quality control procedure.



check_QC <- function(geno.off, geno.par, map, trait, cross.ind, par.per.cross,
                     subcross.ind, par.per.subcross, n.lim, MAF.cr.lim, ABH,
                     het.par, parallel, cluster){
  
  # test the value of minimum cross size
  
  if(n.lim < 15){
    
    stop("It is not allowed/adviced to use cross with less than 15 observations.")
    
  }
  
  # test format of the marker matrices
  
  if(!is.matrix(geno.off)){
    
    stop("The geno.off argument is not a matrix.")
    
  }
  
  if(!is.character(geno.off)){
    
    stop("The geno.off argument is not a character matrix.")
    
  }
  
  if(!is.matrix(geno.par)){
    
    stop("The geno.par argument is not a matrix.")
    
  }
  
  if(!is.character(geno.par)){
    
    stop("The geno.par argument is not a character matrix.")
    
  }
  
  # test if the marker identifiers are the same in the map and in the marker
  # matrices
  
  if(!identical(colnames(geno.off), map[, 1])){
    
    stop(paste("The marker identifier of the offspring marker matrix (geno.off)",
               "and the list of marker (map[, 1]) are not the same."))
    
  }
  
  if(!identical(colnames(geno.par), map[, 1])){
    
    stop(paste("The marker identifier of the parent marker matrix (geno.par)",
               "and the list of marker (map[, 1]) are not the same."))
    
  }
  
  # test if the list of genotype is the same between the offspring marker matrix
  # and the phenotypic values.
  
  if(!identical(rownames(geno.off), trait[, 1])){
    
    stop(paste("The genotype identifiers of the offspring marker matrix",
               "(geno.off) and the list of genotypes (trait[, 1]) are not",
               "the same."))
    
  }
  
  # length cross.ind same as list of genotypes
  
  if(length(cross.ind) != dim(geno.off)[1]){
    
    stop(paste("The cross indicator vector (cross.ind) length does not have",
               "the same length as the genotype list (dim(geno.off)[1])"))
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
    
    stop(paste("The cross identifiers used in cross.ind and",
               "in par.per.cross differ"))
    
  }
  
  # test the similarity of parents list between par.per.cross and
  # rownames(geno.par)
  
  parents <- union(par.per.cross[,2], par.per.cross[,3])
  
  if(sum(!(parents %in% rownames(geno.par)))>0){
    
    list.par <- paste(parents[!(parents %in% rownames(geno.par))])
    
    stop(paste("The following parents indicators:", list.par,
               "(is) are present in par.per.cross object but not in the",
               "rownames of the geno.par matrix"))
    
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
  
  # test if there are as many within cross MAF as cross in the MAF.cr.lim argument
  
  if(!is.null(MAF.cr.lim)){
    
    if(length(MAF.cr.lim) != length(unique(cross.ind))){
      
      stop("MAF.cr.lim must contain one value per cross.")
      
    }
    
  }
  
  # Test if there are some heterozygous parents and the user want to make the
  # ABH assignement without specifying het.par = TRUE
  
  if(ABH && !het.par){
    
    # test if there are some heterozygous parents
    
    par.het <- QC_hetero(mk.mat = geno.par, parallel = parallel,
                         cluster = cluster)
    
    het.mk <- which(par.het != 0)
    
    
    if(length(het.mk) > 0){
      
      message <- paste("The following markers:",
                       paste(colnames(geno.par)[het.mk], collapse = ", "),
                       "are heterozygous for at least one parent. In order",
                       "to proceed to the ABH assignement either remove them",
                       "or use option het.par = TRUE.")
      
      stop(message)
      
    }
    
  }
  
  # test the argument necessary for the ABH assignement par.per.cross and
  # cross.ind.ABH.
  
  
  # test if cluster have been correctly build by the user for parallelization
  
  if (parallel){
    
    if(!inherits(cluster, "cluster")){
      
      stop(paste("You must provide cluster objects to compute this function",
                 "in parallel. Use function makeCluster() from parallel pakage."))
      
    }
    
  }
  
}
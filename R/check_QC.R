############
# check_QC #
############

# check the data format for the quality control procedure.



check_QC <- function(geno.off, geno.par, map, trait, cross.ind, par.per.cross,
                     subcross.ind, par.per.subcross, n.lim, MAF.cr.lim, ABH,
                     het.par, impute, impute.type, map_bp, replace.value,
                     parallel, cluster){
  
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
  
  # test if the user want to perform both ABH assignement and imputation. Not
  # possible.
  
  if((ABH & impute)){
    
    stop(paste("It is not possible to select both ABH assignement and",
               "imputation. Imputation option is made for IBS data that will",
               "be used for the bi-allelic model. For the other (IBD) models,",
               "some imputation of the missing value will be performed during",
               "the IBD value computation in function mppData_form().",
               "Therefore, if you want to get data for IBD models",
               "(cross-specific, parental or ancestral), keep ABH = TRUE,",
               "and set impute = FALSE.", "If you want, to get imputed IBS data",
               "for a bi-allelic, set ABH = FALSE and keep impute = TRUE."))
    
  }
  
  # Test the validity of the arguments introduced if the user want to perform
  # some imputation
  
  if(impute){
    
    if(!(impute.type %in% c("random","family","beagle","beagleAfterFamily",
                            "beagleNoRand", "beagleAfterFamilyNoRand","fix"))){
      
      stop(paste("The impute.type is not correct. Please use one of:",
                 "'random', 'family', 'beagle', 'beagleAfterFamily',",
                 "'beagleNoRand', 'beagleAfterFamilyNoRand', 'fix'."))
      
    }
    
    # check if impute is beagle the right format of the map_bp
    
    if (impute.type %in% c("beagle","beagleAfterFamily",
                           "beagleNoRand", "beagleAfterFamilyNoRand","fix")){
      
      if(!("package:synbreed" %in% search())){
        
        stop(paste("To be able to use Beagle for imputation, please load the",
                   "package synbreed using library(synbreed)."))
        
      }
      
      if(!identical(map_bp[, 1], colnames(geno.off))){
        
        stop(paste("The marker indicators of map_bp should be identical to the",
                   "the marker indicators of geno.off (colnames(geno.off))."))
        
      }
      
      if(!(is.character(map_bp[, 2]) || is.numeric(map_bp[, 2]))){
        
        stop(paste("The chromosome indicator of map_bp (2nd col) should be",
                   "character or numeric."))
        
      }
      
    }
    
    
    # if impute.type is fix. Chekc that a value was provided for replace.value
    
    if((impute.type == "fix") & (is.null(replace.value))){
      
      stop(paste("To use the impute.type = 'fix', you must also provide a value",
                 "the replace.value argument."))
           
    }
    
    if((impute.type == "fix")){
      
      if(!is.numeric(replace.value)){
        
        stop("replace.value must be a numerical value 0, 1 or 2.")
        
      }
      
      if (!(replace.value %in% c(0, 1, 2))){
        
        stop("replace.value must be 0, 1 or 2.")
        
      }
      
    }
    
  }
  
  # test if cluster have been correctly build by the user for parallelization
  
  if (parallel){
    
    if(!inherits(cluster, "cluster")){
      
      stop(paste("You must provide cluster objects to compute this function",
                 "in parallel. Use function makeCluster() from parallel pakage."))
      
    }
    
  }
  
}
#############
# check_IBS #
#############

# Function to check the argument provided to the function IBS.mppData which
# convert the genotype data into 0, 1, 2 format and do an optional
# genotype imputation.

check_IBS <- function(mppData, impute, impute.type, map_bp, replace.value){
  
  # check mppData format
  
  if(!is_mppData(mppData)){
    
    stop(paste('the mppData provided provided is not a mppData object.',
               'Please use function create.mppData().'))
    
  }
  
  # test if correct step in the mppData processing
  
  if(mppData$status != 'QC'){
    
    stop(paste('You have to process the mppData objects in a strict order:',
               'create.mppData(), QC.mppData(), IBS.mppData(), IBD.mppData(),',
               'parent_cluster.mppData(). You can only use IBS.mppData()',
               'after performing create.mppData() and QC.mppData().'))
    
  }
  
  if(impute){
    
    if(!(impute.type %in% c("random","family","beagle","beagleAfterFamily",
                            "beagleNoRand", "beagleAfterFamilyNoRand","fix"))){
      
      stop(paste("The impute.type is not correct. Please use one of:",
                 "'random', 'family', 'beagle', 'beagleAfterFamily',",
                 "'beagleNoRand', 'beagleAfterFamilyNoRand', 'fix'."))
      
    }
    
    # check if impute is beagle the right format of the map_bp
    
    if (impute.type %in% c("beagle","beagleAfterFamily",
                           "beagleNoRand", "beagleAfterFamilyNoRand")){
      
      if(!("package:synbreed" %in% search())){
        
        stop(paste("To be able to use Beagle for imputation, please load the",
                   "package synbreed using library(synbreed)."))
        
      }
      
      if(is.null(map_bp)){
        
        stop(paste("to use imputation with Beagle you must provide marker bp",
                   "position via the argument map_pb."))
        
      }
      
      if(!identical(map_bp[, 1], colnames(mppData$geno.off))){
        
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

  
}
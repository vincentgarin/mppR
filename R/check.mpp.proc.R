##################
# check.mpp.proc #
##################

# function to check the format of all elements introduced in the function
# mpp.proc

check.mpp.proc <- function(mppData, trait, Q.eff, VCOV, plot.gen.eff = FALSE,
                           ref.par = NULL, sum_zero = NULL, n.cores = 1,
                           output.loc){
  
  # 1. test the validity of the provided path to store the results
  
  if(!file.exists(output.loc)){
    
    stop("The path specified in the argument output.loc is not valid.")
    
  }
  
  # 2. check mppData format
  
  
  if(!inherits(mppData, "mppData")) {
    
    stop("The data object provided (argument mppData) is not of class mppData.")
    
  }
  
  # 3. check trait format
  
  check_trait(trait = trait, mppData = mppData)
  
  # 4. check Q.eff argument
  
  
  if (!(Q.eff %in% c("cr", "par", "anc", "biall"))){
    
    stop("The Q.eff argument must be : 'cr', 'par', 'anc' or 'biall'.")
    
  }
  
  # 5. check the VCOV argument
  
  
  if (!(VCOV %in% c("h.err", "h.err.as", "cr.err", "pedigree", "ped_cr.err"))){
    
    stop(paste("The VCOV argument must be : 'h.err', 'h.err.as', 'cr.err',",
               "'pedigree' or 'ped_cr.err'."))
    
  }
  
  # 6. test if the asreml function is present for the compuation of the
  # mixed models
  
    
    if (VCOV != "h.err"){
      
      test <- requireNamespace(package = 'asreml', quietly = TRUE)
      
      if(!test){
        
        stop(paste("To use this type of VCOV, you must have access to the asreml",
                   "function from the asreml-R package."))
        
      }
      
    }
    
  
  # 7. Consistency for parallelization.
  
  
  if ((n.cores > 1) && (VCOV != "h.err")){
    
    stop("Parallelization is only allowed for VCOV = 'h.err'.") 
    
    
  }
  
  
    # 8. Test that if the user wants to fit a bi-allelic model, the
    # plot.gen.eff is not activated
    
    if((Q.eff == "biall") && plot.gen.eff) {
      
      stop("The estimation of the decomposed QTL effect (plot.gen.eff = TRUE)
           per cross or parents can not be performed for the bi-allelic model")
      
    }
  
  # 9. Check the argument ref.par
  
  if(!is.null(ref.par)){
    
    # test that there is only one reference parent and one connected part.
    
    if(length(ref.par) !=1){
      
      stop(paste("You can only specify one reference parent (ref.par) for",
                 "the estimation of the QTL effects."))
      
    }
    
    nb.con.part <- length(design_connectivity(mppData$par.per.cross,
                                              plot_des = FALSE))
    
    if(nb.con.part > 1){
      
      stop(paste("You can only use the ref.par argument if your MPP design",
                 "is composed of a single connected part",
                 "(check with design_connectivity(mppData$par.per.cross))."))
      
    }
    
    # test that reference parent is present in the list of parents.
    
    if(!(ref.par %in% mppData$parents)){
      
      stop(paste("The reference parent you specified in ref.par is not",
                 "contained in the list of parents. Please use one of:",
                 paste(mppData$parents, collapse = ", ")))
      
    }
    
  }
  
  # 10. Check the argument sum_zero
  
  if(!is.null(sum_zero)){
    
    if(sum_zero){
      
      if(!(Q.eff %in% c('par', 'anc'))){
        
        stop(paste('You can only use the sum to zero constraint for the',
                   'parental (Q.eff = "par") or the ancestral (Q.eff = "anc")',
                   'models.'))
        
      }
      
      if(VCOV != 'h.err'){
        
        stop(paste('You can only use the sum to zero constraint with the',
                   'homogeneous error term model (VCOV = "h.err").'))
        
      }
      
    }
    
    
  }
    
}
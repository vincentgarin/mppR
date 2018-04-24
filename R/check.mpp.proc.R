##################
# check.mpp.proc #
##################

# function to check the format of all elements introduced in the function
# mpp.proc

check.mpp.proc <- function(mppData, Q.eff, VCOV, par.clu = NULL,
                           plot.gen.eff = FALSE, ref.par = NULL,
                           sum_zero = NULL, parallel = FALSE, cluster,
                           output.loc){
  
  # 1. test the validity of the provided path to store the results
  
  if(!file.exists(output.loc)){
    
    stop("The path specified in the argument output.loc is not valid.")
    
  }
  
  # 2. check mppData format
  
  
  if(!inherits(mppData, "mppData")) {
    
    stop("The data object provided (argument mppData) is not of class mppData.")
    
  }
  
  # 3. check Q.eff argument
  
  
  if (!(Q.eff %in% c("cr", "par", "anc", "biall"))){
    
    stop("The Q.eff argument must be : 'cr', 'par', 'anc' or 'biall'.")
    
  }
  
  # 4. check the VCOV argument
  
  
  if (!(VCOV %in% c("h.err", "h.err.as", "cr.err", "pedigree", "ped_cr.err"))){
    
    stop(paste("The VCOV argument must be : 'h.err', 'h.err.as', 'cr.err',",
               "'pedigree' or 'ped_cr.err'."))
    
  }
  
  # 5. test if the asreml function is present for the compuation of the
  # mixed models
  
    
    if (VCOV != "h.err"){
      
      test <- requireNamespace(package = 'asreml', quietly = TRUE)
      
      if(!test){
        
        stop(paste("To use this type of VCOV, you must have access to the asreml",
                   "function from the asreml-R package."))
        
      }
      
    }
    

  
  
  # 6. consistency between Q.eff and the type of mppData object
  
  if ((Q.eff=="cr" && mppData$biall) | (Q.eff=="par" && mppData$biall) |
      (Q.eff=="anc" && mppData$biall)){
    
    stop(paste("The mppData object is made for bi-allelic models and not for",
               "cross, parental or ancestral models."))
    
  }
  
  if((Q.eff=="biall" & !mppData$biall)){
    
    stop(paste("The mppData object is made for cross, parental or ancestral",
               "models and not for bi-allelic models."))
    
  }
  
  # 7. Consistency for parallelization.
  
  
  if (parallel){
    
    if(parallel & (!(VCOV == "h.err"))) {
      
      stop("Parallelization is only allowed for VCOV = 'h.err'.") 
      
    }
    
    if(!inherits(cluster, "cluster")){
      
      stop(paste("You must provide cluster objects to compute this function",
                 "in parallel. Use function makeCluster() from parallel pakage."))
      
    }
    
  }
  
  
  
  # 8. if ancestral model, check the format of the par.clu object
  
  
  if(Q.eff == "anc"){
    
    if(is.null(par.clu)){
      
      stop(paste("You need to provide a parent clustering object",
                 "(argument par.clu) for the computation of an ancestral model."))
      
    }
    
    if(!is.integer(par.clu)){
      
      stop("The par.clu argument is not and integer matrix.")
      
    }
    
    if(!identical(mppData$map[, 1], rownames(par.clu))){
      
      stop(paste("The list of markers and in between positions of the par.clu",
                 "object is not the same as the one in the mppData object map."))
      
    }
    
    if(!identical(sort(mppData$parents), sort(colnames(par.clu)))) {
      
      stop("The list of parents of the par.clu object (colnames(par.clu)) is
           not the same as the one of the mppData object.")
      
    }
    
  }
  
  
    # 9. Test that if the user wants to fit a bi-allelic model, the
    # plot.gen.eff is not activated
    
    if((Q.eff == "biall") && plot.gen.eff) {
      
      stop("The estimation of the decomposed QTL effect (plot.gen.eff = TRUE)
           per cross or parents can not be performed for the bi-allelic model")
      
    }
  
  # 10. Check the argument ref.par
  
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
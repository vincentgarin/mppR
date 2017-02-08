################
# check.mpp.cv #
################

# function to check the format of all elements introduced in the function
# mpp_CV


check.mpp.cv <- function(mppData, Q.eff, VCOV, par.clu = NULL,
                         parallel = FALSE, cluster, output.loc){
  
  
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
    
    stop("The Q.eff argument must be : 'cr', 'par', 'anc' or 'biall'")
    
  }
  
  # 4. check the VCOV argument
  
  
  if (!(VCOV %in% c("h.err", "h.err.as", "cr.err", "pedigree", "ped_cr.err"))){
    
    stop(paste("The VCOV argument must be : 'h.err', 'h.err.as', 'cr.err',",
               "'pedigree' or 'ped_cr.err'."))
    
  }
  
  # 5. test if the asreml function is present for the compuation of the
  # mixed models
  
  
  if (VCOV != "h.err"){
    
    if(!((exists("asreml")) && (is.function(asreml)))){
      
      stop(paste("To use this type of VCOV, you must have access to the asreml",
                 "function from the asreml-R package."))
      
    }
    
    message(paste("MESSAGE: Cross-validation",
                  "using mixed models is technically possible but can take",
                  "a hugh amount of time! Due to limited computer power, we",
                  "advice to only use the linear model (VCOV = 'h.err')."))
    
    text <- paste("Press [enter] to continue. If after you still want to stop",
                  "the process, Press Esc.")
    
    readline(prompt = text)
    
  }
  
  
  # 6. consistency between Q.eff and the type of mppData object
  
  if ((Q.eff=="cr" && mppData$biall) | (Q.eff=="par" && mppData$biall) |
      (Q.eff=="anc" && mppData$biall)){
    
    stop(paste("the mppData object is made for bi-allelic models and not for",
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
      
      stop(paste("The list of marker and inbetween position of the par.clu",
                 "object is not the same as the one in the mppData object map."))
      
    }
    
    if(!identical(sort(mppData$parents), sort(colnames(par.clu)))) {
      
      stop(paste("The list of marker and inbetween position of the par.clu",
                 "object is not the same as the one in the mppData object map."))
      
    }
    
  }
  
  
}
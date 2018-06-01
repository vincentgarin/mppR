################
# check.mpp.cv #
################

# function to check the format of all elements introduced in the function
# mpp_CV


check.mpp.cv <- function(mppData, trait, Q.eff, VCOV, n.cores = 1, output.loc,
                         her){
  
  
  # 1. test the validity of the provided path to store the results
  
  if(!file.exists(output.loc)){
    
    stop("The path specified in the argument output.loc is not valid.")
    
  }
  
  # 2. check mppData format
  
  
  if(!inherits(mppData, "mppData")) {
    
    stop("The data object provided (argument mppData) is not of class mppData.")
    
  }
  
  if(mppData$status != 'complete'){
    
    stop(paste('The mppData object is not complete. You must use all processing',
               'functions first in the specified order: QC.mppData, IBS.mppData,',
               'IBD.mppData, and parent_cluster.mppData.'))
    
  }
  
  # 3. check the trait
  
  check_trait(trait = trait, mppData = mppData)
  
  # 4. check Q.eff argument
  
  
  if (!(Q.eff %in% c("cr", "par", "anc", "biall"))){
    
    stop("The Q.eff argument must be : 'cr', 'par', 'anc' or 'biall'")
    
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
    
    message(paste("MESSAGE: Cross-validation",
                  "using mixed models is technically possible but can take",
                  "a hugh amount of time! Due to limited computer power, we",
                  "advice to only use the linear model (VCOV = 'h.err')."))
    
    text <- paste("Press [enter] to continue. If after you still want to stop",
                  "the process, Press Esc.")
    
    readline(prompt = text)
    
  }
  
  # 7. Consistency for parallelization.
  
  
  if ((n.cores > 1) && (VCOV != "h.err")){
    
    stop("Parallelization is only allowed for VCOV = 'h.err'.") 
      
    
  }
  
  
  # 8. check the format of her argument
  
  if((length(her) != 1) & (length(her) != mppData$n.cr)){
    
    stop(paste("the her argument must either be of lenght one representing the",
                "average heritability over cross or of length n.cr representing",
                "the within cross heritabilities"))
    
  }
  
  
}
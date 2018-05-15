#############
# check.MQE #
#############

# function to check the format of the data provided to MQE functions.

# parameters

# mppData:  objects of class mppData.

# Q.eff Character vector specifying the type of QTL effect.

# VCOV Character specify the variance covariance structure.

# n.cores: number of cluster

# cluster cluster to run the function in parallel


check.MQE <- function(mppData = NULL, trait, Q.eff, VCOV, cofactors = NULL,
                      cof.Qeff, n.cores = 1, QTL = NULL, output.loc,
                      fct = "XXX"){
  
  # 1. check mppData format
  #########################
  
  if(is.null(mppData)) {
    
    stop("You must provide a mppData object in argument mppData.")
    
  }
  
  if(!inherits(mppData, "mppData")) {
    
    stop("The object provided in mppData is not of class mppData.")
    
  }
  
  # 2. check trait
  ################
  
  check_trait(trait = trait, mppData = mppData)
  
  # 2. check Q.eff argument
  #########################
  
  if((fct == "proc") | (fct == "forward")){
    
    if (length(Q.eff) <= 1) {
      message("You should provide at least 2 type of QTL effect for the argument Q.eff.")
    }
    
  }
  
  
  test.Qeff <- Q.eff %in% c("cr", "par", "anc", "biall")
  
  if(sum(!test.Qeff) != 0 ){
    
    wrong.Qeff <- Q.eff[!test.Qeff]
    
    message <- paste("The following QTL effects indicators:",
                     paste(wrong.Qeff, collapse = ", "),
                     "are not valid. Please use : 'cr', 'par', 'anc' or 'biall'")
    
    stop(message)
    
  }
  
  
  # 5. check the VCOV argument
  ############################
  
  
  # test if the asreml function is present for the compuation of the mixed models
  
  if(fct != "R2") {
    
    if (!(VCOV %in% c("h.err", "h.err.as", "cr.err", "pedigree", "ped_cr.err"))){
      
      stop(paste("The VCOV argument must be : 'h.err', 'h.err.as', 'cr.err',",
                 "'pedigree' or 'ped_cr.err'."))
      
    }
    
    
    if (VCOV != "h.err"){
      
      test <- requireNamespace(package = 'asreml', quietly = TRUE)
      
      if(!test){
        
        stop(paste("To use this type of VCOV, you must have access to the asreml",
                   "function from the asreml-R package."))
        
      }
      
    }
    
  }
  
  # Warning if the user want to use mixed model for
  
  if((fct == "forward") && (VCOV!="h.err")){
    
    message(paste("MESSAGE: The determination of a MQE model",
                  "using mixed models is technically possible but can take",
                  "a lot of time!"))
    
    text <- paste("Press [enter] to continue. If after you still want to stop",
                  "the process, Press Esc.")
    
    readline(prompt = text)
    
  }
  
  
  
  # 6. Consistency for parallelization.
  ####################################
  
  
  if ((n.cores > 1) && (VCOV != "h.err")){
    
    stop("Parallelization is only allowed for VCOV = 'h.err'.") 
    
    
  }
  
  
  # X. other checks according to specific type of functions
  #########################################################
  
  
  ### test the format of the QTL list introduce for backward elimination, R2 and
  # genetic effects estimation
  
  if(fct %in% c("R2", "QTLeffects")){
    
    if(is.null(QTL)){
      
      stop("No QTL position has been specified for the QTL argument.")
      
    }
    
    if(is.character(QTL)) {
      
      if(sum(!(QTL %in% mppData$map[, 1])) != 0){
        
        wrong.QTL <- QTL[!(QTL %in% mppData$map[, 1])]
        message <- paste("The following QTL positions:",
                         paste(wrong.QTL, collapse = ", "),
                         "are not present in the QTL profile (Qprof).")
        
        stop(message)
        
      }
      
    } else {
      
      stop("The QTL argument is not a list of character.")
      
    }
    
  }
  
  if(fct == "CIM"){
    
    if(is.null(cofactors)){
      
      stop("No cofactors have been provided.")
      
    }
    
  }
  
  if(fct == "proc"){ # test the validity of the path
    
    if(!file.exists(output.loc)){
      
      stop("The path specified in the argument output.loc is not valid.")
      
    }
    
  }
  
  }

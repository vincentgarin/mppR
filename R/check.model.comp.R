####################
# check.model.comp #
####################

# function to check the format of the data provided and the argument consistency
# before QTL model computation (perm, SIM, CIM, QTLEffects, etc.)

# parameters

# mppData object of class mppData.

# trait numeric or character indicator to specify the trait

# Q.eff Character expression specifying the type of QTL effect.

# VCOV Character specify the variance covariance structure.

# plot.gen.eff Logical value specifying if p-value of the single QTL effect
# must be stored.

# parallel logical indicating if fct must be run in parallel

# cluster cluster to run the function in parallel

# cofactors list of cofactors

# QTL list of QTLs

# ref.par character expression indicating a reference parent.

# sum_zero Logical expression indicating if the model should be computed using
# a sum to zero constraint.

# mppData.ts mppData object of the training set

# mppData.vs mppData object of the validation set

# fct specify which type of function to allow specific tests.


check.model.comp <- function(mppData = NULL, trait, Q.eff, VCOV,
                             plot.gen.eff = FALSE,
                             parallel = FALSE, cluster, cofactors = NULL,
                             QTL = NULL, ref.par = NULL, sum_zero = NULL,
                             mppData.ts = NULL, mppData.vs = NULL, fct = "XXX"){
  
  # 1. check mppData format
  #########################
  
  if(fct != "R2_pred"){
    
    if(is.null(mppData)){
      
      stop("no data have been provided for the mppData argument.")
      
    } else {
      
      if(!inherits(mppData, "mppData")) {
        
        stop("The data object provided (argument mppData) is not of class mppData.")
        
      }
      
    }
    
  }
  
  # 2. check trait
  ################
  
  check_trait(trait = trait, mppData = mppData)
  
  # 3. check Q.eff argument
  #########################
  
  if (!(Q.eff %in% c("cr", "par", "anc", "biall"))){
    
    stop("The Q.eff argument must be : 'cr', 'par', 'anc' or 'biall'")
    
  }
  
  # 4. check the VCOV argument
  ############################
  
  
  
  # test if the asreml function is present for the compuation of the mixed models
  
  if (fct != "R2"){
    
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
  
  # Warning if the user want to use mixed model for permutation
  
  if((fct=="perm") && (VCOV!="h.err")){
    
    message(paste("MESSAGE:",
                  "The determination of a significance threshold by permutation",
                  "using mixed models is technically possible but can take",
                  "a hugh amount of time! Due to limited computer power, we",
                  "advice to only use the linear model (VCOV = 'h.err')."))
    
    # ask user if he/she wants to continue
    
    text <- paste("Press [enter] to continue. If after you still want to stop",
                  "the process, Press Esc.")
    
    readline(prompt = text)
    
  }
  
  
  # 5. Consistency for parallelization.
  ####################################
  
  
  if (parallel){
    
    if(parallel & (!(VCOV == "h.err"))) {
      
      stop("Parallelization is only allowed for VCOV = 'h.err'.") 
      
    }
    
    if(!inherits(cluster, "cluster")){
      
      stop(paste("You must provide cluster objects to compute this function",
                 "in parallel. Use function makeCluster() from parallel pakage."))
      
    }
    
  }
  
  
  # 6. other checks according to specific type of functions
  #########################################################
  
  if((fct == "SIM")||(fct == "CIM")) {
    
    ### 5.1 test that if the user wants to fit a bi-allelic model, the
    # plot.gen.eff is not activated
    
    if((Q.eff == "biall") && plot.gen.eff) {
      
      stop(paste("The estimation of the decomposed QTL effect",
                 "(plot.gen.eff = TRUE) per cross or parents can not be performed for",
                 "the bi-allelic model"))
      
    }
    
    if(fct == "CIM"){
      
      if(is.null(cofactors)) {
        
        stop("No cofactors have been introduced.")
        
      }
      
    }
    
  }
  
  
  
  ### test the format of the QTL list introduce for backward elimination, R2 and
  # genetic effects estimation
  
  if(fct %in% c("back", "R2", "R2_pred", "QTLeffects")){
    
    if(is.null(QTL)){
      
      stop("No QTL position has been specified for the QTL argument")
      
    }
    
    if(is.character(QTL)) {
      
      if(sum(!(QTL %in% mppData$map[, 1])) != 0){
        
        wrong.QTL <- QTL[!(QTL %in% mppData$map[, 1])]
        message <- paste("The following QTL positions:",
                         paste(wrong.QTL, collapse = ", "),
                         "are not present in the QTL profile (Qprof).")
        
        stop(message)
        
      }
      
    } else { # the list of QTL is not character test if QTLlist format
      
      if (!inherits(QTL, "QTLlist")){
        
        stop(paste("The format of the argument QTL is not valid.",
                   "Use either a object obtained from the function QTL_select",
                   "or a character vector of markers or in between marker",
                   "positions indicator."))
        
      }
      
    }
    
    
    if (fct == "QTLeffects"){
    
      ### Test the compatibility of the reference parents
        
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
      
      ### Test correct configuration if sum_zero = TRUE
      
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
  
}
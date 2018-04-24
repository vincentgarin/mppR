####################
# check.model.comp #
####################

# function to check the format of the data provided and the argument consistency
# before QTL model computation (perm, SIM, CIM, QTLEffects, etc.)

# parameters

# mppData object of class mppData.

# Q.eff Character expression specifying the type of QTL effect.

# VCOV Character specify the variance covariance structure.

# par.clu parent clustering object.

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


check.model.comp <- function(mppData = NULL, Q.eff, VCOV, par.clu = NULL,
                             plot.gen.eff = FALSE, parallel = FALSE,
                             cluster, cofactors = NULL, QTL = NULL,
                             ref.par = NULL, sum_zero = NULL, mppData.ts = NULL,
                             mppData.vs = NULL, fct = "XXX"){
  
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
  
  
  # 2. check Q.eff argument
  #########################
  
  if (!(Q.eff %in% c("cr", "par", "anc", "biall"))){
    
    stop("The Q.eff argument must be : 'cr', 'par', 'anc' or 'biall'")
    
  }
  
  # 3. check the VCOV argument
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
  
  # 2. consistency between Q.eff and the type of mppData object
  #############################################################
  
  if(fct != "R2_pred"){
    
    if ((Q.eff=="cr" && mppData$biall) | (Q.eff=="par" && mppData$biall) |
        (Q.eff=="anc" && mppData$biall)){
      
      stop(paste("The mppData object is made for bi-allelic models and not for",
                 "cross, parental or ancestral models."))
      
    }
    
    if((Q.eff=="biall" & !mppData$biall)){
      
      stop(paste("The mppData object is made for cross, parental or ancestral",
                 "models and not for bi-allelic models."))
      
    }
    
  } else {
    
    # training set
    
    if ((Q.eff=="cr" && mppData.ts$biall) | (Q.eff=="par" && mppData.ts$biall) |
        (Q.eff=="anc" && mppData.ts$biall)){
      
      stop(paste("The mppData.ts object is made for bi-allelic models and not",
                 "for cross, parental or ancestral models."))
      
    }
    
    if((Q.eff=="biall" & !mppData.ts$biall)){
      
      stop(paste("The mppData.ts object is made for cross, parental or ancestral",
                 "models and not for bi-allelic models."))
      
    }
    
    # validation set
    
    if ((Q.eff=="cr" && mppData.vs$biall) | (Q.eff=="par" && mppData.vs$biall) |
        (Q.eff=="anc" && mppData.vs$biall)){
      
      stop(paste("the mppData.vs object is made for bi-allelic models and not",
                 "for cross, parental or ancestral models."))
      
    }
    
    if((Q.eff=="biall" & !mppData.vs$biall)){
      
      stop(paste("The mppData.vs object is made for cross, parental or",
                 "ancestral models and not for bi-allelic models."))
      
    }
    
  }
  
  
  
  # 3. Consistency for parallelization.
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
  
  
  
  # 4. if ancestral model, check the format of the par.clu object
  ##############################################################
  
  if(Q.eff == "anc"){
    
    if(is.null(par.clu)){
      
      stop(paste("You need to provide a parent clustering object",
                 "(argument par.clu) for the computation of an ancestral model."))
      
    }
    
    if(!is.integer(par.clu)){
      
      stop("The par.clu argument is not and integer matrix.")
      
    }
    
    if(fct != "R2_pred"){
      
      if(!identical(mppData$map[, 1], rownames(par.clu))){
        
        stop(paste("The list of markers and in between positions of the par.clu",
                   "object is not the same as the one in the mppData object map."))
        
      }
      
      
    } else {
      
      if(!identical(mppData.ts$map[, 1], rownames(par.clu))){
        
        stop(paste("The list of markers and in between positions of the par.clu",
                   "object is not the same as the one in the mppData object map."))
        
      }
      
    }
    
  }
  
  # 5. other checks according to specific type of functions
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
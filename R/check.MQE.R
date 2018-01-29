#############
# check.MQE #
#############

# function to check the format of the data provided to MQE functions.

# parameters

# mppData and mppData_bi:  objects of class mppData.

# Q.eff Character vector specifying the type of QTL effect.

# VCOV Character specify the variance covariance structure.

# par.clu parent clustering object.

# parallel logical indicating if fct must be run in parallel

# cluster cluster to run the function in parallel


check.MQE <- function(mppData = NULL, mppData_bi = NULL, Q.eff, VCOV,
                      par.clu = NULL, cofactors = NULL, cof.Qeff,
                      parallel = FALSE, cluster = NULL,
                      QTL = NULL, output.loc, fct = "XXX"){

  # 1. check mppData format
  #########################

  # mppData object for IBD

  if(is.null(mppData) & is.null(mppData_bi)){

    stop("No data have been provided for the mppData and mppData_bi arguments.")

  }

  if(is.null(mppData)) {

    stop("You must provide a mppData object in argument mppData.")

  }

  if(!inherits(mppData, "mppData")) {

    stop("The object provided in mppData is not of class mppData.")

  }

  # mppData object for IBS

  if(fct == "CIM"){

    test.bi <- (("biall" %in% Q.eff) || ("biall" %in% cof.Qeff))

  } else {test.bi <- "biall" %in% Q.eff }



  if(test.bi) {

    # test if the user gave a mppData_bi object made for this purpose

    if(is.null(mppData_bi)) {

      stop(paste("You must provide a mppData for bi-allelic data in",
                 "argument mppData_bi if you want to test for bi-allelic",
                 "QTL effects"))

    }

    if(!inherits(mppData_bi, "mppData")) {

      stop("The object provided in mppData_bi is not of class mppData.")

    }

    if(!mppData_bi$biall){

      stop("The mppData_bi object is not made for bi-allelic data.")

    }

  }

  # test if list of the two dataset are identical

  if((!is.null(mppData)) & (!is.null(mppData_bi))){

    if (!identical(as.character(mppData$map$mk.names),
                   as.character(mppData_bi$map$mk.names))) {



  stop(paste("The list of markers of the two mppData objects are not the same.",
             "For the computation of an MQE model marker list of the two",
             "mppData object must be strictly equivalent. This is potentially",
             "due to the addition of in between position in the formation of the",
             "mppData object. To avoid this problem, in the formation of the",
             "mppData object (mppData_form()), set the argument step with a large",
             "value (bigger than the largest in between marker gap) to avoid",
             "The introduction of in between positions."))

    }

    if (!identical(mppData$trait, mppData_bi$trait)) {

      stop(paste("The phenotypic values (trait) of the two mppData objects",
                 "are not the same."))

    }

  }



  # 2. check Q.eff argument
  #########################

  if((fct == "proc") | (fct == "forward")){

    if (length(Q.eff) <= 1) {
      stop("You must provide at least 2 type of QTL effect for the argument Q.eff.")
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


  # 3. if ancestral model, check the format of the par.clu object
  ##############################################################

  if ((fct == "proc") | (fct == "forward") | (fct == "CIM")) {

    test.anc <-(("anc" %in% Q.eff) || ("anc" %in% cof.Qeff))} else {

      test.anc <- "anc" %in% Q.eff

    }

  if(test.anc){

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

      stop(paste("The list of parents of the par.clu object (colnames(par.clu)) is
           not the same as the one of the mppData object."))

    }

  }

  # 4. consistency between Q.eff and the type of mppData object
  #############################################################

  test.Qeff <- c("cr", "anc", "par") %in% Q.eff

  if ((sum(test.Qeff) >= 1) && mppData$biall){

    stop("the mppData object is made for bi-allelic models and not for
         cross, parental or ancestral models.")

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

      if(!((exists("asreml")) && (is.function(asreml)))){

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


  if (parallel){

    if(parallel & (!(VCOV == "h.err"))) {

      stop("Parallelization is only allowed for VCOV = 'h.err'.")

    }

    if(!inherits(cluster, "cluster")){

      stop(paste("You must provide cluster objects to compute this function",
                 "in parallel. Use function makeCluster() from parallel pakage."))

    }

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

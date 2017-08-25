############
# MQE_proc #
############

#' Multi-QTL effect MPP analysis
#' 
#' This function aim at building multi-QTL models in which different QTL effects
#' (cross-specific, parental, ancestral or bi-allelic) can be assumed at
#' different loci. The possible QTL effect that the user want to allow must be
#' specified in \code{Q.eff}. Once the model has been determined it also computes
#' the QTL genetic effects per cross or parent and global and partial R squared
#' of the detected QTLs.
#' 
#' 
#' The procedure is the following:
#' 
#' \enumerate{
#' 
#' \item{forward regression to determine a multi-QTL model with different
#' possible assumptions for the QTL effect at different loci. The function use
#' \code{\link{MQE_forward}}.}
#' 
#' \item{Optional backward elimination (\code{backward = TRUE}) on the final
#' list of detected QTLs (\code{\link{MQE_BackElim}}).} 
#' 
#' \item{Estimation of the QTL genetic effects and R squared statistics
#' (\code{\link{MQE_genEffects}} and \code{\link{MQE_R2}}).}
#' 
#' \item{Optional plot (\code{plot.MQE = TRUE}) of the last CIM run of the
#' forward regression using the function \code{\link{MQE_plot}.}
#' 
#' }
#' 
#' }
#' 
#' \strong{WARNING!(1)} The computation of \code{MQE_proc()} function using mixed
#' models (all models with \code{VCOV} different than \code{"h.err"})
#' is technically possible but can be irrealistic
#' in practice due to a reduced computer power. Since a mixed model is computed at
#' each single position it can take a lot of time. From our estimation it can take
#' between 20 to 50 times more time than for linear models. If the number of
#' detected QTL is supposed to be small (until 5) it could still be feasible.
#' 
#' \strong{WARNING! (2)} The estimation of the random pedigree models
#' (\code{VCOV = "pedigree" and "ped_cr.err"}) can be unstable. Sometimes the
#' \code{asreml()} function fails to produce a results and returns the following
#' message: \strong{\code{GIV matrix not positive definite: Singular pivots}}.
#' So far we were not able to identify the reason of this problem and to
#' reproduce this error because it seems to happen randomly. From our
#' experience, trying to re-run the function one or two times should allow
#' to obtain a result.
#' 
#' @param pop.name \code{Character} name of the studied population.
#' Default = "MPP_MQE".
#' 
#' @param trait.name \code{Character} name of the studied trait.
#' Default = "trait1".
#' 
#' @param mppData An IBD object of class \code{mppData}
#' See \code{\link{mppData_form}} for details. Default = NULL.
#'
#' @param mppData_bi Required IBS object of class \code{mppData} if the user
#' wants to allow QTLs with a bi-allelic effect. \strong{The list of marker must
#' be strictly the same as the one of \code{mppData}.} Default = NULL.
#' 
#' @param Q.eff \code{Character} vector of possible QTL effects the user want to
#' test. Elements of Q.eff can be "cr", "par", "anc" or "biall". For details
#' look at \code{\link{mpp_SIM}}.
#' 
#' @param par.clu Required argument if the user wants to allow QTLs with an
#' ancestral effect. \code{interger matrix} representing the results of a parents
#' genotypes
#' clustering. The columns represent the parental lines and the rows
#' the different markers or in between positions. \strong{The columns names must
#' be the same as the parents list of the mppData object. The rownames must be
#' the same as the map marker list of the mppData object.} At a particular
#' position, parents with the same value are assumed to inherit from the same
#' ancestor. for more details, see \code{\link{USNAM_parClu}} and
#' \code{\link{parent_cluster}}. Default = NULL.
#'
#' @param VCOV \code{Character} expression defining the type of variance
#' covariance structure used: 1) "h.err" for an homogeneous variance residual term
#' (HRT) linear model; 2) "h.err.as" for a HRT model fitted by REML using
#' \code{ASReml-R}; 3) "cr.err" for a cross-specific variance residual terms
#' (CSRT) model; 4) "pedigree" for a random pedigree term and HRT model;
#' and 5) "ped_cr.err" for random pedigree and CSRT model.
#' For more details see \code{\link{mpp_SIM}}. Default = "h.err".
#' 
#' @param threshold \code{Numeric} value representing the -log10(p-value)
#' threshold above which a position can be considered as significant.
#' Significance threshold values can be obtained by permutation using
#' \code{\link{mpp_perm}} function. Default = 3.
#' 
#' @param window \code{Numeric} value in centi-Morgan representing the distance
#' on the left an right of a cofactor position where it is not included in the
#' model. Default value = 20.
#' 
#' @param backward \code{Logical} value. If \code{backward = TRUE},
#' the function performs
#' a backward elimination on the list of selected QTLs. Default = TRUE.
#' 
#' @param alpha.bk \code{Numeric} value indicating the significance level for
#' the backward elimination. Terms with p-values above this value will
#' iteratively be removed. Default = 0.05.
#' 
#' @param plot.MQE \code{Logical} value. If \code{plot.MQE = TRUE},
#' the function will make a plot of the last run of the MQE model
#' determination using function \code{\link{MQE_plot}}. Default = FALSE.
#' 
#' @param parallel \code{Logical} value specifying if the function should be
#' executed in parallel on multiple cores. To run function in parallel user must
#' provide cluster in the \code{cluster} argument. \strong{Parallelization is
#' only available for HRT (linear) models \code{VCOV = "h.err"}}.
#' Default = FALSE.
#' 
#' @param cluster Cluster object obtained with the function \code{makeCluster()}
#' from the \code{parallel} package. Default = NULL.
#' 
#' @param verbose \code{Logical} value indicating if the steps of MQE_proc should
#' be printed. It will not affect the printing of the other functions called by
#' \code{MQE_proc()}, especially the printing of \code{asreml()}.
#' Default = TRUE.
#' 
#' @param output.loc Path where a folder will be created to save the results.
#' By default the function uses the current working directory.
#' 
#' @return Return:
#'
#' \code{List} containing the following items:
#' 
#' \item{n.QTL}{Number of detected QTLs.}
#' 
#' \item{QTL}{\code{Data.frame} with QTL positions.}
#' 
#' \item{R2}{\code{list} containing R squared statistics of the QTL effects.
#' for details see \code{\link{QTL_R2}} and \code{\link{MQE_R2}} output
#' sections.}
#' 
#' \item{QTL.effects}{\code{List} of genetic effects per QTL. For details see
#' \code{\link{MQE_genEffects}} output section.}
#' 
#'
#' Some output files are also saved at the location specified
#' (\code{output.loc}):
#' 
#' \enumerate{
#' 
#' \item{A QTL report (QTL_REPORT.txt) with: 1) the number of detected QTLs;
#' 2) the global R squared statistics; 3) for each QTL, position information
#'  and estimated QTL genetic effect per cross or parents (for details see
#'  \code{\link{MQE_genEffects}}).}
#' 
#' \item{The list of QTLs (QTL.txt).}
#' 
#' \item{The QTL R squared statistics (QTL_R2.txt) (for details see
#' \code{\link{MQE_R2}} or \code{\link{QTL_R2}}).}
#' 
#' \item{General results of the QTL detection process: Number of QTL and
#' global adjusted and non-adjusted R squared statistics. (QTL_genResults.txt).}
#' 
#' \item{if \code{plot.MQE = TRUE}, a plot of the last QTL detection run profile
#' (plot_MQE.pdf) obtained with the function \code{\link{MQE_plot}}.}
#' 
#' 
#' }
#'
#' @author Vincent Garin
#' 
#' @seealso \code{\link{mppData_form}}, \code{\link{mpp_perm}},
#' \code{\link{mpp_SIM}},
#' \code{\link{MQE_forward}}, \code{\link{MQE_BackElim}},
#' \code{\link{MQE_genEffects}}, \code{\link{MQE_plot}}, \code{\link{MQE_R2}},
#' \code{\link{parent_cluster}}, \code{\link{QTL_R2}},
#' \code{\link{USNAM_parClu}}
#' 
#' @examples
#' 
#' data(USNAM_mppData)
#' data(USNAM_mppData_bi)
#' data(USNAM_parClu)
#' 
#' mppData <- USNAM_mppData
#' mppData_bi <- USNAM_mppData_bi
#' par.clu <- USNAM_parClu
#' 
#' 
#' # Equalize the list of markers of the two mppData objects and par.clu
#' 
#' com.mk.list <- intersect(mppData$map$mk.names, mppData_bi$map$mk.names)
#' 
#' mppData <- mppData_subset(mppData = mppData, mk.list = com.mk.list)
#' mppData_bi <- mppData_subset(mppData = mppData_bi, mk.list = com.mk.list)
#' par.clu <- par.clu[rownames(par.clu) %in% com.mk.list, ]
#' 
#' 
#' \dontrun{
#' 
#' # Specify a location where your results will be saved
#' my.loc <- "C:/.../..."
#' 
#' MQE <- MQE_proc(pop.name = "USNAM", trait.name = "ULA", mppData = mppData,
#'                 mppData_bi = mppData_bi, Q.eff = c("par", "anc", "biall"),
#'                 par.clu = par.clu, output.loc = my.loc)
#'                  
#' 
#' # Using parallel
#' 
#' library(parallel)
#' n.cores <- detectCores()
#' cluster <- makeCluster(n.cores-1)
#' 
#' MQE <- MQE_proc(pop.name = "USNAM", trait.name = "ULA", mppData = mppData,
#'                 mppData_bi = mppData_bi, Q.eff = c("par", "anc", "biall"),
#'                 par.clu = par.clu, parallel = TRUE, cluster = cluster,
#'                 output.loc = my.loc)
#'                  
#'                 
#' }
#' 
#' 
#' @export
#' 



MQE_proc <- function(pop.name = "MPP_MQE", trait.name = "trait1",
                     mppData = NULL, mppData_bi = NULL, Q.eff, par.clu = NULL,
                     VCOV = "h.err", threshold = 3, window = 20, backward = TRUE,
                     alpha.bk = 0.05, plot.MQE = FALSE, parallel = FALSE,
                     cluster = NULL, verbose = TRUE,
                     output.loc = getwd()) {
  
  # 1. check the format of the data
  #################################
  
  check.MQE(mppData = mppData, mppData_bi = mppData_bi, Q.eff = Q.eff,
            VCOV = VCOV, par.clu = par.clu, parallel = parallel,
            cluster = cluster, output.loc = output.loc, fct = "proc")
  
  
  # 2. Create a directory to store the results
  ############################################
  
  # create a directory to store the results of the QTL analysis
  
  end.char <- substr(output.loc, nchar(output.loc), nchar(output.loc))
  
  if(end.char == "/"){
    
    folder.loc <- paste0(output.loc, paste("MQE",pop.name, trait.name,
                                           VCOV, sep = "_"))
    
  } else {
    
    folder.loc <- paste0(output.loc, "/", paste("MQE",pop.name, trait.name,
                                                VCOV, sep = "_"))
    
  }
  
  dir.create(folder.loc)
  
  # 3. forward QTL detection
  ##########################
  
  if(verbose){
    
    cat("\n")
    cat("QTL detection")
    cat("\n")
    cat("\n")
    
  }
    
    QTL <- MQE_forward(mppData = mppData, mppData_bi = mppData_bi,
                       Q.eff = Q.eff, par.clu = par.clu, VCOV = VCOV,
                       threshold = threshold, window = window,
                       parallel = parallel, cluster = cluster,
                       verbose = verbose)
    
    if(backward){
      
      if(verbose){
        
        cat("\n")
        cat("Backward elimination")
        cat("\n")
        cat("\n")
        
      }
      
      QTL <- MQE_BackElim(mppData = mppData, mppData_bi = mppData_bi,
                          QTL = QTL[, 1], Q.eff = QTL[, 5], par.clu = par.clu,
                          VCOV = VCOV, alpha = alpha.bk)
      
      if (is.null(QTL)) { # test if QTL have been selected
        
        stop("No QTL position stayed in the model after the backward elimination.
         This is probably due to an error in the computation of the model
         in asreml() function.")
        
      }
      
    }
    
    # save the list of QTL
  
    file.QTL <- paste(folder.loc, "/", "QTL.txt", sep = "")
    
    write.table(QTL, file = file.QTL, quote = FALSE, row.names = FALSE, 
                sep = "\t")
  
  # 4. QTL effects and R2
  #######################
    
    if(verbose){
      
      cat("\n")
      cat("Computation QTL genetic effects and R2")
      cat("\n")
      cat("\n")
      
    }  
    
  QTL_effect <- MQE_genEffects(mppData = mppData, mppData_bi = mppData_bi,
                               QTL = QTL[, 1], Q.eff = QTL[, 5],
                               par.clu = par.clu, VCOV = VCOV)
  
  R2 <- MQE_R2(mppData = mppData, mppData_bi = mppData_bi, QTL = QTL[, 1],
               Q.eff = QTL[, 5], par.clu = par.clu, glb.only = FALSE)
  
  
  # save R2 results
  
  QTL.R2 <- data.frame(QTL[, 1:5], round(R2[[3]], 2), round(R2[[4]], 2),
                       round(R2[[5]], 2), round(R2[[6]], 2),
                       stringsAsFactors = FALSE)
  
  colnames(QTL.R2)[6:9] <- c("R2.diff", "adj.R2.diff", "R2.sg", "adj.R2.sg")
  
  write.table(QTL.R2, file = paste0(folder.loc, "/", "QTL_R2.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  # 5. Optional plot
  ##################
  
  if(plot.MQE){
    
    if(verbose){
      
      cat("\n")
      cat("Plot MQE last run profile")
      cat("\n")
      cat("\n")
      
    }
    
    CIM <- MQE_CIM(mppData = mppData, mppData_bi = mppData_bi, par.clu = par.clu,
                   VCOV = VCOV, cofactors = QTL[, 1], cof.Qeff = QTL[, 5],
                   chg.Qeff = TRUE, window = window, parallel = parallel,
                   cluster = cluster)
    
    main.plot <- paste("MQE", pop.name, trait.name, VCOV)
    
    pdf(paste0(folder.loc, "/", "plot_MQE.pdf"), height = 10, width = 16)
    
    print(MQE_plot(mppData = mppData, Qprof = CIM, QTL = QTL, window = window,
                   threshold = threshold, main = main.plot))
    
    dev.off()
    
  }
  
  # 6. results processing
  #######################
  
  if(verbose){
    
    cat("\n")
    cat("Results processing")
    cat("\n")
    cat("\n")
    
  }
  
  
  # QTL report
  
  QTL_report(out.file = paste0(folder.loc, "/", "QTL_REPORT.txt"),
             main = paste(pop.name, trait.name, VCOV), QTL.info = QTL,
             QTL.effects = QTL_effect, R2 = R2)
  
  # save general results
  
  gen.res <- c(dim(QTL)[1], round(R2[[1]][1], 2), round(R2[[2]][1], 2))
  names(gen.res) <- c("nb.QTL", "glb.R2", "glb.adj.R2")
  
  write.table(gen.res, file = paste0(folder.loc, "/", "QTL_genResults.txt"),
              quote = FALSE, sep = "\t", col.names = FALSE)
  
  
  
  # form the R object to be returned
  
  results <- list(n.QTL = dim(QTL)[1], QTL = QTL, R2 = R2,
                  QTL.effects = QTL_effect)
  
  return(results)
  
}
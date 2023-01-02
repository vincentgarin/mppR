#################
# print.summary #
#################

#' Print summary.mppData summary.QeffRes or summary.QR2Res object
#' 
#' @param x object of class \code{summary.mppData}, \code{summary.QeffRes} or,
#' \code{summary.QR2Res}
#' 
#' @param ... Ignored. 
#' 
#' @examples 
#' 
#' data(mppData)
#' sum.mppData <- summary(mppData)
#' print(sum.mppData)
#' 
#' data(mppData)
#' SIM <- mpp_SIM(mppData)
#' QTL <- QTL_select(SIM)
#' QTL.effects <- QTL_gen_effects(mppData = mppData, QTL = QTL, Q.eff = "cr")
#' sum.QeffRes <- summary(QTL.effects)
#' print(sum.QeffRes)
#' 
#' data(mppData)
#' SIM <- mpp_SIM(mppData)
#' QTL <- QTL_select(SIM)
#' Q_R2 <- QTL_R2(mppData = mppData, QTL = QTL, Q.eff = "cr")
#' summary(Q_R2)
#' 
#' @method print summary.mppData
#' @export
print.summary.mppData <- function(x, ...){
  
  cat("object of class 'mppData'", "\n \n")
  cat("\t Type of population: ", x$typePop, "\n \n")
  cat("\t No. Genotypes: ", x$Ngeno, "\n")
  print(x$par.per.cross, quote = FALSE)
  cat("\n")
  cat("\t Phenotype(s): ", paste(x$phenoName), "\n")
  cat("\t Percent phenotyped: ", paste(x$phenoPer), "\n \n")
  cat("\t Total marker: ", x$mkNb, "\n")
  cat("\t No. markers:", x$mkChr)
  
}

#' @method print summary.QeffRes
#' @export
print.summary.QeffRes <- function(x, ...){
  
  # write the title
  
  cat("QTL effects")
  cat("\n")
  cat(rep("*", 11), sep = "")
  cat("\n")
  cat("\n")
  
  # write the number of QTL and total R squared
  
  cat(paste("Number of QTL(s):", length(x$QTL.effects)))
  cat("\n")
  
  for (i in 1:length(x$QTL.effects)){
    
    cat("\n")
    cat("\n")
    cat(paste("QTL", x$Q_ind[i]))
    cat("\n")
    cat(rep("-",5), sep = "")
    cat("\n")
    cat("\n")
    
    # QTL position information (with or without CI information)
    
    print.data.frame(x$QTL.info[i, ], row.names = FALSE)
    
    
    cat("\n")
    cat("\n")
    cat("QTL effect per cross or parent:")
    cat("\n")
    cat("\n")
    
    print(x$QTL.effects[[i]])
    
    
  }
  
}

#' @method print summary.QR2Res
#' @export
print.summary.QR2Res <- function(x, ...){
  
  # write the title
  
  cat("QTL R2")
  cat("\n")
  cat(rep("*", 6), sep = "")
  cat("\n")
  cat("\n")
  
  if(x[[1]] == 'glb.only'){
    
    cat(paste("Global R2:", round(x$R2[[1]], 1)))
    cat("\n")
    cat("\n")
    cat(paste("Global adjusted R2:", round(x$R2[[2]], 1)))
    
  } else {
    
    cat(paste("Global R2:", round(x$R2[[1]], 1)))
    cat("\n")
    cat("\n")
    cat(paste("Global adjusted R2:", round(x$R2[[2]], 1)))
    cat("\n")
    cat("\n")
    
    # individual QTL results
    
    ind_Q_res <- data.frame(round(x$R2$part.adj.R2.sg, 1),
                            round(x$R2$part.adj.R2.diff, 1))
    
    colnames(ind_Q_res) <- c('part. adj. sg. R2', 'part. adj. diff. R2')
    
    print.data.frame(ind_Q_res)
    
    
  }
  
}

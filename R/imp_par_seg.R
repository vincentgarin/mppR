###############
# imp_par_seg #
###############

#' Impute parent marker scores
#' 
#' Impute the parent marker scores based on offpring segregation. The function
#' impute parent marker score in situation where one parent is missing and
#' the other is homozygous. In that case, it looks at the cross segregation.
#' If there is segratation, the missing parent get the allele marker score
#' different than other parent. If there is no segregation, the missing parent
#' gets the same score as the homozygous parent.
#' 
#' @param mppData  An object of class \code{mppData} formed with
#' \code{\link{create.mppData}}.
#' 
#' @return
#' 
#' parent marker score matrix with some missing marker imputed
#' 
#' @author Vincent Garin
#' 
#' @export


imp_par_seg <- function(mppData){
  
  # Restore data
  geno.off <- mppData$geno.off
  geno.par <- mppData$geno.par
  map <- mppData$map
  cross.ind <- mppData$cross.ind
  pheno <- mppData$pheno
  par.per.cross <- mppData$par.per.cross
  
  
  # Remove markers with genotyping error
  prob.mk <- QC_GenotypingError(mk.mat = rbind(geno.par, geno.off),
                                parallel = FALSE, cluster = 1)
  
  if(!is.null(prob.mk)){
    
    ind.prob <- which(colnames(geno.off) %in% prob.mk)
    geno.par <- geno.par[, -ind.prob]
    geno.off <- geno.off[, -ind.prob]
    
  }
  
  
  # initial percentage of NA
  miss_per <- sum(is.na(c(geno.par)))/(dim(geno.par)[1] * dim(geno.par)[2])
  
  allele.ref <- apply(X = rbind(geno.par, geno.off), MARGIN = 2,
                      FUN = allele.sc)[1:2, ]
  
  # detect the situation where it is impossible to get the reference alleles
  test_seg <- apply(allele.ref, 2, FUN = function(x) sum(is.na(x)))
  prob.mk <- (test_seg == 2)
  geno.par <- geno.par[, !prob.mk]
  geno.off <- geno.off[, !prob.mk]
  allele.ref <- allele.ref[, !prob.mk]
  
  geno.par.clu <- c()
  crosses <- unique(cross.ind)
  
  geno.par.imp <- c()
  
  for (k in seq_along(crosses)) {
    
    # select parents from the cross
    parents.line <- par.per.cross[k, 2:3]
    
    parent1.gen <- geno.par[which(rownames(geno.par) == parents.line[1]), ]
    parent2.gen <- geno.par[which(rownames(geno.par) == parents.line[2]), ]
    parents.gen <- rbind(parent1.gen, parent2.gen)
    
    off.gen <- subset(x = geno.off, subset = cross.ind == crosses[k], drop = FALSE)
    
    # Segregation pattern of the parents
    par.seg <- apply(X = parents.gen, MARGIN = 2, FUN = par.segretation)
    
    for (i in 1:length(par.seg)) {
      
      if(par.seg[i] == 2){
        
        # test if some info
        if(sum(!is.na(off.gen[, i])) != 0){
          
          # test if segregation
          if(length(unique(off.gen[, i])) == 1){
            
            # no segregation: P1 = P2
            parents.gen[is.na(parents.gen[, i]), i] <- parents.gen[!is.na(parents.gen[, i]), i]
            
            # segregation
          } else {
            
            # segregation: Pi -> other allele than (Pj)
            all_Pj <- parents.gen[!is.na(parents.gen[, i]), i]
            parents.gen[is.na(parents.gen[, i]), i] <- allele.ref[allele.ref[, i] != all_Pj, i]
            
          }
          
        }
      }
      
    }
    
    rownames(parents.gen) <- c(parents.line[1], parents.line[2])
    
    geno.par.imp <- rbind(geno.par.imp, parents.gen)
    
    
  }
  
  geno.par.imp <- geno.par.imp[rownames(geno.par), ]
  
  # percentage of NA after imputation
  miss_per2 <- sum(is.na(c(geno.par.imp)))/(dim(geno.par)[1] * dim(geno.par)[2])
  
  print(paste('missing percent before imputation:', paste0(miss_per, ';'),
              'missing percent after imputation', miss_per2))
  
  return(list(geno.par.imp = geno.par.imp, allele.ref = allele.ref))
  
}
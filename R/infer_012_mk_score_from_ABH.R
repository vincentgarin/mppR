###############################
# infer_012_mk_score_from_ABH #
###############################

#' Infer 012 marker score from ABH information
#' 
#' For each genotype marker scores, the function look at the most probable
#' parental origin (A, B or H), translate this information into an ATCG marker
#' score looking at the parent scores, and finally into 012 score. Such function
#' should be called if parents have non missing information. Marker imputation
#' of the parents can be done using \code{\link{impute_beagle}}.
#'
#' @param mppData An object of class \code{mppData}.
#' 
#' @param geno.par Parents genotype marker \code{matrix} with parent as row
#' and markers as columns. column names should be the marker identifiers
#' Row names should be the parents identifiers.
#' 
#' @return
#' 
#' Offspring genotype marker matrix in 012 format. 0: homozygous major allele,
#' 1: heterozygous, 2: homozygous minor allele.
#' 
#' @export
#' 

infer_012_mk_score_from_ABH <- function(mppData, geno.par){
  
  # restrict to the selection of markers present in mppData
  geno.par <- geno.par[, mppData$map$mk.names]
  
  par_nm <- rownames(geno.par)
  
  alleles <- mppData$allele.ref
  
  # convert the parents scores major and minor allele
  
  for(m in 1:dim(geno.par)[2]){
    
    geno.par[geno.par[, m] == alleles[1, m], m] <- 2
    geno.par[geno.par[, m] == alleles[2, m], m] <- 0
    geno.par[geno.par[, m] %in% alleles[3:4, m], m] <- 1
    
  }
  
  geno.par <- apply(geno.par, 2, as.numeric)
  rownames(geno.par) <- par_nm
  
  n.zigo <- mppData$n.zigo
  
  n_ind <- dim(mppData$pheno)[1]
  
  mk_mat <- matrix(NA, dim(mppData$geno.IBS)[1], dim(mppData$geno.IBS)[2])
  mk_count <- 1
  
  # Reconstruct a global ABH matrix (AA -> 1; BB -> 2; AB ->3)
  
  for(i in 1:length(mppData$geno.IBD$geno)){
    
    n_mk <- dim(mppData$geno.IBD$geno[[i]]$prob)[2]
    
    for(j in 1:n_mk){
      
      if(n.zigo == 3){
        
        pAA <- mppData$geno.IBD$geno[[i]]$prob[, j, 1]
        pAB <- mppData$geno.IBD$geno[[i]]$prob[, j, 2]
        pBB <- mppData$geno.IBD$geno[[i]]$prob[, j, 3]
        
      } else if (n.zigo == 2){
        
        pAA <- mppData$geno.IBD$geno[[i]]$prob[, j, 1]
        pAB <- rep(0, n_ind)
        pBB <- mppData$geno.IBD$geno[[i]]$prob[, j, 3]
        
      }
      
      mk_mat[, mk_count] <- apply(X = cbind(pAA, pBB, pAB), MARGIN = 1,
                                  FUN = function(x) which.max(x))
      
      mk_count <- mk_count + 1
      
    }
    
  }
  
  # iterate over the different crosses
  
  cross_ind <- mppData$cross.ind
  
  mk_mat_empt <- matrix(NA, dim(mk_mat)[1], dim(mk_mat)[2])
  
  for(c in 1:mppData$n.cr){
    
    parents.line <- mppData$par.per.cross[c, 2:3]
    
    parent1.gen <- geno.par[which(rownames(geno.par) == parents.line[1]), ]
    parent2.gen <- geno.par[which(rownames(geno.par) == parents.line[2]), ]
    parents.gen <- rbind(parent1.gen, parent2.gen)
    
    # select individual from cross c
    
    mk_mat_c <- mk_mat[cross_ind %in% mppData$par.per.cross[c, 1], ]
    
    for(k in 1:dim(mk_mat_c)[2]){
      
      mk_mat_c[mk_mat_c[, k] == 1, k] <- parents.gen[1, k]
      mk_mat_c[mk_mat_c[, k] == 2, k] <- parents.gen[2, k]
      mk_mat_c[mk_mat_c[, k] == 3, k] <- 0.5 * parents.gen[1, k] + 0.5 * parents.gen[2, k]
      
      # problem with the AB (H) situation
      
    }
    
    mk_mat_empt[cross_ind %in% mppData$par.per.cross[c, 1], ] <- mk_mat_c
    
  }
  
  rownames(mk_mat_empt) <- rownames(mppData$geno.off)
  colnames(mk_mat_empt) <- colnames(mppData$geno.off)
  
  return(mk_mat_empt)
  
}
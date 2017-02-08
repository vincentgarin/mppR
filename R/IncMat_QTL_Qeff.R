###################
# IncMat_QTL_Qeff #
###################

# Function to form QTL incidence matrix with a list of QTL that has different
# type of QTL effects


IncMat_QTL_Qeff <- function(x, QTL, mppData, Q.eff, par.clu, ref.all.most) {
  
  
  # 2. for parental or ancestral model change the order
  
    if(Q.eff == "par"){
      
      # determine the connected parts
      
      con.part <- design_connectedness(par.per.cross = mppData$par.per.cross,
                                       plot.des = FALSE)
      
      allele_order <- c()
      allele_ref <- c()
      
      for(i in seq_along(con.part)){
        
        con.part_i <- con.part[[i]]
        
        # subset the par.per.cross object
        
        index <- apply(X = mppData$par.per.cross[, c(2,3)], MARGIN = 1,
                       FUN = function(x, ref) sum(x %in% ref) > 0,
                       ref = con.part_i)
        
        par.per.cross_i <- mppData$par.per.cross[index, , drop = FALSE]
        
        # susbet the cross indicator according to the cross retained
        
        cross.ind_i <- mppData$cross.ind[mppData$cross.ind %in% par.per.cross_i[, 1]]
        
        allele_ord_i <- most.used.allele(par.per.cross_i, cross.ind_i,
                                         most.used = !ref.all.most)
        
        allele_order <- c(allele_order, allele_ord_i)
        allele_ref <- c(allele_ref, allele_ord_i[length(allele_ord_i)])
        
      }
      
      # order the parental alleles in QTL incidence matrices
      
      QTL <- QTL[, allele_order]
      
    } else if (Q.eff == "anc"){
      
      # determine the connected parts
      
      par.clu_x <- par.clu[x, ]
      par.clu_x <- paste0("A.allele", par.clu_x)
      names(par.clu_x) <- mppData$parents
      
      all.p1 <- par.clu_x[mppData$par.per.cross[, 2]]
      all.p2 <- par.clu_x[mppData$par.per.cross[, 3]]
      
      par.per.cross_x <- cbind(mppData$par.per.cross[, 1], all.p1, all.p2)
      
      con.part <- design_connectedness(par.per.cross = par.per.cross_x,
                                       plot.des = FALSE)
      
      # within connected part determine the most (least) used allele
      
      allele_order <- c()
      allele_ref <- c()
      
      for(i in seq_along(con.part)){
        
        con.part_i <- con.part[[i]]
        
        # subset the par.per.cross object
        
        index <- apply(X = par.per.cross_x[, c(2,3)], MARGIN = 1,
                       FUN = function(x, ref) sum(x %in% ref) > 0,
                       ref = con.part_i)
        
        par.per.cross_i <- par.per.cross_x[index, , drop = FALSE]
        
        # susbet the cross indicator according to the cross retained
        
        cross.ind_i <- mppData$cross.ind[mppData$cross.ind %in% par.per.cross_i[, 1]]
        
        allele_ord_i <- most.used.allele(par.per.cross_i, cross.ind_i,
                                         most.used = !ref.all.most)
        
        allele_order <- c(allele_order, allele_ord_i)
        allele_ref <- c(allele_ref, allele_ord_i[length(allele_ord_i)])
        
      }
      
      # order the parental alleles in QTL incidence matrices
      
      QTL <- QTL[, allele_order]
      
    }
  
  return(QTL)
  
}
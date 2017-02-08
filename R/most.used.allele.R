####################
# most.used.allele #
####################

# function to determine which parental (ancestral) allele is the most used
# according to first the number of cross where it segregate and second the
# number of genotype where it could be present

# Return the odered list of parental (ancestral) most used alleles 

most.used.allele <- function(par.per.cross, cross.ind, most.used){
  
  par.list <- unique(c(par.per.cross[, 2], par.per.cross[, 3]))
  n.par <- length(par.list)
  n.cr <- dim(par.per.cross)[1]
  
  # 1. determine in how many cross the allele is used ?
  
  n.cr.p <- table(c(par.per.cross[, 2], par.per.cross[, 3]))
  n.cr.p <- n.cr.p[par.list]
  
  ######### problem there
  
  n.cr.p <- as.data.frame(n.cr.p)
  
  if(dim(n.cr.p)[2] == 1) {colnames(n.cr.p)[1] <- "n.cr.p"
  
  } else if(dim(n.cr.p)[2] == 2) {colnames(n.cr.p)[2] <- "n.cr.p"}
  
  ####### stop
  
  # explanation: In linux (Other versuion of R than 3.3.2) when the table is
  # subseted it becomes a vector (1 column). On window it sort of keep the
  # 'table' structure with two columns. 
  
  # 2. determine in how many genotypes an allele potentially segregates ?
  
  n.gen <- table(factor(cross.ind, levels = unique(cross.ind)))
  n.gen <- n.gen[par.per.cross[, 1]]
  
  
  n.off.p <- rep(0, n.par)
  
  for (i in 1:n.par){
    
    par_i <- par.list[i]
    
    n.off.pi <- c()
    
    for(j in 1:n.cr){
      
      if(par_i %in% par.per.cross[j, 2:3]) {
        
        n.off.pi <- c(n.off.pi, n.gen[j])
        
      }
      
    }
    
    n.off.p[i] <- sum(n.off.pi)
    
  }
  
  # 3. Order the parent list
  
  data.par <- data.frame(par.list, n.cr.p, n.off.p, stringsAsFactors = FALSE)
  data.par <- data.par[order(data.par[, 'n.cr.p'], data.par[, 'n.off.p'],
                             decreasing = most.used), ]
  parents.ord <- data.par[, 1]
  
  return(parents.ord)
  
}
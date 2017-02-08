##################
# mppData_subset #
##################

#' Subset mppData object
#' 
#' Subsets \code{mppData} objects by markers (\code{mk.list}) or by genotypes
#' (\code{gen.list}).
#' 
#' @param mppData An object of class \code{mppData}.
#' See \code{\link{mppData_form}} for details.
#' 
#' @param mk.list \code{Character} vector of marker or in betwee position names,
#' \code{numeric} vector of marker or inbetween positions indicator,
#' or \code{logical} vector. Default = NULL.
#' 
#' @param gen.list \code{Character} vector of genotypes names, \code{numeric}
#' vector of genotypes, or \code{logical} vector. Default = NULL.
#' 
#' @return Return:
#' 
#' \item{mppData_red}{Subseted \code{mppData} object with selected markers
#' or selected genotypes.}
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{mppData_form}}, \code{\link{mppData_chgPheno}}
#' 
#' @examples
#' 
#' ### marker subset
#' 
#' data(USNAM_mppData)
#' 
#' # Random selection of markers
#' mk.list <-  sample(USNAM_mppData$map[, 1], 50)
#' mppData_sub <- mppData_subset(mppData = USNAM_mppData, mk.list = mk.list)
#' 
#' # Selection marker of chromosome 1
#' mk.list <-  (USNAM_mppData$map[, 2] == 1)
#' mppData_sub <- mppData_subset(mppData = USNAM_mppData, mk.list = mk.list)
#' 
#' ### genotype subset
#' 
#' # Random selection of genotypes
#' gen.list <-  sample(USNAM_mppData$geno.id, 200)
#' mppData_sub <- mppData_subset(mppData = USNAM_mppData, gen.list = gen.list)
#' 
#' # Selection of genotype from cross 2 and 5
#' crosses <- unique(USNAM_mppData$cross.ind)
#' gen.list <-  USNAM_mppData$geno.id[USNAM_mppData$cross.ind %in% crosses[c(2, 5)]]
#' mppData_sub <- mppData_subset(mppData = USNAM_mppData, gen.list = gen.list)
#' 
#' @export
#'  


mppData_subset <- function(mppData, mk.list = NULL, gen.list = NULL) {
  
  # 1. Check data.format and arguments
  ####################################
  
  stopifnot(inherits(mppData, "mppData"))
  
  # user must at least specify one argument (mk.list or gen.list)
  
  if(is.null(mk.list) && is.null(gen.list)){
    
    stop("None of the argument specifying what should be substituted has been
         specified. Please mk.list or gen.list")
    
  }
  
  # user must only specify one argument (mk.list or gen.list)
  
  if(!is.null(mk.list) && !is.null(gen.list)){
    
    stop("You can only specify one argument at a time. Either subsetting by
         marker (mk.list) or by genotype (gen.list).")
    
  }
  
  # test the format of the marker and genotype list given by the user
  
  
  if (!is.null(mk.list)){
    
    # test format marker list
    
    if (!(is.character(mk.list) || is.logical(mk.list) || is.numeric(mk.list))) {
      
      stop("The marker list must be a character, logical or numeric vector
           specifying the list of marker to subset.")    
      
    }
    
    if(is.logical(mk.list)){
      
      if(length(mk.list) != dim(mppData$map)[1])
        
        stop("The mk.list does not have the same length as the map")
      
    }
    
    # Convert different type of marker list into a character vector.
    
    if (is.logical(mk.list) || is.numeric(mk.list)) {
      
      mk.list <- mppData$map[mk.list, 1]
      
    } 
    
    
  } else if (!is.null(gen.list)){
    
    # test format marker list
    
  if (!(is.character(gen.list) || is.logical(gen.list) || is.numeric(gen.list))) {
      
      stop("The genotypes list (gen.list) must be a character, logical or numeric
            vector specifying the list of genotypes to subset.")    
      
    }
    
    if(is.logical(gen.list)){
      
      if(length(gen.list) != length(mppData$geno.id))
        
    stop("The genotype list (gen.list) does not have the same length as the 
         genotypes list")
      
    }
    
    # Convert different type of marker list into a character vector. Then in a
    # logical list to keep the same order.
    
    if (is.logical(gen.list) || is.numeric(gen.list)) {
      
      gen.list <- mppData$geno.id[gen.list]
      
      geno.ind <- mppData$geno.id %in% gen.list 
      
    } else if(is.character(gen.list)) {
      
    geno.ind <- mppData$geno.id %in% gen.list
        
    }
    
  }
  
  
  # 2. cross, parental, ancestral mppData object
  ##############################################
  
  if (!mppData$biall) {
    
  
      ### 2.1 substitute by marker
    
    if (!is.null(mk.list)){
      
    
      # modify the genotype marker matrix
      
      for (i in 1:length(mppData$geno$geno)) {
        
        # get the list of marker of the considered chromosome
        
        chr.mk.names <- attr(mppData$geno$geno[[i]]$prob, "dimnames")[[2]]
        
        # indicate which marker of chromosome i is in the reduced list
        
        indicator <- chr.mk.names %in% mk.list
        
        mppData$geno$geno[[i]]$prob <- mppData$geno$geno[[i]]$prob[, indicator, ]
        
        
      }
      
      # modify the marker geno.par argument (if present)
      
      
      if (!is.null(mppData$geno.par)) {
        
        mppData$geno.par <- mppData$geno.par[(mppData$geno.par[, 1] %in% mk.list), ]
        
        
      }
      
      # modify the map
      
      mppData$map <- mppData$map[(mppData$map[, 1] %in% mk.list), ]
      
      # recalculate the position indicators
      
      mppData$map[, 3] <- sequence(table(mppData$map[, 2]))
      
      
      return(mppData)
      
      
      ### 2.2 substitute by genotype      

    } else if (!is.null(gen.list)){

      # modify genotype matrix
      
      
      for (i in 1:length(mppData$geno$geno)) {
        
        
        mppData$geno$geno[[i]]$prob <- mppData$geno$geno[[i]]$prob[geno.ind, , ]
        
      }
      
      # subset geno.id
      
      mppData$geno.id <- mppData$geno.id[geno.ind]
      
      # subset trait
      
      mppData$trait <- subset(x = mppData$trait, subset = geno.ind, drop = FALSE)
      
      # subset cross indicator
      
      mppData$cross.ind <- mppData$cross.ind[geno.ind]
      
      # subset ped.mat
      
      ped.temp <- as.matrix(mppData$ped.mat)
      
      ped.mat.found <- ped.temp[ped.temp[, 1] == "founder", ]
      ped.mat.off <- ped.temp[ped.temp[, 1] == "offspring", ]
      ped.mat.off <- ped.mat.off[geno.ind, ]
      
      found.sub <- unique(c(ped.mat.off[, 3], ped.mat.off[, 4]))
      ped.mat.found <- ped.mat.found[ped.mat.found[, 2] %in% found.sub, ,
                                     drop = FALSE]
      
      pedigree.new <- rbind(ped.mat.found, ped.mat.off)
      
      mppData$ped.mat <- data.frame(pedigree.new, stringsAsFactors = FALSE)
      
      # review the par.per.cross, parents, n.cr and n.par objects
      
      list.cr <- unique(mppData$cross.ind)
      
      ppc <- mppData$par.per.cross
      ppc <- ppc[ppc[, 1] %in% list.cr, ]
      mppData$parents <- union(ppc[, 2], ppc[, 3])
      mppData$n.cr <- length(list.cr)
      mppData$n.par <- length(mppData$parents)
      mppData$par.per.cross <- ppc
      
      return(mppData)

    }



    # 3. Bi-allelic mppData object
    ##############################


  } else {

    ### 2.1 substitute by marker

    if (!is.null(mk.list)){


      # modify the genotype marker matrix

      mppData$geno <- mppData$geno[, mppData$map[, 1] %in% mk.list]
      
      # modify the reference alleles scores
      
      mppData$allele.ref <- mppData$allele.ref[, mppData$map[, 1] %in% mk.list]
      
      # modify the marker geno.par argument (if present)


      if (!is.null(mppData$geno.par)) {

      mppData$geno.par <- mppData$geno.par[(mppData$geno.par[, 1] %in% mk.list), ]


      }

      # modify the map

      mppData$map <- mppData$map[(mppData$map[, 1] %in% mk.list), ]

      # recalculate the position indicators

      mppData$map[, 3] <- sequence(table(mppData$map[, 2]))


      return(mppData)


      ### 2.2 substitute by genotype

    } else if (!is.null(gen.list)){

      # modify genotype matrix
      
      mppData$geno <- mppData$geno[geno.ind, ]

      # subset geno.id

      mppData$geno.id <- mppData$geno.id[geno.ind]

      # subset trait

      mppData$trait <- subset(x = mppData$trait, subset = geno.ind, drop = FALSE)

      # subset cross indicator

      mppData$cross.ind <- mppData$cross.ind[geno.ind]

      # subset ped.mat

      ped.mat.found <- mppData$ped.mat[mppData$ped.mat[, 1] == "founder", ]
      ped.mat.off <- mppData$ped.mat[mppData$ped.mat[, 1] == "offspring", ]
      ped.mat.off <- ped.mat.off[geno.ind, ]
      mppData$ped.mat <- rbind(ped.mat.found, ped.mat.off)
      
      # review the par.per.cross, parents, n.cr and n.par objects
      
      list.cr <- unique(mppData$cross.ind)
      
      ppc <- mppData$par.per.cross
      ppc <- ppc[ppc[, 1] %in% list.cr, ]
      mppData$parents <- union(ppc[, 2], ppc[, 3])
      mppData$n.cr <- length(list.cr)
      mppData$n.par <- length(mppData$parents)
      mppData$par.per.cross <- ppc
      

      return(mppData)

    }


   }
  
} # end function
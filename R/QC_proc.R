###########
# QC_proc #
###########

#' Quality control procedure
#'
#' Perform different operation of quality control on the marker matrix and
#' return marker matrix, map and phenotype data in format to form a
#' \code{mppData} object.
#' 
#' The different operation of the quality control arethe following:
#' 
#' \enumerate{
#' 
#' \item{Remove markers with genotyping error (more than 2 possible alleles)
#' \code{\link{QC_GenotypingError}}.}
#' 
#' \item{Remove markers that are completely monomorphic or missing
#' \code{\link{QC_MAF}}.}
#' 
#' \item{Remove marker with a population minor allele frequency (MAF) below the
#' value provided in \code{MAF.pop.lim}.}
#' 
#' \item{Remove marker with a missing rate higher than \code{mk.miss} at the
#' population level \code{\link{QC_missing}}.}
#' 
#' \item{If the map contains some markers at the same position, the function
#' keep only the most polymorphic position.}
#' 
#' \item{Remove crosses with less than \code{n.lim} observations
#' \code{\link{QC_minCrSize}}. \code{n.lim} can not take a value below 15.}
#' 
#' \item{Remove genotypes with a missing rate higher than \code{gen.miss}
#' \code{\link{QC_missing}}.}
#' 
#' \item{Determine markers having a problematic MAF within cross
#' \code{\link{QC_MAF}} and \code{\link{QC_tagMAFCr}}. The critical within
#' cross MAF can be specified by the user via the argument \code{MAF.cr.lim}.
#' By default, the critical within cross MAF are defined by the following
#' function of the cross-size (n.cr): MAF(n.cr) = (5/n.cr) + 0.05.
#' This means that for small cross sizes (n.cr = \code{n.lim} = 15),
#' the crosses must have at least a bit more than five genotypes that segreate.
#' When the number of genotypes per cross increases, then the within cross MAF
#' tend to 0.05. If the within cross MAF is below the limit in at least one
#' cross, then marker scores of the problematic cross are either put as missing
#' (\code{MAF.cr.miss = TRUE} default) or the whole marker is discared
#' (\code{MAF.cr.miss = FALSE}).}
#' 
#' \item{If \code{rem.NA.cr = TRUE}, remove the markers that are
#' completely missing in at least one cross. Default = FALSE}
#' 
#' \item{If \code{ABH = TRUE}, convert offspring genotype data into ABH format
#' (\code{\link{cross_ABH}}).
#' If parents have heterozygous marker scores, user must specify it
#' using \code{het.par = TRUE} (\code{\link{cross_ABH_het}}).}
#' 
#' }
#' 
#' 
#' @param geno.off \code{Character} marker score \code{matrix} of the offspring
#' with genotypes as row and markers as column.
#' \strong{Rows names must be the offspring genotypes identifiers similar to
#' the one used in \code{trait}. The columns names must be the marker names
#' similar to the one used in \code{map}. Marker scores must be coded using one
#' letter per allele. For example, AA, CC, GG, TT, AC, AG, AT, CA, CG, CT,
#' GA, GC, GT, TA, TC, TG. Missing values must be coded \code{NA}.}
#' 
#' @param geno.par \code{Character} marker score \code{matrix} of the parents
#' with genotypes as row and markers as column.
#' \strong{Rows names must be the parents genotypes identifiers similar to
#' the one used in \code{par.per.cros}. The columns names must be the marker
#' names similar to the one used in \code{map}. Marker scores must be coded
#' using one letter per allele. For example, AA, CC, GG, TT, AC, AG, AT, CA,
#' CG, CT, GA, GC, GT, TA, TC, TG. Missing values must be coded \code{NA}.}
#' 
#' @param map Three columns \code{data.frame} with: 1) marker or in between
#' position identifiers; 2) chromosome; 3) positions in centi-Morgan.\strong{
#' The marker identifiers must be identical to the column names of the maker
#' matrices (\code{geno.off} and \code{geno.par}).}
#' 
#' @param trait two columns \code{data.frame} with : 1) \code{character}
#' genotypes identifiers; 2) \code{numeric} trait values. \strong{The genotypes
#' identifiers must be identical to the rownames of the offspring marker matrix
#' (\code{geno.off}).}
#' 
#' @param cross.ind \code{Character} vector indicating to which cross each
#' genotype belong.
#' 
#' @param par.per.cross Three columns \code{Character matrix} specifying :
#' 1) the cross indicators; 2) the parents 1 identifiers
#' of the crosses; 3) the parents 2 identifiers of the crosses.
#' \strong{The list of cross must contain the same cross indicators as in
#' \code{cross.ind} and they must appear in the same order.
#' The list of parent identifiers must be the same to the rownames of
#' the argument \code{geno.par}}.
#' 
#' @param subcross.ind Optional \code{character} vector specifying to which
#' sub-cross each genotype belong. \strong{This argument is only necessary if
#' the user want to make ABH assignement per sub-cross. In that case argument
#' \code{par.per.subcross must also be provided.}}. Default = NULL.
#' 
#' @param par.per.subcross Optional three columns \code{Character matrix}
#' specifying : 1) the sub-cross indicators; 2) the parents 1 identifiers
#' of the sub-crosses; 3) the parents 2 identifiers of the sub-crosses.
#' \strong{This argument is only necessary if the user want to make ABH
#' assignement per sub-cross. In that case, the list of cross must contain the
#' same cross indicators as in \code{subcross.ind} and they must appear in the
#' same order. The list of parent identifiers must be the same to the rownames of
#' the argument \code{geno.par}}. Default = NULL.
#' 
#' @param n.lim \code{Numeric} value specifying the minimum cross size.
#' Default = 15. \strong{it is not allowed to use cross with less than 15
#' observation (\code{n.lim < 15}).}
#' 
#' @param MAF.pop.lim \code{Numeric} value specifying the minimum minor allele
#' frequency for a marker at the population level. Default = 0.05.
#' 
#' @param MAF.cr.lim \code{Numeric vector} specifying the critical within cross
#' MAF. Default = NULL. Marker with a problematic segregation rate in at least
#' one cross is either set as missing within the problematic cross
#' (\code{MAF.cr.miss = TRUE}), or remove from the marker matrix
#' (\code{MAF.cr.miss = FALSE}). Default = NULL.
#' 
#' @param MAF.cr.miss Logical value specifying if maker with a too low
#' segregation rate within cross should be put as missing or discarded. If
#' \code{MAF.cr.miss = TRUE}, marker scores with a problematic segregation rate
#' in at least one cross will be put as missing in the problematic cross(es).
#' If \code{MAF.cr.miss = FALSE}, the whole marker will be discarded.
#' Default = TRUE.
#' 
#' @param rem.NA.cr \code{Logical} value specifying if the marker that are
#' completely missing in at least one cross should be remored. Default = FALSE.
#' 
#' @param mk.miss \code{Numerical} value comprised between 0 and 1 indicating
#' the missingness rate at the population level above which a marker will be
#' removed. Default = 0.1.
#' 
#' @param gen.miss \code{Numerical} value comprised between 0 and 1 indicating
#' the missingness rate above which a genotype will be removed. Default = 0.25.
#' 
#' @param ABH \code{Logical} value specifying if the marker matrix data should
#' be converted into ABH format within cross. \strong{If user wants to convert
#' the markers to the ABH format and parent have heterozygous marker score.
#' Then the argument \code{het.par} must be set to TRUE}.
#' If (\code{ABH = FALSE}), the marker data will be let in the same format as
#' the one used in \code{geno.off}. Default = TRUE.
#' 
#' @param het.par \code{Logical} value. if \code{het.par = TRUE}, the function
#' will use the offspring segregation to try to infer the allele that was
#' transmitted by the heterozygous parent at a particular locus in order to
#' make the ABH conversion. Default = FALSE.
#' 
#' @param parallel \code{Logical} value specifying if the function should be
#' executed in parallel on multiple cores. To run function in parallel user must
#' provide cluster in the \code{cluster} argument. Default = FALSE.
#' 
#' @param cluster Cluster object obtained with the function \code{makeCluster()}
#' from the parallel package. Default = NULL.
#' 
#' @return 
#' 
#' list containing the following objects
#' 
#' \item{geno.off}{Subseted offspring marker martrix. if \code{ABH = TRUE}
#' markers scores transformed into ABH format.}
#' 
#' \item{geno.par}{Subseted parents marker martrix.}
#' 
#' \item{geno.par.clu}{Parent marker matrix with only problematic markers,
#' marker fully monomorphic or fully missing removed. This element can be
#' used for the argument \code{marker.data} of the function
#' \code{\link{parent_cluster}}.}
#' 
#' \item{map.par.clu}{Genetic map corresponding to the list of marker of the
#' \code{geno.par.clu} object. This map can be used for the argument
#' \code{haplo.map} in function \code{\link{parent_cluster}}.}
#' 
#' \item{cross.ind}{Corresponding cross indicator vector.}
#' 
#' \item{par.per.cross}{Matrix specifying the crosses and for each cross the
#' parent 1 of the cross and the parent 2.} 
#' 
#' \item{trait}{Corresponding phenotypic values.}
#' 
#' \item{map}{Corresponding genotypic map.}
#' 
#' \item{rem.mk}{Vector of markers that have been removed.}
#' 
#' @author Vincent Garin
#' 
#' @seealso
#' 
#' \code{\link{cross_ABH}}, \code{\link{cross_ABH_het}}
#' \code{\link{QC_GenotypingError}}, \code{\link{QC_MAF}},
#' \code{\link{QC_minCrSize}}, \code{\link{QC_missing}},
#' \code{\link{QC_tagMAFCr}}, \code{\link{parent_cluster}}
#' 
#' @examples
#' 
#' data(USNAM_geno)
#' data(USNAM_pheno)
#' data(USNAM_map)
#' 
#' cross.ind <- substr(USNAM_pheno[, 1], 1, 4)
#' geno.off <- USNAM_geno[7:506, ]
#' geno.par <- USNAM_geno[1:6, ]
#' 
#' map <- USNAM_map
#' trait <- USNAM_pheno
#' par.per.cross <- cbind(unique(cross.ind), rep("B73", 5),
#'                        rownames(geno.par)[2:6])
#' 
#' # Quality control procedure
#' 
#' data <- QC_proc(geno.off = geno.off, geno.par = geno.par, map = map,
#'                 trait = trait, cross.ind = cross.ind,
#'                 par.per.cross = par.per.cross, n.lim = 15,
#'                 MAF.pop.lim = 0.05, mk.miss = 0.1,
#'                 gen.miss = 0.25, ABH = TRUE, het.par = TRUE,
#'                 parallel = FALSE)
#' 
#' # QC_proc() outputs can be directly used to form the mppData object
#' # (gathering of all data necessary for the QTL analysis) and to form
#' # the parent clustering object necessary for the ancestral model
#' 
#' # formation of mppData object
#' 
#' \dontrun{
#' 
#' mppData <- mppData_form(geno = data$geno.off, geno.par = data$geno.par,
#'                         biall = FALSE, type = "F", nb.gen = 6,
#'                         type.mating = "selfing", map = data$map,
#'                         trait = data$trait, cross.ind = data$cross.ind,
#'                         par.per.cross = data$par.per.cross,
#'                         step = 5, dir = getwd())
#' 
#' # Parent clustering
#' 
#' library(clusthaplo)
#'   
#' cluster <- parent_cluster(haplo.map = data$map.par.clu,
#'                           consensus.map = mppData$map[, -3],
#'                           marker.data = t(data$geno.par.clu), na.strings = NA,
#'                           step.size = 1000, window = 25)
#' 
#' par.clu <- cluster[[1]]
#' 
#' # then these data can be used for QTL detection. Let us for example compute
#' # a QTL detection based on the ancestral model
#' 
#' QTL <- mpp_proc(pop.name = "USNAM", trait.name = "ULA", mppData = mppData,
#'                 Q.eff = "anc", par.clu = par.clu, output.loc = getwd())
#' 
#' }
#'
#' @export
#' 


QC_proc <- function(geno.off, geno.par, map, trait, cross.ind, par.per.cross,
                    subcross.ind = NULL, par.per.subcross = NULL,
                    n.lim = 15, MAF.pop.lim = 0.05, MAF.cr.lim = NULL,
                    MAF.cr.miss = TRUE, rem.NA.cr = FALSE, mk.miss = 0.1,
                    gen.miss = 0.25, ABH = TRUE, het.par = FALSE,
                    parallel = FALSE, cluster = NULL){
  
  
  # 1. check the format of the data
  #################################
  
  check_QC(geno.off = geno.off, geno.par = geno.par, map = map, trait = trait,
           cross.ind = cross.ind, par.per.cross = par.per.cross,
           subcross.ind = subcross.ind, par.per.subcross = par.per.subcross,
           n.lim = n.lim, MAF.cr.lim = MAF.cr.lim, ABH = ABH, het.par = het.par,
           parallel = parallel, cluster = cluster)
  
  
  # 2. Remove markers with genotyping error
  #########################################
  
  prob.mk.list <- c()
  
  prob.mk <- QC_GenotypingError(mk.mat = rbind(geno.par, geno.off),
                                parallel = parallel, cluster = cluster)
  
  if(!is.null(prob.mk)){
    
    prob.mk.list <- c(prob.mk.list, prob.mk)
    
    ind.prob <- which(colnames(geno.off) %in% prob.mk)
    geno.par <- geno.par[, -ind.prob]
    geno.off <- geno.off[, -ind.prob]
    
  }
  
  # 3. remove completely monomorphic or missing markers
  ######################################################
  
  ### 3.1 monomorphic or missing in the parents
  
  parent.MAF <- QC_MAF(mk.mat = geno.par, parallel = parallel,
                       cluster = cluster)
  
  mono <- which(parent.MAF == 0)
  miss <- which(is.na(parent.MAF))
  
  prob.mk.id <- c(mono, miss)
  
  if(length(prob.mk.id) > 0){
    
    prob.mk.list <- c(prob.mk.list, colnames(geno.par)[prob.mk.id])
    
    geno.par <- geno.par[, -prob.mk.id]
    geno.off <- geno.off[, -prob.mk.id]
    
  }
  
  ### 3.2 keep parent genotype and corresponding map to be used for
  # clustering. later
  
  geno.par.clu <- geno.par
  
  map.par.clu <- QC_matchMarker(mk.mat = geno.par.clu, map = map)[[2]]
  
  
  ### 3.3 monomorphic or missing in the offsprings
  
  off.MAF <- QC_MAF(mk.mat = geno.off, parallel = parallel,
                    cluster = cluster)
  
  mono <- which(off.MAF == 0)
  miss <- which(is.na(off.MAF))
  
  prob.mk.id <- c(mono, miss)
  
  if(length(prob.mk.id) > 0){
    
    prob.mk.list <- c(prob.mk.list, colnames(geno.par)[prob.mk.id])
    
    geno.par <- geno.par[, -prob.mk.id]
    geno.off <- geno.off[, -prob.mk.id]
    
  }
  
  # 4. Remove marker with problematic MAF at the population level
  ###############################################################
  
  MAF.pop <- QC_MAF(mk.mat = rbind(geno.par, geno.off), parallel = parallel,
                    cluster = cluster)
  
  prob.mk.id <- which(MAF.pop < MAF.pop.lim)
  
  if(length(prob.mk.id) > 0){
    
    prob.mk.list <- c(prob.mk.list, colnames(geno.par)[prob.mk.id])
    
    geno.par <- geno.par[, -prob.mk.id]
    geno.off <- geno.off[, -prob.mk.id]
    
  }
  
  # 5. Remove markers with too high missing rate at the population level
  ######################################################################
  
  miss.ind.mk <- QC_missing(mk.mat = geno.off, threshold = mk.miss)
  
  if(dim(miss.ind.mk)[1] > 0){
    
    prob.mk.list <- c(prob.mk.list, colnames(geno.par)[miss.ind.mk[, 2]])
    
    geno.par <- geno.par[, -miss.ind.mk[, 2]]
    geno.off <- geno.off[, -miss.ind.mk[, 2]]
    
  }
  
  # 6. Remove less polymorphic marker(s) if some markers are at the same position
  ###############################################################################
  
  # select maximum 1 marker per position
  
  map <- QC_matchMarker(mk.mat = geno.off, map = map)[[2]]
  
  difference <- diff(map[, 3])
  difference <- c(1, difference) # add 1 for the first position.
  
  if(sum(difference == 0) > 0){
    
    # Identify blocks of marker that are at the same position.
    
    list.pos <- list()
    i <- 1
    max <- length(difference)
    list.pos.id <- 1
    
    while(i <= max){
      
      # test if the value is zero -> same position
      
      if(difference[i] == 0){
        
        vec <- c(i-1, i) # start the vector
        i <- i + 1
        
        if(i > max){
          list.pos[[list.pos.id]] <- vec
          break()
        }
        
        while(difference[i] == 0){ # continue as long as the values are zero
          
          vec <- c(vec, i)
          i <- i + 1
          
          if(i > max){
            list.pos[[list.pos.id]] <- vec
            break()
          }
          
        } # see if there are other positions
        
        list.pos[[list.pos.id]] <- vec # store the result
        list.pos.id <- list.pos.id + 1
        
      } else {
        
        i <- i + 1
        
      }
      
    }
    
    # identify the most polymorphic marker within the pre-selected blocks and
    # therefore the positions to remove.
    
    ref.MAF <- MAF.pop[colnames(geno.off)]
    
    prob.mk.id <- c()
    
    for(i in 1:length(list.pos)){
      
      rem.mk_i <- list.pos[[i]][-which.max(ref.MAF[list.pos[[i]]])]
      
      prob.mk.id <- c(prob.mk.id, rem.mk_i)
      
    }
    
    prob.mk.list <- c(prob.mk.list, colnames(geno.par)[prob.mk.id])
    
    geno.par <- geno.par[, -prob.mk.id]
    geno.off <- geno.off[, -prob.mk.id]
    map <- map[-prob.mk.id, ]
    
  }
  
  
  # 7. Remove genotypes with too high missing rate at the population level
  ########################################################################
  
  geno.ref <- rownames(geno.off) # make a reference list of genotypes
  
  miss.ind.gen <- QC_missing(mk.mat = geno.off, MARGIN = 1,
                             threshold = gen.miss)
  
  if(dim(miss.ind.gen)[1] > 0){
    
    geno.off <- geno.off[-miss.ind.gen[, 2], ]
    
    # adapt the other arguments which depend on the genotype list
    
    ind.geno <- geno.ref %in% rownames(geno.off)
    cross.ind <- cross.ind[ind.geno]
    trait <- trait[ind.geno, ]
    
    if(!is.null(subcross.ind)) {
      
      subcross.ind <- subcross.ind[ind.geno]
      
    }
    
  }
  
  # 8. Remove cross with a too small size
  #######################################
  
  geno.ref <- rownames(geno.off) # make a reference list of genotypes
  
  geno.off <- QC_minCrSize(mk.mat = geno.off, cross.ind = cross.ind,
                           n.lim = n.lim)
  
  # adapt the other arguments which depend on the genotype list
  
  ind.geno <- geno.ref %in% rownames(geno.off)  
  cross.ind <- cross.ind[ind.geno]
  trait <- trait[ind.geno, ]
  
  if(!is.null(subcross.ind)) {
    
    subcross.ind <- subcross.ind[ind.geno]
    
  }
  
  # 9. Remove marker with problematic MAF at the cross level
  ##########################################################
  
  MAF.pop <- QC_MAF(mk.mat = geno.off, cross.ind = cross.ind,
                    parallel = parallel, cluster = cluster)
  
  within.cr.MAF <- MAF.pop[[2]]
  
  # functions to determine the MAF limit within crosses
  
  if(is.null(MAF.cr.lim)){
    
    MAF.lim <- function(floor, n.cr){
      
      (5/n.cr) + floor
      
    }
    
    # determine the number of observation per cross. First transform into
    # factor with specified order
    
    n.cr <- table(factor(cross.ind, levels = unique(cross.ind)))
    
    lim <- MAF.lim(floor = 0.05, n.cr = n.cr)
    
  } else {
    
    lim <- MAF.cr.lim
    
  }
  
  MAF.cr.ind <- QC_tagMAFCr(MAF = MAF.pop, MAF.lim = lim, tag.mono = FALSE,
                            parallel = parallel, cluster = cluster)
  
  # two options to manage the marker with problementic MAF within cross.
  # 1.: Put these markers as missing within the cross; 2.: remove the marker 
  
  if(MAF.cr.miss){ # put NA markers with prob. within cross MAF in at least 1 cross
    
    cr.id <- unique(cross.ind)
    
    for(i in 1:dim(geno.off)[2]){
      
      test <- (within.cr.MAF[, i] < lim) & (within.cr.MAF[, i] != 0)
      test[is.na(test)] <- FALSE
      
      if(sum(test) > 0){ # at least one cross has a problematic MAF
        
        geno.off[cross.ind %in% cr.id[test], i] <- NA
        
      }
      
    }
    
    if (rem.NA.cr){
      
      prob.mk.id <- which(is.na(MAF.cr.ind))
      
      if(length(prob.mk.id) > 0){
        
        prob.mk.list <- c(prob.mk.list, colnames(geno.par)[prob.mk.id])
        
        geno.par <- geno.par[, -prob.mk.id]
        geno.off <- geno.off[, -prob.mk.id]
        
      }
      
    } 
    
  } else { # remove markers with problematic within cross MAF in at least 1 cross
    
    
    if (rem.NA.cr){
      
      miss.cr <- which(is.na(MAF.cr.ind))
      prob.mk <- which(MAF.cr.ind)
      
      prob.mk.id <- c(miss.cr, prob.mk)
      
      if(length(prob.mk.id) > 0){
        
        prob.mk.list <- c(prob.mk.list, colnames(geno.par)[prob.mk.id])
        
        geno.par <- geno.par[, -prob.mk.id]
        geno.off <- geno.off[, -prob.mk.id]
        
      }
      
    } else {
      
      prob.mk.id <- which(MAF.cr.ind)
      
      if(length(prob.mk.id) > 0){
        
        prob.mk.list <- c(prob.mk.list, colnames(geno.par)[prob.mk.id])
        
        geno.par <- geno.par[, -prob.mk.id]
        geno.off <- geno.off[, -prob.mk.id]
        
      } 
      
    }
    
  }
  
  
  # 10. equalize the list of markers and the list of genotypes
  ############################################################
  
  # adjust all object that must be modified when the list of
  # genotype decrease: cross.ind (subcross.ind), par.per.cross
  # geno.par (if one cross is removed the parent that appear only in
  # this cross must also be removed.
  
  ### 10.1 markers
  
  match <- QC_matchMarker(mk.mat = geno.off, map = map)
  map <- match$new.map
  
  
  
  
  ### 10.2 genotypes (equalize also cross.ind and subcross.ind)
  
  if(!is.null(subcross.ind)){
    
    trait.aug <- data.frame(trait, cross.ind, subcross.ind,
                            stringsAsFactors = FALSE)
    
    rownames(trait.aug) <- trait.aug[, 1]
    
    match <- QC_matchGeno(mk.mat = geno.off, pheno = trait.aug)
    trait <- match$new.pheno[, 1:2]
    cross.ind <- match$new.pheno[, 3]
    subcross.ind <- match$new.pheno[, 4]
    rm(match)
    
  } else {
    
    trait.aug <- data.frame(trait, cross.ind, stringsAsFactors = FALSE)
    
    rownames(trait.aug) <- trait.aug[, 1]
    
    match <- QC_matchGeno(mk.mat = geno.off, pheno = trait.aug)
    trait <- match$new.pheno[, 1:2]
    cross.ind <- match$new.pheno[, 3]
    rm(match)  
    
  }
  
  ### 10.3 modify the par.per.cross argument and geno.par
  
  cr.list <- unique(cross.ind)
  par.per.cross <- par.per.cross[par.per.cross[, 1] %in% cr.list, ]
  
  par.list <- union(par.per.cross[, 2], par.per.cross[, 3])
  geno.par <- geno.par[rownames(geno.par) %in% par.list, ]
  
  # Modify also the geno.par.clu by removing unused parents
  
  geno.par.clu <- geno.par.clu[rownames(geno.par.clu) %in% par.list,]
  
  # 11. optional ABH assignement
  #############################
  
  if(ABH){
    
    # check if a cross was completely removed. Then remove it from the
    # par.per.cross argument.
    
    if(!is.null(subcross.ind)){
      
      cr.list <- unique(subcross.ind)
      par.per.cross.ABH <- par.per.subcross[par.per.subcross[, 1] %in% cr.list, ]
      cross.ind.ABH <- subcross.ind
      
    } else {
      
      par.per.cross.ABH <- par.per.cross 
      cross.ind.ABH <- cross.ind
      
    }
    
    ### 11.1 case with heterozygous markers scores.
    
    if(het.par){
      
      # ABH assignement with heterogeneous parents
      
      geno.off <- cross_ABH_het(par.sc = geno.par, off.sc = geno.off,
                                cross.ind = cross.ind.ABH,
                                par.per.cross = par.per.cross.ABH)
      
      
    } else {
      
      geno.off <- cross_ABH(par.sc = geno.par, off.sc = geno.off,
                            cross.ind = cross.ind.ABH,
                            par.per.cross = par.per.cross.ABH)
      
    }
    
  }
  
  
  # 12. Results
  #############
  
  results <- list(geno.off = geno.off, geno.par = geno.par,
                  geno.par.clu = geno.par.clu, map.par.clu = map.par.clu,
                  cross.ind = cross.ind, par.per.cross = par.per.cross,
                  trait = trait, map = map,  rem.mk = prob.mk.list)
  
  return(results)
  
}
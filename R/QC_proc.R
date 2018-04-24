###########
# QC_proc #
###########

#' Quality control procedure
#'
#' Perform different operations of quality control (QC) on the marker matrix and
#' return the marker matrix, the map and the phenotype data in format to form a
#' \code{mppData} object.
#' 
#' The different operations of the quality control are the following:
#' 
#' \enumerate{
#' 
#' \item{Remove markers with more than two alleles
#' (\code{\link{QC_GenotypingError}}).}
#' 
#' \item{Remove markers that are monomorphic or fully missing in the parents
#' (\code{\link{QC_MAF}}).}
#' 
#' \item{Remove markers with a missing rate higher than \code{mk.miss} across
#' the entire MPP (\code{\link{QC_missing}}).}
#' 
#' \item{Remove genotypes with more missing markers than \code{gen.miss}
#' (\code{\link{QC_missing}}).}
#' 
#' \item{Remove crosses with less than \code{n.lim} genotypes
#' (\code{\link{QC_minCrSize}}).}
#' 
#' \item{Keep only the most polymorphic marker when multiple markers map at the
#' same position.} 
#' 
#' \item{Remove markers with a minor allele frequency (MAF) below a threshold
#' given by \code{MAF.pop.lim} (\code{\link{QC_MAF}}).}
#' 
#' \item{Determine markers having a too low within cross MAF
#' (\code{\link{QC_MAF}} and \code{\link{QC_tagMAFCr}}). The user can give the
#' critical values for MAF within cross using \code{MAF.cr.lim}.
#' 
#' By default, the within cross MAF values are defined by the following function
#' of the cross-size n.ci: MAF(n.ci) = 0.5 if n.ci c [0, 10] and MAF(n.ci) =
#' (4.5/n.ci) + 0.05 if n.ci > 10. This means that up to 10 genotypes,
#' the critical within cross MAF is set to 50%. Then it decreases when the
#' number of genotype increases until 5% set as a lower bound. 
#' 
#' If the within cross MAF is below the limit in at least one cross, then marker
#' scores of the problematic cross are either put as missing
#' (\code{MAF.cr.miss = TRUE}) or the whole marker is discarded
#' (\code{MAF.cr.miss = FALSE}). By default, \code{MAF.cr.miss = TRUE} which
#' allows to include a larger number of markers and to cover a wider genetic
#' diversity.}
#' 
#' \item{If \code{ABH = TRUE}, convert offspring genotype data into ABH format
#' (\code{\link{cross_ABH}} or \code{\link{cross_ABH_het}}).}
#' 
#' \item{If \code{impute = TRUE}, imputation of the missing values for the
#' IBS bi-allelic model using the \code{codeGeno()} function from synbreed
#' (Wimmer et al., 2012). The imputation can only be done if data are not
#' converted into ABH format (\code{ABH = FALSE}).}
#' 
#' }
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
#' the one used in \code{par.per.cross}. The columns names must be the marker
#' names similar to the one used in \code{map}. Marker scores must be coded
#' using one letter per allele. For example, AA, CC, GG, TT, AC, AG, AT, CA,
#' CG, CT, GA, GC, GT, TA, TC, TG. Missing values must be coded \code{NA}.}
#' 
#' @param map Three columns \code{data.frame} with: 1) marker or in between
#' position identifiers; 2) chromosome; 3) positions in centi-Morgan.\strong{
#' The marker identifiers must be identical to the column names of the maker
#' matrices (\code{geno.off} and \code{geno.par}).}
#' 
#' @param trait Two columns \code{data.frame} with : 1) \code{character}
#' genotypes identifiers; 2) \code{numeric} trait values. \strong{The genotypes
#' identifiers must be identical to the rownames of the offspring marker matrix
#' (\code{geno.off}).}
#' 
#' @param cross.ind \code{Character} vector indicating to which cross each
#' genotype belongs.
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
#' \code{par.per.subcross} must also be provided.}. Default = NULL.
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
#' Default = 15.
#' 
#' @param MAF.pop.lim \code{Numeric} value specifying the minimum minor allele
#' frequency for a marker at the population level. Default = 0.05.
#' 
#' @param MAF.cr.lim \code{Numeric vector} specifying the critical within cross
#' MAF. Marker with a problematic segregation rate in at least
#' one cross is either set as missing within the problematic cross
#' (\code{MAF.cr.miss = TRUE}), or remove from the marker matrix
#' (\code{MAF.cr.miss = FALSE}). For default value see details.
#' 
#' @param MAF.cr.miss Logical value specifying if maker with a too low
#' segregation rate within cross should be put as missing or discarded. If
#' \code{MAF.cr.miss = TRUE}, marker scores with a problematic segregation rate
#' in at least one cross will be put as missing in the problematic cross(es).
#' If \code{MAF.cr.miss = FALSE}, the whole marker will be discarded.
#' Default = TRUE.
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
#' the markers to the ABH format and parent have heterozygous or missing marker
#' scores, the argument \code{het_miss_par} must be set to TRUE}.
#' If (\code{ABH = FALSE}), the marker data will be let in the same format as
#' the one used in \code{geno.off}. Default = TRUE.
#' 
#' @param het_miss_par \code{Logical} value. if \code{het_miss_par = TRUE},
#' the function will use the offspring segregation to try to infer the allele
#' that was transmitted by the heterozygous or missing parent at a particular
#' locus in order to make the ABH conversion. Default = FALSE.
#' 
#' @param verbose \code{Logical} value indicating if the steps of the QC should
#' be printed. Default = TRUE.
#' 
#' @param impute \code{Logical} value. if \code{impute = TRUE}, the function
#' will impute missing values using the \code{codeGeno()} function from the
#' synbreed package. Default = FALSE.
#' 
#' @param impute.type \code{character} with one out of \code{"fix"},
#' \code{"random"}, \code{"family"}, \code{"beagle"}, \code{"beagleAfterFamily"},
#' \code{"beagleAfterFamilyNoRand"}. For details see synbreed package
#' documentation. \strong{To be able to use Beagle for imputation, Please load
#' the synbreed package using \code{library(synbreed)}} Default = "random".
#' 
#' @param map_bp \code{data.frame} with three columns specifying for each marker
#' position the marker identifier, the \code{numeric} or \code{character}
#' chromosome the and physical bp position. This argument is necessary for
#' imputation using Beagle. Default = NULL.
#' 
#' @param replace.value \code{numeric} scalar to replace missing value in case
#' \code{impute.type = fix}. Only 0, 1, 2. Should be chosen. Default = NULL.
#' 
#' @param label.heter This is either a scalar or vector of characters to identify
#' heterozygous genotypes or a function returning TRUE if an element of the
#' marker matrix is the heterozygous genotype. Defining a function is useful,
#' if number of unique heterozygous genotypes is large, i.e. if genotypes are
#' coded by alleles. If the heterozygous genotype is coded like
#' "A/T","G/C", ..., "AG", "CG", ..., "T:C", "G:A", ... or "G|T", "A|C", ...
#' then label.heter="alleleCoding" can be used. Note that
#' heterozygous values must be identified unambiguously by label.heter. Use
#' label.heter=NULL if there are only homozygous genotypes, i.e. in DH lines,
#' to speed up computation and restrict imputation to values 0 and 2.
#' Default = "alleleCoding".
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
#' \item{rem.geno}{Vector of genotypes that have been removed.}
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
#' @references
#' 
#' Wimmer, V., Albrecht, T., Auinger, H. J., & Schon, C. C. (2012). synbreed: a
#' framework for the analysis of genomic prediction data using R.
#' Bioinformatics, 28(15), 2086-2087.
#' 
#' Browning, B. L., & Browning, S. R. (2013). Improving the accuracy and
#' efficiency of identity-by-descent detection in population data. Genetics,
#' 194(2), 459-471.
#' 
#' @examples
#' 
#' data(USNAM_geno)
#' data(USNAM_pheno)
#' data(USNAM_map)
#' 
#' cross.ind <- substr(rownames(USNAM_pheno), 1, 4)
#' geno.off <- USNAM_geno[7:506, ]
#' geno.par <- USNAM_geno[1:6, ]
#' 
#' map <- USNAM_map
#' trait <- data.frame(rownames(USNAM_pheno), USNAM_pheno[, 1],
#'                    stringsAsFactors = FALSE)
#' colnames(trait) <- c('genotypes', 'ULA')
#' rownames(trait) <- rownames(USNAM_pheno)
#' par.per.cross <- cbind(unique(cross.ind), rep("B73", 5),
#'                        rownames(geno.par)[2:6])
#' 
#' # QC IBD model and analysis
#' ###########################
#' 
#' data <- QC_proc(geno.off = geno.off, geno.par = geno.par, map = map,
#'                 trait = trait, cross.ind = cross.ind,
#'                 par.per.cross = par.per.cross, n.lim = 15,
#'                 MAF.pop.lim = 0.05, mk.miss = 0.1,
#'                 gen.miss = 0.25, ABH = TRUE, het_miss_par = TRUE,
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
#'                         IBS = FALSE, type = "F", nb.gen = 6,
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
#' # QC IBS model with imputation
#' ##############################
#' 
#' data(USNAM_geno)
#' data(USNAM_pheno)
#' data(USNAM_map)
#' 
#' cross.ind <- substr(USNAM_pheno[, 1], 1, 4)
#' geno.off <- USNAM_geno[7:506, ]
#' 
#' # It is a RIL population. So set the heterozygous score to missing.
#' 
#' geno.off[geno.off == "AC"] <- NA
#' geno.off[geno.off == "AG"] <- NA
#' geno.off[geno.off == "AT"] <- NA
#' geno.off[geno.off == "CG"] <- NA
#' geno.off[geno.off == "CT"] <- NA
#' geno.off[geno.off == "GT"] <- NA
#' 
#' geno.par <- USNAM_geno[1:6, ]
#' 
#' map <- USNAM_map
#' trait <- data.frame(rownames(USNAM_pheno), USNAM_pheno[, 1],
#'                    stringsAsFactors = FALSE)
#' colnames(trait) <- c('genotypes', 'ULA')
#' rownames(trait) <- rownames(USNAM_pheno)
#' par.per.cross <- cbind(unique(cross.ind), rep("B73", 5),
#'                        rownames(geno.par)[2:6])
#' 
#' 
#' data <- QC_proc(geno.off = geno.off, geno.par = geno.par, map = map,
#'                 trait = trait, cross.ind = cross.ind,
#'                 par.per.cross = par.per.cross, ABH = FALSE, impute = TRUE,
#'                 impute.type = "family")
#'
#' @export
#' 


QC_proc <- function(geno.off, geno.par, map, trait, cross.ind, par.per.cross,
                    subcross.ind = NULL, par.per.subcross = NULL,
                    n.lim = 15, MAF.pop.lim = 0.05, MAF.cr.lim = NULL,
                    MAF.cr.miss = TRUE, mk.miss = 0.1,
                    gen.miss = 0.25, ABH = TRUE, het_miss_par = FALSE,
                    verbose = TRUE, impute = FALSE, impute.type = "random",
                    map_bp = NULL, replace.value = NULL,
                    label.heter = "alleleCoding", parallel = FALSE,
                    cluster = NULL){
  
  # 1. check the format of the data
  #################################
  
  check_QC(geno.off = geno.off, geno.par = geno.par, map = map, trait = trait,
           cross.ind = cross.ind, par.per.cross = par.per.cross,
           subcross.ind = subcross.ind, par.per.subcross = par.per.subcross,
           n.lim = n.lim, MAF.pop.lim = MAF.pop.lim, mk.miss = mk.miss,
           gen.miss = gen.miss, MAF.cr.lim = MAF.cr.lim, ABH = ABH,
           het_miss_par = het_miss_par, impute = impute,
           impute.type = impute.type, map_bp = map_bp,
           replace.value = replace.value, parallel = parallel,
           cluster = cluster)
  
  
  # 2. Remove markers with genotyping error
  #########################################
  
  if(verbose){
    
    message("The quality control procedure can take few minutes!")
    
  }
  
  
  init.nb.mk <- dim(geno.off)[2]
  nb.mk.rem <- c()
  init.nb.gen <- dim(geno.off)[1]
  nb.gen.rem <- c()
  
  prob.mk.list <- c()
  prob.gen.list <- c()
  
  prob.mk <- QC_GenotypingError(mk.mat = rbind(geno.par, geno.off),
                                parallel = parallel, cluster = cluster)
  
  if(is.null(prob.mk)) {rem.mk_i <- 0 } else {rem.mk_i <- length(prob.mk)}
  
  if(!is.null(prob.mk)){
    
    prob.mk.list <- c(prob.mk.list, prob.mk)
    
    ind.prob <- which(colnames(geno.off) %in% prob.mk)
    geno.par <- geno.par[, -ind.prob]
    geno.off <- geno.off[, -ind.prob]
    nb.mk.rem <- c(nb.mk.rem, rem.mk_i)
    
  }
  
  if(verbose){
    
    cat("\n")
    cat(paste("Check genotyping error                        :",
              rem.mk_i, "markers removed", "\n"))
    
  }
  
  
  
  # 3. remove monomorphic markers in the parents
  ##############################################
  
  parent.MAF <- QC_MAF(mk.mat = geno.par, parallel = parallel,
                       cluster = cluster)
  
  mono <- which(parent.MAF == 0)
  miss <- which(is.na(parent.MAF))
  
  prob.mk.id <- c(mono, miss)
  
  if(is.null(prob.mk.id)) {rem.mk_i <- 0 } else {rem.mk_i <- length(prob.mk.id)}
  
  if(length(prob.mk.id) > 0){
    
    prob.mk.list <- c(prob.mk.list, colnames(geno.par)[prob.mk.id])
    
    geno.par <- geno.par[, -prob.mk.id]
    geno.off <- geno.off[, -prob.mk.id]
    nb.mk.rem <- c(nb.mk.rem, rem.mk_i)
    
  }
  
  if(verbose){
    
    cat(paste("Remove monomorphic/missing marker in parents  :",
              rem.mk_i, "markers removed", "\n"))
    
  }
  
  ### keep parent genotype and corresponding map to be used for
  # clustering. later
  
  geno.par.clu <- geno.par
  
  map.par.clu <- QC_matchMarker(mk.mat = geno.par.clu, map = map)[[2]]
  
  
  # 4. Remove markers with too high missing rate at the population level
  ######################################################################
  
  miss.ind.mk <- QC_missing(mk.mat = geno.off, threshold = mk.miss)
  
  if(dim(miss.ind.mk)[1] > 0) {rem.mk_i <- dim(miss.ind.mk)[1]
  } else {rem.mk_i <- 0}
  
  if(dim(miss.ind.mk)[1] > 0){
    
    prob.mk.list <- c(prob.mk.list, colnames(geno.par)[miss.ind.mk[, 2]])
    
    geno.par <- geno.par[, -miss.ind.mk[, 2]]
    geno.off <- geno.off[, -miss.ind.mk[, 2]]
    nb.mk.rem <- c(nb.mk.rem, rem.mk_i)
    
  }
  
  if(verbose){
    
    cat(paste("Remove marker with missing rate <", mk.miss, "        :",
              rem.mk_i, "markers removed", "\n"))
    
  }
  
  # 5. Remove genotypes with too high missing rate at the population level
  ########################################################################
  
  geno.ref <- rownames(geno.off) # make a reference list of genotypes
  
  miss.ind.gen <- QC_missing(mk.mat = geno.off, MARGIN = 1,
                             threshold = gen.miss)
  
  if(dim(miss.ind.gen)[1] > 0) {rem.gen_i <- dim(miss.ind.gen)[1]
  } else {rem.gen_i <- 0}
  
  if(dim(miss.ind.gen)[1] > 0){
    
    geno.off <- geno.off[-miss.ind.gen[, 2], ]
    nb.gen.rem <- c(nb.gen.rem, rem.gen_i)
    prob.gen.list <- c(prob.gen.list, as.character(miss.ind.gen[, 1]))
    
    # adapt the other arguments which depend on the genotype list
    
    ind.geno <- geno.ref %in% rownames(geno.off)
    cross.ind <- cross.ind[ind.geno]
    trait <- trait[ind.geno, ]
    
    if(!is.null(subcross.ind)) {
      
      subcross.ind <- subcross.ind[ind.geno]
      
    }
    
  }
  
  if(verbose){
    
    cat(paste("Remove genotype with missing rate <", gen.miss, "     :",
              rem.gen_i, "genotypes removed", "\n"))
    
  }
  
  
  # 6. Remove cross with a too small size
  #######################################
  
  geno.ref <- rownames(geno.off) # make a reference list of genotypes
  
  geno.off <- QC_minCrSize(mk.mat = geno.off, cross.ind = cross.ind,
                           n.lim = n.lim)
  
  rem.gen_i <- length(geno.ref) - dim(geno.off)[1]
  
  if(rem.gen_i > 0){
    
    nb.gen.rem <- c(nb.gen.rem, rem.gen_i)
    gen.removed <- geno.ref[!(geno.ref %in% rownames(geno.off))]
    prob.gen.list <- c(prob.gen.list, gen.removed)
    
    # adapt the other arguments which depend on the genotype list
    
    ind.geno <- geno.ref %in% rownames(geno.off)  
    cross.ind <- cross.ind[ind.geno]
    trait <- trait[ind.geno, ]
    
    if(!is.null(subcross.ind)) {
      
      subcross.ind <- subcross.ind[ind.geno]
      
    }
    
  }
  
  if(verbose){
    
    cat(paste("Remove crosses with less than", n.lim, "observations :",
              rem.gen_i, "markers removed", "\n"))
    
  }
  
  
  # 7. Remove less polymorphic marker(s) if some markers are at the same position
  ###############################################################################
  
  # select maximum 1 marker per position
  
  map <- QC_matchMarker(mk.mat = geno.off, map = map)[[2]]
  
  difference <- diff(map[, 3])
  difference <- c(1, difference) # add 1 for the first position.
  
  rem.mk_j <- sum(difference == 0)
  
  if(rem.mk_j > 0){
    
    MAF.pop <- QC_MAF(mk.mat = geno.off, parallel = parallel, 
                      cluster = cluster) 
    
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
    
    nb.mk.rem <- c(nb.mk.rem, rem.mk_j)
    
  }
  
  if(verbose){
    
    cat(paste("Remove markers at the same position           :",
              rem.mk_j, "markers removed", "\n"))
    
  }
  
  # 8. Remove markers with too low MAF at the population level
  ############################################################
  
  off.MAF <- QC_MAF(mk.mat = geno.off, cross.ind = cross.ind,
                    parallel = parallel, cluster = cluster)
  
  
  MAF.pop <- off.MAF[[1]]
  MAF.cr <- off.MAF[[2]]
  
  prob.mk.id <- which(MAF.pop < MAF.pop.lim)
  
  if(is.null(prob.mk.id)) {rem.mk_i <- 0 } else {rem.mk_i <- length(prob.mk.id)}
  
  
  if(length(prob.mk.id) > 0){
    
    prob.mk.list <- c(prob.mk.list, colnames(geno.par)[prob.mk.id])
    
    geno.par <- geno.par[, -prob.mk.id]
    geno.off <- geno.off[, -prob.mk.id]
    MAF.cr <- MAF.cr[, -prob.mk.id]
    MAF.pop <- MAF.pop[-prob.mk.id]
    nb.mk.rem <- c(nb.mk.rem, rem.mk_i)
    
  }
  
  if(verbose){
    
    cat(paste("Remove markers with MAF <", MAF.pop.lim,"               :",
              rem.mk_i, "markers removed", "\n"))
    
  }
  
  # 9. Remove marker with problematic within cross MAF
  ####################################################
  
  
  MAF.pop <-  list(MAF.pop = MAF.pop, MAF.cr = MAF.cr)
  class(MAF.pop) <- c("list", "mafRes")
  
  # functions to determine the MAF limit within crosses
  
  if(is.null(MAF.cr.lim)){
    
    MAF.lim <- function(x, floor){
      
      if(x <= 10){ 0.5 } else { (4.5/x) + floor }
      
    }
    
    # determine the number of observation per cross. First transform into
    # factor with specified order
    
    n.cr <- table(factor(cross.ind, levels = unique(cross.ind)))
    
    lim <- unlist(lapply(X = n.cr, FUN = MAF.lim, floor = 0.05))
    
  } else {
    
    lim <- MAF.cr.lim
    
  }
  
  MAF.cr.ind <- QC_tagMAFCr(MAF = MAF.pop, MAF.lim = lim, tag.mono = FALSE,
                            parallel = parallel, cluster = cluster)
  
  # two options to manage the marker with problementic MAF within cross.
  # 1.: Put these markers as missing within the cross; 2.: remove the marker 
  
  prob.mk.id <- which(MAF.cr.ind)
  rem.mk_i <- length(prob.mk.id)
  
  if(MAF.cr.miss){ # put NA markers with prob. within cross MAF in at least 1 cross
    
    cr.id <- unique(cross.ind)
    
    for(i in 1:dim(geno.off)[2]){
      
      test <- (MAF.cr[, i] < lim) & (MAF.cr[, i] != 0)
      test[is.na(test)] <- FALSE
      
      if(sum(test) > 0){ # at least one cross has a problematic MAF
        
        geno.off[cross.ind %in% cr.id[test], i] <- NA
        
      }
      
    }
    
    
  } else { # remove markers with problematic within cross MAF in at least 1 cross
    
    
    if(rem.mk_i > 0){
      
      prob.mk.list <- c(prob.mk.list, colnames(geno.par)[prob.mk.id])
      
      geno.par <- geno.par[, -prob.mk.id]
      geno.off <- geno.off[, -prob.mk.id]
      nb.mk.rem <- c(nb.mk.rem, rem.mk_i)
      
    }
    
    if(verbose){
      
      cat(paste("Remove markers critical witin cross MAF            :",
                rem.mk_i, "markers removed", "\n"))
      
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
  
  ### 10.3 verify that there are no cross with too few data
  
  freq <- table(cross.ind)
  
  if(sum(freq < 15)){
    
    warning(paste("It is still possible to perform a MPP QTL analysis with",
                  "crosses containing less that 15 genotypes. However, we",
                  "advice to use minimum 15 individuals per cross to have",
                  "enough information to estimate within crosses QTL effects."))
    
  }
  
  ### 10.4 modify the par.per.cross argument and geno.par
  
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
    
    if(het_miss_par){
      
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
  
  # 12. Optional marker imputation
  ################################
  
  if(impute){
    
    if(verbose){
      
      cat("\n")
      cat("Missing values imputation \n")
      
    }
    
    # separate imputation for fixed value
    
    if(impute.type == "fix"){
      
      geno012 <- geno_012(mk.mat = geno.off)[[1]]
      geno012[is.na(geno012)] <- replace.value
      geno.off <- geno012
      
    } else {
      
      family <- data.frame(cross.ind)
      rownames(family) <- rownames(geno.off)
      
      if(impute.type %in% c("random", "family", "fix")) {
        
        map_gp <- map[, 2:3]
        map.unit <- "cM"
        n.cores <- 1
        
      } else { # Cases with Beagle
        
        message("The imputation using Beagle can take several minutes!")
        
        map_gp <- map_bp[map_bp[, 1] %in% map[, 1], 2:3]
        map.unit <- "bp"
        if(!is.null(cluster)){n.cores <- length(cluster)} else {n.cores <- 1}
        
      }
      
      colnames(map_gp) <- c("chr", "pos")
      rownames(map_gp) <- map[, 1]
      
      
      gp <- create.gpData(geno = geno.off, map = map_gp, family = family,
                          map.unit = map.unit)
      
      gp.imp <- tryCatch(expr = codeGeno(gpData = gp, impute = TRUE,
                                         impute.type = impute.type,
                                         replace.value = replace.value, maf = 0,
                                         nmiss = 1, label.heter = label.heter,
                                         verbose = FALSE, cores = n.cores),
                         error = function(e) NULL)
      
      
      if(is.null(gp.imp)){
        
        if(map.unit == "cM"){ warning("The missing values imputation failed!")
          
        } else{ warning("The missing values imputation using Beagle failed!") }
        
      } else {geno.off <- gp.imp$geno }
      
    }
    
    
    
  }
  
  # 13. Results
  #############
  
  results <- list(geno.off = geno.off, geno.par = geno.par,
                  geno.par.clu = geno.par.clu, map.par.clu = map.par.clu,
                  cross.ind = cross.ind, par.per.cross = par.per.cross,
                  trait = trait, map = map,  rem.mk = prob.mk.list,
                  rem.geno = prob.gen.list)
  
  ###### final message
  
  if(verbose){
    
    tot.rem.mk <- init.nb.mk - sum(nb.mk.rem)
    tot.rem.gen <- init.nb.gen - sum(nb.gen.rem)
    
    cat("\n")
    cat("   End             :", tot.rem.mk ,
        "marker(s) remain after the check\n")
    cat("                    ",
        tot.rem.gen , "genotypes(s) remain after the check\n")
    
    
  }
  
  return(results)
  
}
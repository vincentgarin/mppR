################
# mppData_form #
################

#' Data for MPP QTL analysis
#' 
#' Form a single data object to perform MPP QTL anlyses.
#' 
#' \code{mppData_form} combines into a single data object the all the necessary
#' data to run QTL analyses except the parent clustering results
#' (\code{\link{parent_cluster}}). Two types of \code{mppData} objects can be
#' produced. The first one (\code{IBS = FALSE}), for cross-specific, parental and
#' ancestral models calculates identical by descent (IBD) probabilities using
#' \code{read.cross()} and \code{calc.genoprob()} functions from the R/qtl package
#' (Broman et al. 2009). F-type (F), back-cross (BC), double
#' haploid (DH) and recombinant imbred lines (RIL) are the allowed population type.
#' DH and RIL populations are read as back-cross by R/qtl. For these two
#' population types, no heterozygosity is allowed. It must be converted into missing.
#' \code{\link{USNAM_mppData}} is an example of an IBD \code{mppData} object.
#' 
#' The second type of \code{mppData} objects (\code{IBS = TRUE}) uses identical
#' by state (IBS) genetic predictors for bi-allelic model. In such case,
#' genotype marker score are transformed into 0, 1, 2 format representing
#' the copy number of the least frequent allele. For IBS \code{mppData} objects,
#' the user can provide imputed data for the argument \code{geno.off}. Marker
#' score imputation can be done using \code{\link{QC_proc}}. The markers scores
#' can be provided to \code{geno.off} using letters (1 letter per allele)
#' setting (\code{IBS.format = "ATCG"}) or already in 0, 1, 2 format
#' (\code{IBS.format = "012"}). \code{\link{USNAM_mppData_bi}} is
#' an example of an IBS \code{mppData} object.
#' 
#' 
#' @param geno.off \code{Character} or \code{Numeric matrix} representing
#' genotype marker scores with genotypes as row and markers as column. For
#' cross-specific, parental and ancestral models (\code{IBS = FALSE}) marker
#' data must be in ABH format (for details see \code{\link{cross_ABH}}).
#' 
#' For bi-allelic models (\code{IBS = TRUE}), the marker scores in
#' \code{geno.off} can be coded using letters using one letter for each allele.
#' For example, AA, CC, GG, TT, AC, AG, AT, CA, CG, CT, GA, GC, GT, TA, TC, TG
#' (\code{IBS.format = "ATCG"}). \code{geno.off} marker score can also be
#' already coded in 0, 1, 2 format (\code{IBS.format = "012"}). 
#' 
#' For both type of mppData object (bi-allelic or not), \strong{the row names
#' must be the genotypes identifiers similars to the one of
#' trait argument. The column names must be the marker
#' identifiers similar to the first column of the map. Missing value
#' must be coded \code{NA}.}
#' 
#' @param geno.par Optional argument used as indication in the results for the
#' QTL genetic effects. \code{geno.par} is a \code{character matrix}
#' representing the marker scores of the parents with genotypes as row and
#' markers as column. \strong{The column names must be similar as the marker
#' list in \code{map} and \code{geno.off} arguments. The rownames must represent
#' the parents identifiers and be identical to the parent list given in
#' \code{par.per.cross}}. If \code{IBS.format = "012"}, \code{geno.par} is
#' set to NULL. Default = NULL.
#' 
#' @param IBS \code{Logical} value. Put \code{IBS = TRUE} if you want to
#' obtain an IBS \code{mppData} object for the bi-allelic model.
#' Default = FALSE.
#' 
#' @param IBS.format \code{Character} indicator for the format of \code{geno.off}
#' marker score for an IBS \code{mppData} object. \code{IBS.format = "ATCG"}, if
#' marker scores are coded with one letter per allele. \code{IBS.format = "012"},
#' if marker scores are coded in 0, 1, 2 format. Default = "ATCG". 
#' 
#' @param type \code{Character} indicator for the type of population analysed:
#' type = "F" for Fn (F cross n gernerations); type = "BC" for BCn (backcross
#' n generations); type = "BCsFt" for backcross followed by selfing;
#' type = "DH" for double haploids; and type = "RIL"
#' for recombinant inbred lines. For RIL type specify if the population was
#' obtain using selfing or sibling mating using \code{type.mating}.
#' If type = "RIL" or "DH", heterozygous marker scores of the geno.off argument
#' must be set as missing (NA).
#' 
#' @param F.gen \code{Numeric} integer representing the number of F generations.
#' For example F.gen = 2 for F2. Default = NULL.
#' 
#' @param BC.gen \code{Numeric} integer representing the number of
#' backcross generations. For example BC.gen = 1 for single backcross.
#' Default = NULL.
#' 
#' @param type.mating \code{Character} specifying for a RIL population if it was
#' obtained by selfing ("selfing") or by sibling mating ("sib.mat").
#' Default = NULL.
#' 
#' @param map Three columns \code{data.frame} with: 1) marker or in between
#' position identifiers; 2) chromosome; 3) positions in centi-Morgan.
#' \strong{The marker identifiers must be identical to the column names of
#' \code{geno.off}}.
#' 
#' @param trait Two columns \code{data.frame} with : 1) \code{character}
#' genotypes identifiers; 2) \code{numeric} trait values. \strong{The genotypes
#' identifiers must be identical to the rownames of \code{geno.off}).}
#' 
#' @param cross.ind \code{Character} vector with the same length as the number of
#' genotypes which specifies to which cross each genotype belongs.
#' 
#' @param par.per.cross Three columns \code{Character matrix} specifying :
#' 1) the cross indicators (\strong{The cross indicators must be  similar to
#' the one used in \code{cross.ind} and appear in the same order}); 2) the
#' parents 1 identifiers of the crosses; 3) the parents 2 identifiers of the
#' crosses. \strong{The list of parent identifiers must be similar to the
#' rownames of \code{geno.par}}.
#' 
#' @param step \code{Numeric} value representing the maximum distance (in cM)
#' between positions at which the genotype probabilities are calculated.
#' Default value = 10000 (the function do not add position in between markers).
#'  
#' @param error.prob \code{Numeric} value for assumed genotyping error rate
#' used in the calculation of the penetrance Pr(observed genotype | true genotype).
#' Default = 0.0001.
#' 
#' @param map.function \code{Character} expression specifying the type of map
#' function used to infer the IBD probabilities. possibility to choose
#' between "haldane", "kosambi","c-f","morgan".
#' Parameter for the function \code{calc.genoprob()}. Default = "haldane".
#' 
#' @param stepwidth \code{Character} expression indicating if the intermediate
#' points should be fixed ("fixed") or if the algorithm should insert a minimal
#' number of points in between positions ("max"). If stepwidth = "max" the
#' maximum distance between points is step. Default = "max".
#' 
#' @param dir Path where a .csv file Cross_object.csv will be stored to be
#' used for cross object formation. By default, the function uses the current
#' working directory.
#'
#' @return 
#' 
#' Return: a {list} of class \code{mppData} containing the following items :
#' 
#' \item{geno}{For a cross-specific, parental and ancestral model
#' (\code{IBS = FALSE}), cross object obtained with R/qtl functions
#' The main elements of the cross object are the IBD probabilities. For a
#' bi-allelic model (\code{IBS = TRUE}), marker score matrix recoded
#' as 0, 1, 2 according to the copies number of the most frequent allele.}
#' 
#' \item{allele.ref}{If \code{IBS = TRUE}, \code{matrix} with reference allele
#' scores. The first row represents the minor allele (lowest frequency),
#' the second the one represent the major allele (largest frequency) and the
#' two others the heterozygous scores.}
#' 
#' \item{geno.id}{ \code{Character} vector of genotpes identifiers.}
#' 
#' \item{geno.par}{Genotypic map and parent genotype marker matrix if argument
#' \code{geno.par} was provided by the user. For a cross-specific, parental
#' and ancestral model, if positions in between markers were added during
#' the IBD computation, these position are filled as "-".}
#' 
#' \item{map}{\code{Data.frame} with four columns: 1) Marker or in between
#' position names; 2) chromosomes; 3) position numbers on the chromosome;
#' 4) positions in cM.}
#' 
#' \item{trait}{\code{Data.frame} with phenotypic values.}
#' 
#' \item{cross.ind}{\code{Character} vector specifying to which cross each
#' genotype belongs.}
#' 
#' \item{ped.mat}{Four columns \code{data.frame}: 1) the type of genotype:
#' "offspring" for the last genration and "founder" for the genotypes above
#' the offspring in the pedigree; 2) the genotype indicator; 3-4) the parent 1
#' (2) of each line.}
#' 
#' \item{par.per.cross}{Same object as the one given by the user in
#' \code{par.per.cross}.} 
#' 
#' \item{parents}{\code{Character} vector of parents.}  
#' 
#' \item{n.cr}{Number of crosses.}
#' 
#' \item{n.par}{Number of parents.}
#' 
#' \item{type}{\code{Character} expression indicating the type of population.}
#' 
#' \item{n.zigo}{\code{Numeric} value Indicating the number of differenct
#' genotypes: 2 (AA/BB) or 3 (AA/AB/BB)}
#' 
#' \item{biall}{\code{Logical} value specifying if data are made for a
#' bi-allelic model}
#' 
#'  
#' @author Vincent Garin
#' 
#' @seealso \code{\link{cross_ABH}}, \code{\link{mppData_subset}},
#' \code{\link{pedigree_update.mppData}},
#' \code{\link{mpp_SIM}}, \code{\link{parent_cluster}},
#' \code{\link{USNAM_mppData}}, \code{\link{USNAM_mppData_bi}}
#' 
#' @references 
#' 
#' Broman KW, Wu H, Sen S, Churchill GA (2003) R/qtl: QTL mapping
#' in experimental crosses. Bioinformatics 19:889-890.
#' 
#' Broman, K. W., & Sen, S. (2009). A Guide to QTL Mapping with R/qtl (Vol. 46).
#' New York: Springer.
#' 
#'  
#' @examples
#'  
#' # Mpp data for cross, parental or ancestral model
#' #################################################
#' 
#' # Genotypes
#' data(USNAM_geno)
#' data(USNAM_genoABH)
#' # For details about ABH data transformation see cross_ABH() examples
#' 
#' geno.par <- USNAM_geno[1:6, colnames(USNAM_genoABH)]
#' 
#' # Phenotypes
#' data(USNAM_pheno)
#' 
#' trait <- data.frame(rownames(USNAM_pheno), USNAM_pheno[, 1],
#'                    stringsAsFactors = FALSE)
#' colnames(trait) <- c('genotypes', 'ULA')
#' rownames(trait) <- rownames(USNAM_pheno)
#' 
#' # Map
#' data(USNAM_map)
#' 
#' # Same list of markers as in the genotype matrix
#' map <- USNAM_map[USNAM_map[, 1] %in% colnames(USNAM_genoABH),]
#' 
#' cross.ind <- substr(rownames(USNAM_pheno), 1, 4)
#' 
#' par.per.cross <- cbind(unique(cross.ind), rep("B73", 5), rownames(geno.par)[-1])
#' 
#' \dontrun{
#' 
#' # Specify path to a directory
#' my.dir <- "C:/.../"
#' 
#' # IBD mppData for cross, parental or ancestral model
#' data <- mppData_form(geno.off = USNAM_genoABH, geno.par = geno.par,
#'                      IBS = FALSE, type = "F", F.gen = 6, map = map,
#'                      trait = trait, cross.ind = cross.ind,
#'                      par.per.cross = par.per.cross, step = 5,
#'                      map.function = "haldane",  stepwidth = "max",
#'                      dir = my.dir)
#' }
#'
#' 
#' # Mpp data for bi-allelic model
#' ###############################
#' 
#' # Genotypes
#' data(USNAM_geno)
#' 
#' geno <- USNAM_geno[, QC_MAF(USNAM_geno) != 0] # remove monomorphic positions
#' geno.par <- geno[1:6, ]
#' 
#' # Phenotypes
#' data(USNAM_pheno)
#' 
#' trait <- data.frame(rownames(USNAM_pheno), USNAM_pheno[, 1],
#'                    stringsAsFactors = FALSE)
#' colnames(trait) <- c('genotypes', 'ULA')
#' rownames(trait) <- rownames(USNAM_pheno)
#' 
#' # Map
#' data(USNAM_map)
#' 
#' # Same list of markers as in the genotype matrix
#' map <- USNAM_map[USNAM_map[, 1] %in% colnames(geno),]
#' 
#' cross.ind <- substr(USNAM_pheno[, 1], 1, 4)
#' 
#' par.per.cross <- cbind(unique(cross.ind), rep("B73", 5), rownames(geno.par)[-1])
#' 
#' 
#' \dontrun{
#' 
#' # Specify path to a directory
#' my.dir <- "C:/.../"
#' 
#' # IBS mppData for bi-allelic model
#' data_IBS <- mppData_form(geno.off = geno[7:506,], geno.par = geno.par,
#'                            type = "F", F.gen = 6, IBS = TRUE,
#'                            IBS.format = "ATCG", map = map,
#'                            trait = trait, cross.ind = cross.ind,
#'                            par.per.cross = par.per.cross, dir = my.dir)
#' 
#' }
#' 
#' @export
#'


mppData_form <- function(geno.off, geno.par = NULL, IBS = FALSE,
                         IBS.format = "ATCG",type, F.gen = NULL,
                         BC.gen = NULL, type.mating = NULL, map, trait,
                         cross.ind, par.per.cross, step = 10000,
                         error.prob = 1e-04, map.function = "haldane",
                         stepwidth = "max", dir = getwd()) {
  
  
  # 1. chech of the data format
  #############################
  
  check.mppData(geno = geno.off, geno.par = geno.par, biall = IBS,
                IBS.format = IBS.format,  type = type,
                type.mating = type.mating, BC.gen = BC.gen, F.gen = F.gen,
                trait = trait, map = map, cross.ind = cross.ind,
                par.per.cross = par.per.cross, dir = dir)
  
  
  # 2. Elements that are similar to the two types of mppData object
  #################################################################
  
  ### 2.1 keep the geno names
  
  geno.names <- rownames(geno.off)
  
  ### 2.2 phenotypic values
  
  trait.val <- subset(x = trait, select = 2, drop = FALSE)
  
  # we keep only the phenotypic values
  
  ### 2.3 type of population
  
  if (type == "F") {
    
    type.pop <- paste0("F", "(n = ", F.gen,")")
    
  } else if (type == "BC") {
    
    type.pop <- paste0("Back-cross ", "(n = ", BC.gen,")")
    
  } else if (type == "DH") {
    
    type.pop <- "Double haploid"
    
  } else if (type == "RIL") {
    
    type.pop <- "Recombinant inbred line"
    
    if(type.mating == "selfing") {
      
      type.pop <- paste (type.pop, "by selfing")
      
    } else {
      
      type.pop <- paste (type.pop, "by sibling mating")
      
    }
    
  } else if (type == 'BCsFt'){
    
    type.pop <- paste0('Back-cross followed by selfing ', '(',
                       paste0('BC', BC.gen, 'F', F.gen), ')')
    
  }
  
  ### 2.4 number of allele class
  
  if ((type == "BC") | (type == "DH") | (type == "RIL")) {
    
    n.zigo <- 2
    
  } else if ((type == "F")| (type == "BCsFt")) {
    
    n.zigo <- 3
    
  }
  
  ### 2.5 pedigree data.frame
  
  p1 <- par.per.cross[, 2]
  p2 <- par.per.cross[, 3]
  
  names(p2) <- names(p1) <- par.per.cross[, 1]
  
  ped.mat <- data.frame(rep("offspring", length(geno.names)), trait[, 1],
                        p1[cross.ind], p2[cross.ind], stringsAsFactors = FALSE)
  colnames(ped.mat) <- c("type" ,"genotypes", "parent1", "parent2")
  
  ### 2.6 list of parents number of parents and cross
  
  parents <- union(par.per.cross[, 2], par.per.cross[, 3])
  
  n.par <- length(parents)
  
  n.cr <- dim(par.per.cross)[1]
  
  
  # 3. Elements that are not similar for the different mppData objects
  ####################################################################
  
  ### 3.1 mppData object for cross, parental and ancestral model
  
  if(!IBS){
    
    # 3.1.1 Compute IBD probabilities
    
    # format data to form a cross object
    
    chr.info <- t(map[, 2:3])
    colnames(chr.info) <- map[, 1]
    
    geno.aug <- rbind(chr.info, geno.off)
    
    trait.aug <- c(c("", ""), trait[, 2])
    
    geno.aug <- cbind(trait.aug, geno.aug)
    colnames(geno.aug)[1] <- colnames(trait)[2] 
    
    # Export the data in a .csv file in the specified directory
    
    file.name <- paste(dir, "/", "Cross_object.csv", sep = "")
    
    write.csv(geno.aug, file = file.name, row.names = FALSE)
    
    # form a R/qtl cross object reading the data using the specifiec type of
    # population
    
    if (type == "F") {
      
      cross.object <- read.cross("csv", , file.name, F.gen = F.gen, 
                                 crosstype = "bcsft")
      
      
    } else if (type == "BC") {
      
      cross.object <- read.cross("csv", , file.name, BC.gen = BC.gen, 
                                 crosstype = "bcsft")
      
      
    } else if (type == "RIL") {
      
      # need to read the object as a backcross
      
      cross.object <- read.cross("csv", file = file.name, genotypes = c("A", "B"),
                                 alleles = c("A", "B"))
      
      # then convert it following the type of mating
      
      if (type.mating == "selfing") {
        
        cross.object <- convert2riself(cross.object)
        
      }
      
      if (type.mating == "sib.mat") {
        
        cross.object <- convert2risib(cross.object)
        
      }
      
    } else if (type == "DH") {
      
      # need to read the object as a backcross
      
      cross.object <- read.cross("csv", file = file.name, genotypes = c("A", "B"),
                                 alleles = c("A", "B"))
      
      class(cross.object)[1] <- "dh"
      
    } else if (type == "BCsFt"){
      
      cross.object <- read.cross("csv", , file.name, F.gen = F.gen,
                                 BC.gen = BC.gen, crosstype = "bcsft")
      
    }
    
    # computation of the IBD probabilities
    
    cross.object <- calc.genoprob(cross.object, step = step,
                                  error.prob = error.prob, 
                                  stepwidth = stepwidth,
                                  map.function = map.function)
    
    # modify the names of the added position adding a chromosome indicator
    # (1_, 2_, etc.), include genotype indicators and make a new map
    
    mark.names <- c()
    positions <- c()
    chr.ind <- c()
    
    for (i in 1:nchr(cross.object)) {
      
      names <- dimnames(cross.object$geno[[i]]$prob)[[2]] 
      
      sel.mark <- which(names %in% map[, 1] == FALSE)
      
      names[sel.mark] <- paste(i, names[sel.mark], sep = "_")
      
      dimnames(cross.object$geno[[i]]$prob)[[1]] <- geno.names
      
      dimnames(cross.object$geno[[i]]$prob)[[2]] <- names
      
      names(attributes(cross.object$geno[[i]]$prob)$map) <- names
      
      mark.names <- c(mark.names, names)
      
      positions <- c(positions, attributes(cross.object$geno[[i]]$prob)$map)
      
      chr.ind <- c(chr.ind, rep(i, length(names)))
      
    }
    
    # 3.1.2 map
    
    new.map <- data.frame(mark.names, chr.ind, sequence(table(chr.ind)),
                          positions, stringsAsFactors = FALSE)
    colnames(new.map) <- c("mk.names","chr","pos.ind","pos.cM")
    
    
    # 3.1.3 parent genotypes (if provided)  
    
    if(!is.null(geno.par)){
      
      geno.par <- geno.par[parents, ] # re-order the parents rows
      
      #  if positions were added fill them in the parent genotypes
      
      if (!identical(new.map[, 1], colnames(geno.par))) {
        
        match.id <- match(new.map[, 1], colnames(geno.par))
        
        geno.par.new <- c()
        
        for (i in 1:length(new.map[, 1])) {
          
          if (is.na(match.id[i])) {
            
            geno.par.new <- rbind(geno.par.new, rep("-", n.par))
            
          } else {
            
            geno.par.new <- rbind(geno.par.new,
                                  as.character(geno.par[, match.id[i]]))
            
          }
          
        }
        
      } else { geno.par.new <- t(geno.par) } # no extra position. geno.par stay the same
      
      # combine with the new map information
      
      geno.par.new <- data.frame(new.map, geno.par.new, stringsAsFactors = FALSE)
      colnames(geno.par.new)[5:dim(geno.par.new)[2]] <- parents
      
      mppData <- list(geno = cross.object, geno.id = geno.names,
                      geno.par = geno.par.new, map = new.map, trait = trait.val,
                      cross.ind = cross.ind, ped.mat = ped.mat,
                      par.per.cross = par.per.cross, parents = parents,
                      n.cr = n.cr, n.par = n.par, type = type.pop, n.zigo =
                        n.zigo, biall = IBS)
      
      class(mppData) <- c("mppData", "list")
      
      return(mppData)  
      
    } else {
      
      mppData <- list(geno = cross.object, geno.id = geno.names, geno.par = NULL,
                      map = new.map,
                      trait = trait.val, cross.ind = cross.ind, ped.mat = ped.mat,
                      par.per.cross = par.per.cross, parents = parents,
                      n.cr = n.cr, n.par = n.par, type = type.pop, n.zigo = n.zigo,
                      biall = IBS)
      
      class(mppData) <- c("mppData", "list")
      
      return(mppData)
      
    } 
    
    ### 3.2 mppData object for bi-allelic model      
    
  } else {
    
    # 3.2.1 Transform marker score at 0, 1, 2 format (0 = AA, 1 = AT, 2 = TT)
    
    
    if(IBS.format == "ATCG"){
      
      geno.trans <- geno_012(mk.mat = geno.off)
      
      geno012 <- geno.trans[[1]]
      
      allele.ref <- geno.trans[[2]]
      
    } else if (IBS.format == "012") {
      
      geno012 <- geno.off
      
      allele.ref <- NULL
      
      geno.par <- NULL
      
    }
    
    
    # 3.2.2 map
    
    pos.ind <- sequence(table(map[, 2]))
    
    new.map <- data.frame(map[, 1], map[, 2], pos.ind, map[, 3],
                          stringsAsFactors = FALSE)
    
    colnames(new.map) <- c("mk.names", "chr", "pos.ind", "pos.cM")
    
    # 3.2.3 parent genotypes (if provided)  
    
    if (!is.null(geno.par)) {
      
      geno.par <- geno.par[parents, ]
      
      geno.par <- data.frame(new.map, t(geno.par), stringsAsFactors = FALSE)
      colnames(geno.par)[5:dim(geno.par)[2]] <- parents
      
      mppData <- list(geno = geno012, allele.ref = allele.ref,
                      geno.id = geno.names, geno.par = geno.par,
                      map = new.map, trait = trait.val, cross.ind = cross.ind,
                      ped.mat = ped.mat, par.per.cross = par.per.cross,
                      parents = parents, n.cr = n.cr, n.par = n.par,
                      type = type.pop, n.zigo = n.zigo, biall = IBS)
      
      class(mppData) <- c("mppData", "list")
      
      return(mppData)
      
    } else {
      
      mppData <- list(geno = geno012, allele.ref = allele.ref,
                      geno.id = geno.names, geno.par = NULL, map = new.map,
                      trait = trait.val, cross.ind = cross.ind, ped.mat = ped.mat,
                      par.per.cross = par.per.cross, parents = parents,
                      n.cr = n.cr, n.par = n.par, type = type.pop,
                      n.zigo = n.zigo, biall = IBS)
      
      class(mppData) <- c("mppData", "list")
      
      return(mppData)
      
    }
    
  }
  
}
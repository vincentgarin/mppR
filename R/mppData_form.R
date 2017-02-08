################
# mppData_form #
################

#' Data for MPP QTL analysis
#' 
#' Form a single data object containing (almost) all data necessary for the MPP
#' QTL anlyses.
#' 
#' \code{mppData_form} is a function that gather into a single data object the
#' necessary data for QTL analyses. Two type of mppData object can be produced.
#' The first one, for cross-specific, parental and ancestral models
#' (for details see \code{\link{mpp_SIM}}), calculates identical by descent
#' (IBD) probabilities using
#' \code{read.cross} and \code{calc.genoprob} functions from the R/qtl package
#' (Broman et al. 2009). The results are saved in a cross object. Several types
#' of populations are allowed: F-type (F), back-cross (bc), double
#' haploid (dh) and recombinant imbred lines (RIL). DH and RIL populations are
#' read as back-cross by R/qtl. For these two population types,
#' no heterozygosity is allowed. It must be converted into missing.
#' \code{\link{USNAM_mppData}} is an example of an IBD \code{mppData} object.
#' 
#' The second type of mppData object for identical by state (IBS) genetic
#' predictors for bi-allelic model. In such case
#' genotype marker score as simply translated into 0, 1, 2 format representing
#' the copy number of the less frequent allele. An homozygous marker score
#' with two alleles having the highest frequency is therefore set as reference
#' (0).
#' \code{\link{USNAM_mppData_bi}} is an example of an IBS \code{mppData} object.
#' 
#' @param geno.off \code{Character matrix} representing genotype marker scores
#' with genotypes as row and markers as column. For cross-specific, parental and
#' ancestral models (\code{biall = FALSE}) marker matrix in ABH format.
#' \strong{The ABH assignement must be done per cross. "A"("B") score means that
#' at the considered position the genotype received its allele from parent 1 (2)
#' of the cross. "H" means that the marker is heterozygous and NA or "-" that it
#' is missing or undertermined. This can be done using \code{\link{cross_ABH}}.}
#' 
#' For bi-allelic models (\code{biall = TRUE}) \strong{marker scores must be
#' coded using one letter for each allele. For example, AA, CC, GG, TT, AC,
#' AG, AT, CA, CG, CT, GA, GC, GT, TA, TC, TG. Missing values must be coded NA.}
#' 
#' For both type of mppData object (bi-allelic or not), \strong{the row names
#' (rownames(geno.off)) must be the genotypes identifiers similars to the one of
#' trait argument. The column names (colnames(geno.off)) must be the marker
#' identifiers similar to the first column of the map.}
#' 
#' @param geno.par \code{Character matrix} representing marker scores of the
#' parents with genotypes as row and markers as column.\strong{ The column
#' names must be the as the marker list in map and geno.off arguments.
#' The rownames must represent the parents identifiers and be identical to
#' the parent list
#' given in par.per.cross argument}. This argument is optional and is only
#' use as complementary information for the results. It is however
#' strongly recommanded to include it. Default = NULL.
#' 
#' @param biall \code{Logical} value. Put \code{biall = TRUE} if you want to
#' obtain an IBS \code{mppData} object for the bi-allelic model.
#' Default = FALSE.
#' 
#' @param type \code{Character} indicator for the type of population analysed.
#' type = "F" for Fn (F cross n gernerations), type = "bc" for BCn (backcross
#' n generations).  Use nb.gen argument to secify the number of generations
#' of selfing or back-crossing. Type = "dh" for double haploids and type = "RIL"
#' for recombinant inbred lines. For RIL type specify if the population was
#' obtain using selfing or sibling mating using type.mating argument.
#' If type = "RIL" or "dh" heterozygous marker score of the geno.off argument must
#' be set as missing ("-" or NA).
#' 
#' @param nb.gen \code{Numeric} integer representing the number of generations
#' for F population and backcross. For example nb.gen = 2 for F2 or nb.gen = 1 for
#' single backcross BC. Default = NULL.
#' 
#' @param type.mating Character specifying the type of mating for a RIL population.
#' If type.mating = "selfing" selfing, if type.mating = "sib.mat" sibling mating.
#' Default = NULL.
#' 
#' @param map Three columns \code{data.frame} with: 1) marker or in between
#' position identifiers; 2) chromosome; 3) positions in centi-Morgan.
#' The marker identifiers must be identical to the column names of the maker
#' matrix (argument geno.off).
#' 
#' @param trait two columns \code{data.frame} with : 1) \code{character}
#' genotypes identifiers; 2) \code{numeric} trait values. \strong{The genotypes
#' identifiers must be identical to the rownames of the  marker matrix
#' (argument \code{geno.off}).}
#' 
#' @param cross.ind \code{Character} vector with the same length as the number of
#' genotypes which specifies to which cross each genotype belong.
#' 
#' @param par.per.cross Three columns \code{character matrix} specifying :
#' 1) the cross indicators (\strong{In the same order as they appear in
#' \code{cross.ind}}); 2) the parents 1 identifiers of the crosses;
#' 3) the parents 2 identifiers of the crosses. \strong{The alleles coming from
#' parent 1 (2) within each cross correspond to the scores A (B) in argument
#' \code{geno.off} for a non bi-allelic model. The cross indicators must be 
#' similar to the one used in the cross.ind argument. The list of parent
#' identifiers must be the same to the rownames of par.sc argument.}
#' 
#' @param step \code{Numeric} value representing the maximum distance (in cM)
#' between positions at which the genotype probabilities are calculated.
#' Parameter for the function calc.genoprob().Definition of this parameter
#' from R/qtl package (Broman et al. 2009)). Default value = 10000 (the function
#' do not add position in between markers).
#'  
#' @param error.prob \code{Numeric} value for assumed genotyping error rate
#' used in
#' the calculation of the penetrance Pr(observed genotype | true genotype).
#' Parameter for the function \code{calc.genoprob()}. Definition of this parameter
#' from R/qtl package (Broman et al. 2009). Default = 0.0001.
#' 
#' @param map.function \code{Character} expression specifying the type of map
#' function used to infer the IBD probabilities. possibility to choose
#' between "haldane", "kosambi","c-f","morgan".
#' Parameter for the function \code{calc.genoprob()}. Default = "haldane".
#' 
#' @param stepwidth \code{Character} expression indicating if the intermediate
#' points should be fixed ("fixed") or if the algorithm should insert a minimal
#' number of points in between positions ("max"). If stepwidth = "max" the
#' maximum distance between points is step. Stepwidth correspond to the
#' parameter passed to the function \code{calc.genoprob()}.
#' Default = "max".
#' 
#' @param dir Path where a .csv file Cross_object.csv will be stored to be
#' used for cross object formation. By default the function uses the current
#' working directory.
#'
#' @return 
#' 
#' Return: a {list} of class \code{mppData} containing the following items :
#' 
#' \item{geno}{For a cross-specific, parental and ancestral model
#' (\code{biall = FALSE}), cross object obtained with R/qtl functions
#' The main elements of the cross object are the IBD probabilities. For a
#' bi-allelic model (\code{biall = TRUE}), marker score matrix recoded
#' as 0, 1, 2 according to the copies number of the most frequent allele.}
#' 
#' \item{allele.ref}{If \code{biall = TRUE}, \code{character matrix} with
#' reference allele scores. The first row represent the allele with the highest
#' MAF, the second, the one with the lowest and the two other lines represent
#' the heterozygous scores.}
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
#' 4) Positions in cM.}
#' 
#' \item{trait}{\code{Data.frame} with phenotypic values.}
#' 
#' \item{cross.ind}{\code{Character} vector specifying to which cross each
#' genotype belong.}
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
#' \item{biall}{\code{Logic} value specifying if data are made for a
#' bi-allelic model}
#' 
#'  
#' @author Vincent Garin
#' 
#' @seealso \code{\link{cross_ABH}}, \code{\link{mppData_subset}},
#' \code{\link{mppData_chgPheno}}, \code{\link{mppData_chgPedigree}},
#' \code{\link{mpp_SIM}},
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
#' @examples
#'  
#' # Mpp data for cross parental or ancestral model
#' ################################################
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
#' trait <- USNAM_pheno[, 1:2]
#' colnames(trait) <- c("genotype", "ULA")
#' 
#' # Map
#' data(USNAM_map)
#' 
#' # Same list of markers as in the genotype matrix
#' map <- USNAM_map[USNAM_map[, 1] %in% colnames(USNAM_genoABH),]
#' 
#' cross.ind <- USNAM_pheno[, 3]
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
#'                      biall = FALSE, type = "F", nb.gen = 6, map = map,
#'                      trait = trait, cross.ind = cross.ind,
#'                      par.per.cross = par.per.cross, step = 5,
#'                      map.function = "haldane",  stepwidth = "max",
#'                      dir = my.dir)
#' }
#'
#' 
#' # Mpp data for bi-allelic model
#' ##############################
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
#' trait <- USNAM_pheno[, 1:2]
#' colnames(trait) <- c("genotype", "ULA")
#' 
#' # Map
#' data(USNAM_map)
#' 
#' # Same list of markers as in the genotype matrix
#' map <- USNAM_map[USNAM_map[, 1] %in% colnames(geno),]
#' 
#' cross.ind <- USNAM_pheno[, 3]
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
#' data_biall <- mppData_form(geno.off = geno[7:506,], geno.par = geno.par,
#'                            type = "F", nb.gen = 6, biall = TRUE, map = map,
#'                            trait = trait, cross.ind = cross.ind,
#'                            par.per.cross = par.per.cross, dir = my.dir)
#' 
#' }
#' 
#' @export
#'


mppData_form <- function(geno.off, geno.par = NULL, biall = FALSE, type,
                         nb.gen = NULL, type.mating = NULL, map, trait,
                         cross.ind, par.per.cross, step = 10000,
                         error.prob = 1e-04, map.function = "haldane",
                         stepwidth = "max", dir = getwd()) {
  
  
  # 1. chech of the data format
  #############################
  
  check.mppData(geno = geno.off, geno.par = geno.par, biall = biall, type = type,
                type.mating = type.mating, nb.gen = nb.gen, trait = trait,
                map = map, cross.ind = cross.ind, par.per.cross = par.per.cross,
                dir = dir)
  
  
  # 2. Elements that are similar to the two types of mppData object
  #################################################################
  
  ### 2.1 keep the geno names
  
  geno.names <- rownames(geno.off)
  
  ### 2.2 phenotypic values
  
  trait.val <- subset(x = trait, select = 2, drop = FALSE)
  
  # we keep only the phenotypic values
  
  ### 2.3 type of population
  
  if (type == "F") {
    
    type.pop <- paste0("F", "(n = ", nb.gen,")")
    
  } else if (type == "bc") {
    
    type.pop <- paste0("Back-cross ", "(n = ", nb.gen,")")
    
  } else if (type == "dh") {
    
    type.pop <- "Double haploid"
    
  } else if (type == "RIL") {
    
    type.pop <- "Recombinant inbred line"
    
    if(type.mating == "selfing") {
      
      type.pop <- paste (type.pop, "by selfing")
      
    } else {
      
      type.pop <- paste (type.pop, "by sibling mating")
      
    }
    
  }
  
  ### 2.4 number of allele class
  
  if ((type == "bc") | (type == "dh") | (type == "RIL")) {
    
    n.zigo <- 2
    
  } else if ((type == "F")) {
    
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
  
  if(!biall){
    
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
      
      cross.object <- read.cross("csv", , file.name, F.gen = nb.gen, 
                                 crosstype = "bcsft")
      
      
    } else if (type == "bc") {
      
      cross.object <- read.cross("csv", , file.name, BC.gen = nb.gen, 
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
      
    } else if (type == "dh") {
      
      # need to read the object as a backcross
      
      cross.object <- read.cross("csv", file = file.name, genotypes = c("A", "B"),
                                 alleles = c("A", "B"))
      
      class(cross.object)[1] <- "dh"
      
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
                      n.zigo, biall = biall)
      
    class(mppData) <- c("mppData", "list")
  
    return(mppData)  
    
  } else {
    
    mppData <- list(geno = cross.object, geno.id = geno.names, map = new.map,
                    trait = trait.val, cross.ind = cross.ind, ped.mat = ped.mat,
                    par.per.cross = par.per.cross, parents = parents,
                    n.cr = n.cr, n.par = n.par, type = type.pop, n.zigo = n.zigo,
                    biall = biall)
    
    class(mppData) <- c("mppData", "list")
    
    return(mppData)
    
  } 
      
  ### 3.2 mppData object for bi-allelic model      
    
  } else {
    
  # 3.2.1 Transform marker score at 0, 1, 2 format (0 = AA, 1 = AT, 2 = TT)

    geno.trans <- geno_012(mk.mat = geno.off)
    
    geno012 <- geno.trans[[1]]
    
    allele.ref <- geno.trans[[2]]
    
  # 3.2.2 map
    
    pos.ind <- sequence(table(map[, 2]))

    new.map <- data.frame(map[, 1], map[, 2], pos.ind, map[, 3],
                          stringsAsFactors = FALSE)

    colnames(new.map) <- c("mk.names", "chr", "pos.ind", "pos.cM")
  
  # 3.2.3 parent genotypes (if provided)  
  
    if (!is.null(geno.par)) {
    
      geno.par <- geno.par[parents, ]
      
      geno.par <- data.frame(new.map, t(geno.par), stringsAsFactors = FALSE)
      
      mppData <- list(geno = geno012, allele.ref = allele.ref,
                      geno.id = geno.names, geno.par = geno.par,
                      map = new.map, trait = trait.val, cross.ind = cross.ind,
                      ped.mat = ped.mat, par.per.cross = par.per.cross,
                      parents = parents, n.cr = n.cr, n.par = n.par,
                      type = type.pop, n.zigo = n.zigo, biall = biall)
      
      class(mppData) <- c("mppData", "list")
      
      return(mppData)
        
    } else {
      
      mppData <- list(geno = geno012, allele.ref = allele.ref,
                      geno.id = geno.names, map = new.map,
                      trait = trait.val, cross.ind = cross.ind, ped.mat = ped.mat,
                      par.per.cross = par.per.cross, parents = parents,
                      n.cr = n.cr, n.par = n.par, type = type.pop,
                      n.zigo = n.zigo, biall = biall)
      
      class(mppData) <- c("mppData", "list")
      
      return(mppData)
      
    }
  
  }
    
} # End function
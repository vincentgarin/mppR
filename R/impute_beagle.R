########################################
# Imputation genetic data using Beagle #
########################################

#' Impute maker data using Beagle
#'
#' Impute maker data using Beagle (Browning & Browning, 2018).
#' 
#' @param mk.mat \code{Character} genotype  marker score \code{matrix}.
#' \strong{Marker scores must be coded using one letter for each allele.
#' For example, AA, CC, GG, TT, AC, AG, AT, CA, CG, CT, GA, GC, GT, TA, TC, TG.
#' Missing values must be coded NA.}
#' 
#' @param map Four columns \code{data.frame} with: 1) marker names, 2) chromosome,
#' 3) position in cM, 4) physical bp position (must be in ascending order).
#' 
#' @param alleles 4xN_marker matrix with SNP marker alleles of dimensions
#' 4 x N_marker: homozygous score minor allele (first row), 
#' homozygous score major allele (second row), heterozygous score 1, heterozygous
#' score 2.
#' 
#' @param beagle.loc \code{Character} path specifying the location where the
#' Beagle software is located.
#' 
#' @param undiferentiate.het.score \code{Logical} value specifying if in the imputed
#' value the two different phased heterozygous scores (e.g. AT and TA) should be
#' considered as similar (e.g. AT). Default = TRUE.
#' 
#' @return
#' 
#' imputed genotype marke matrix
#' 
#' @author Vincent Garin
#' 
#' @references 
#' 
#' B L Browning, Y Zhou, and S R Browning (2018). A one-penny imputed genome
#' from next generation reference panels. Am J Hum Genet 103(3):338-348.
#' doi:10.1016/j.ajhg.2018.07.015
#' 
#' @export
#' 

impute_beagle <- function(mk.mat, map, alleles, beagle.loc,
                          undiferentiate.het.score = TRUE){
  
# control that the list of markers is the same in the map and in the mk.mat

if(!identical(map[, 1], colnames(mk.mat))){
    
  stop('The list of marker in the map and in mk.mat are not identical')
    
}
  
  if(is.unsorted(map[, 4], strictly = TRUE)){
    
    stop('The map physical position are not strictly ascending.')
    
  }

# write input file - vcf gt data fix part
#########################################  
  
# start from an existing object and try to insert my values
  
# read vcf object
# vcf_obj <- read.vcfR(file = 'G:/Beagle/ref.18May20.d20.vcf.gz')
  
data(vcfR_example)

# prepare fixed information

mk.mat <- t(mk.mat)

ref <- substring(alleles[2, ], 1, 1)
alt <- substring(alleles[1, ], 1, 1)
homo_ref <- paste0(ref, ref)
homo_alt <- paste0(alt, alt)
het1 <- paste0(ref, alt)
het2 <- paste0(alt, ref)

d_fix <- data.frame(chr = map[, 2], pos = map[, 4], id = map[, 1],
                    ref = ref, alt = alt, qual = '.', filter = 'PASS',
                    info = '.')

d_fix$chr <- as.character(d_fix$chr)
d_fix$pos <- as.character(d_fix$pos)
d_fix <- as.matrix(d_fix)
colnames(d_fix) <- colnames(vcf@fix)
rownames(d_fix) <- NULL

# change in the vcf object
vcf@fix <- d_fix

##########

# write input file - vcf gt data marker scores part
###################################################
  
# convert the marker scores in 0/0, 0/1 and 1/1 format

mk.mat.new <- matrix('./.', nrow = dim(mk.mat)[1], ncol = dim(mk.mat)[2])

for(i in 1:dim(mk.mat)[1]){
  
  mk.mat.new[i, mk.mat[i, ] == homo_ref[i]] <- '0/0'
  mk.mat.new[i, mk.mat[i, ] == homo_alt[i]] <- '1/1'
  mk.mat.new[i, mk.mat[i, ] %in% c(het1[i], het2[i])] <- '1/0'
  
  
}

colnames(mk.mat.new) <- colnames(mk.mat)

vcf_gt_mine <- mk.mat.new

n_0 <- max(nchar(1:dim(mk.mat)[2])) + 1
s_nb <- str_pad(1:dim(mk.mat)[2], n_0, pad = "0")
sample_id <- paste0('SAMP', s_nb)

colnames(mk.mat.new) <- sample_id

vcf_gt <- cbind('GT', mk.mat.new)
colnames(vcf_gt)[1] <- 'FORMAT'

# change in the vcf object
vcf@gt <- vcf_gt

########

# write vcf object
write.vcf(x = vcf, file = file.path(beagle.loc, 'geno.vcf.gz'))

# write input file - plink map information
##########################################

# map_plk <- data.frame(map[, 2], map[, 1], map[, 3], map[, 4],
#                        stringsAsFactors = FALSE)

map_plk <- data.frame(map[, 2], map[, 1], map[, 3], map[, 4],
                      stringsAsFactors = FALSE)

gt_file <- file.path(beagle.loc, 'plink_map.map.gz')
gz1 <- gzfile(gt_file, "w")

write.table(map_plk, file = gz1, row.names = FALSE, col.names = FALSE,
            quote = FALSE)

close(gz1)


### 

# run Beagle from terminal

setwd(beagle.loc)

cmd <- 'java -jar beagle.18May20.d20.jar gt=geno.vcf.gz map=plink_map.map.gz out=out.gt'

system(command = 'beagle_execute.sh')

system(cmd)

# get the output

vcf_out <- read.vcfR(file = file.path(beagle.loc ,'out.gt.vcf.gz'))

# remove the created files

file.remove(c('geno.vcf.gz', 'plink_map.map.gz', 'out.gt.vcf.gz', 'out.gt.log'))

# get the gt object with the allele

gt <- extract.gt(vcf_out, return.alleles = TRUE)

colnames(gt) <- colnames(mk.mat)

# set the marker into AA, TT, GT, etc format

gt <- Beagle_mk_code_convert(gt, undiferentiate.het.score)

return(gt)

}
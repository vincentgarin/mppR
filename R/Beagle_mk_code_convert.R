##########################
# Beagle_mk_code_convert #
##########################

Beagle_mk_code_convert <- function(data, undiferentiate.het.score){
  
  if(undiferentiate.het.score){
    
    data[data == 'A|A'] <- 'AA'
    data[data == 'C|C'] <- 'CC'
    data[data == 'T|T'] <- 'TT'
    data[data == 'G|G'] <- 'GG'
    
    data[data == 'A|C'] <- 'AC'
    data[data == 'C|A'] <- 'AC'
    
    data[data == 'A|G'] <- 'AG'
    data[data == 'G|A'] <- 'AG'
    
    data[data == 'C|T'] <- 'CT'
    data[data == 'T|C'] <- 'CT'
    
    data[data == 'C|G'] <- 'CG'
    data[data == 'G|C'] <- 'CG'
    
    data[data == 'A|T'] <- 'AT'
    data[data == 'T|A'] <- 'AT'
    
    data[data == 'G|T'] <- 'GT'
    data[data == 'T|G'] <- 'GT'
    
    
  } else {
    
    data[data == 'A|A'] <- 'AA'
    data[data == 'C|C'] <- 'CC'
    data[data == 'T|T'] <- 'TT'
    data[data == 'G|G'] <- 'GG'
    
    data[data == 'A|C'] <- 'AC'
    data[data == 'C|A'] <- 'CA'
    
    data[data == 'A|G'] <- 'AG'
    data[data == 'G|A'] <- 'GA'
    
    data[data == 'C|T'] <- 'CT'
    data[data == 'T|C'] <- 'TC'
    
    data[data == 'C|G'] <- 'CG'
    data[data == 'G|C'] <- 'GC'
    
    data[data == 'A|T'] <- 'AT'
    data[data == 'T|A'] <- 'TA'
    
    data[data == 'G|T'] <- 'GT'
    data[data == 'T|G'] <- 'TG'
    
  }
  
  
  return(data)
  
}
#
# rqtl.to.qtab.utils.R
#
# Copyright (c) 2011 Danny Arends, Pjotr Prins, Karl Broman and Ritsert C. Jansen
# last modified Jan, 2012
# first written Jan, 2012
#
# rqtl.symbol.toN
# rqtl.symbol.toTG

get.chr <- function(cross,m=5){
  chr <- 1
  nmar <- length(pull.map(cross)[[chr]])
  while(m > nmar){
    m <- m - nmar
    chr <- chr + 1
    nmar <- length(pull.map(cross)[[chr]])
  }
  chr;
}

get.loc <- function(cross,m=1){
  unlist(pull.map(cross))[m]
}

get.phenotype.type <- function(cross, phenotype){
 if(class(cross$pheno[,phenotype])=="numeric") return("Float")
 if(class(cross$pheno[,phenotype])=="character") return("Char")
 if(class(cross$pheno[,phenotype])=="factor") return("Int")
 return("Float")
}

# end of rqtl.to.qtab.utils.R

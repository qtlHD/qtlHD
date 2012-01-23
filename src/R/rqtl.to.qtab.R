#
# rqtl.to.qtab.R
#
# Copyright (c) 2011 Danny Arends, Pjotr Prins, Karl Broman and Ritsert C. Jansen
# last modified Jan, 2012
# first written Jan, 2012
#
# getID
# rqtl.to.qtab.symbols, rqtl.to.qtab.location, rqtl.to.qtab.genotypes, rqtl.to.qtab.phenotypes
# rqtl.to.qtab, rqtl.to.qtab.test

getID <- function(){
  VER <- "0.1"
  ID <- paste("qtlHD-in-", VER, sep="")
  ID
}

rqtl.to.qtab.symbols <- function(cross, filename="symbols.qtab",descr = "My cross from r/qtl"){
  cat(file=filename, "# --- ",getID()," Symbol ",descr, "\n",sep="");
  cat(file=filename, "# --- Symbol Genotype begin", "\n",sep="", append=TRUE)
  mysymbols <- unique(as.numeric((pull.geno(cross))))
  for(symbol in mysymbols) {
    cat(file=filename, rqtl.symbol.toN(cross, symbol)," as ", sep="", append=TRUE)
    cat(file=filename, rqtl.symbol.toTG(cross, symbol), "\n", sep="", append=TRUE)
  }
  cat(file=filename, "# --- Symbol Genotype end", "\n", sep="", append=TRUE)
}

rqtl.to.qtab.location <- function(cross, filename="locations.qtab", descr = "My cross from r/qtl"){
  cat(file=filename, "# --- ",getID()," Location ", descr, "\n", sep="")
  cat(file=filename, "# --- Data Location begin", "\n", sep="", append=TRUE)
  genotypes <- pull.geno(cross)
  cnt <- 1
  for(n in colnames(genotypes)){
    cat(file=filename, n,"\t", get.chr(cross, cnt),"\t", get.loc(cross, cnt), "\n", sep="", append=TRUE)
    cnt <- cnt+1
  }
  cat(file=filename, "# --- Data Location end", "\n", sep="", append=TRUE)
}

rqtl.to.qtab.genotypes <- function(cross, filename="genotypes.qtab", descr = "My cross from r/qtl"){
  cat(file=filename, "# --- ",getID()," Genotype ", descr, "\n", sep="")
  cat(file=filename, "# --- Data Genotype begin", "\n", sep="", append=TRUE)
  genotypes <- pull.geno(cross)
  cnt <- 1
  for(n in colnames(genotypes)){
    cat(file=filename, n, "\n", sep="", append=TRUE)
    cnt <- cnt+1
  }
  for(i in 1:nrow(genotypes)){
    cat(file=filename, "Ind", i, sep="", append=TRUE)
    for(p in 1:ncol(genotypes)){
      cat(file=filename, "\t", rqtl.symbol.toTG(cross,genotypes[i,p]), sep="", append=TRUE)
    }
    cat(file=filename, "\n", sep="", append=TRUE)
  }
  cat(file=filename, "# --- Data Genotype end", "\n", sep="", append=TRUE)
}

rqtl.to.qtab.phenotypes <- function(cross, filename="phenotypes.qtab", descr = "My cross from r/qtl"){
  cat(file=filename, "# --- ",getID()," Phenotype ", descr, "\n", sep="")
  cat(file=filename, "# --- Type Phenotype begin", "\n", sep="", append=TRUE)
  for(n in colnames(cross$pheno)) {
    cat(file=filename, n,"\t", get.phenotype.type(cross,n), "\n", sep="", append=TRUE)
  }
  cat(file=filename, "# --- Type Phenotype end", "\n", sep="", append=TRUE)

  cat(file=filename, "# --- Data Phenotype begin", "\n", sep="", append=TRUE)
  phenotypes <- pull.pheno(cross)
  for(n in colnames(phenotypes)){
    cat(file=filename, "#\t", n, "\n", sep="", append=TRUE)
  }
  for(i in 1:nrow(phenotypes)){
    cat(file=filename, "Ind", i, sep="", append=TRUE)
    for(p in 1:ncol(phenotypes)){
      cat(file=filename, "\t", cross$pheno[i,p], sep="", append=TRUE)
    }
    cat(file=filename, "\n", sep="", append=TRUE)
  }
  cat(file=filename, "# --- Data Phenotype end", "\n", sep="", append=TRUE)
}

rqtl.to.qtab <- function(cross, filestem="cross", descr = "My cross from r/qtl", verbose = TRUE){
  if(verbose) cat("Writing symbols\n")
  rqtl.to.qtab.symbols(cross,paste(filestem,"_symbols.qtab",sep=""),descr=descr)
  if(verbose) cat("Writing genetic map\n")
  rqtl.to.qtab.location(cross,paste(filestem,"_location.qtab",sep=""),descr=descr)
  if(verbose) cat("Writing genotypes\n")
  rqtl.to.qtab.genotypes(cross,paste(filestem,"_genotypes.qtab",sep=""),descr=descr)
  if(verbose) cat("Writing phenotypes\n")
  rqtl.to.qtab.phenotypes(cross,paste(filestem,"_phenotypes.qtab",sep=""),descr=descr)
}

rqtl.to.qtab.test <- function(){
  library(qtl)
  setwd("e:/github/qtlHD/src/R/")
  source("rqtl.to.qtab.symbols.R")
  source("rqtl.to.qtab.utils.R")
  data(multitrait)
  data(listeria)
  data(hyper)
  setwd("e:/qtab")
  rqtl.to.qtab(multitrait,"multi")
  rqtl.to.qtab(listeria,"listeria")
  rqtl.to.qtab(hyper,"hyper")
}

# end of rqtl.to.qtab.R

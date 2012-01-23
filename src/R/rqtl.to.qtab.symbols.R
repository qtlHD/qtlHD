#
# rqtl.to.qtab.symbols.R
#
# Copyright (c) 2011 Danny Arends, Pjotr Prins, Karl Broman and Ritsert C. Jansen
# last modified Jan, 2012
# first written Jan, 2012
#
# rqtl.symbol.toN
# rqtl.symbol.toTG

rqtl.symbol.toN <- function(cross,symbol){
  if(is.na(symbol)) return("NA -")
  if(class(cross)[1]=="bc"){
    if(symbol==1) return("A AA")
    if(symbol==2) return("H AB")
  }
  if(class(cross)[1]=="riself" || class(cross)[1]=="risib"){
    if(symbol==1) return("A AA")
    if(symbol==2) return("B BB")
  }
  if(class(cross)[1]=="f2"){
    if(symbol==1) return("A AA")
    if(symbol==2) return("H")
    if(symbol==3) return("B BB")
    if(symbol==4) return("HorA D")
    if(symbol==5) return("HorB C")
  }
}

rqtl.symbol.toTG <- function(cross,symbol){
  if(is.na(symbol)) return("None")
  if(class(cross)[1]=="bc"){
    if(symbol==1) return("0,0")
    if(symbol==2) return("0,1")
  }
  if(class(cross)[1]=="riself" || class(cross)[1]=="risib"){
    if(symbol==1) return("0,0")
    if(symbol==2) return("1,1")
  }
  if(class(cross)[1]=="f2"){
    if(symbol==1) return("0,0")
    if(symbol==2) return("0,1")
    if(symbol==3) return("1,1")
    if(symbol==4) return("0,0 0,1")
    if(symbol==5) return("0,1 1,1")
  }
}

# end of rqtl.to.qtab.symbols.R

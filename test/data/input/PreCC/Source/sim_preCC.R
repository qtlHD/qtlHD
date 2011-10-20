library(qtl)
library(xoi)

#create.parent

set.seed(28280951)

# read and reorganize map
map <- read.csv("../map.csv")
map$pos <- map$pos/2e6
chr <- c(18:19, "X")
nmap <- split(map$pos, map$chr)[chr]
for(i in chr) {
  names(nmap[[i]]) <- map[map$chr==i,1]
  class(nmap[[i]]) <- ifelse(i=="X", "X", "A")
}
class(nmap) <- c("map", "list")
map <- nmap
rm(nmap)

# read and reorganize parental strains' genotypes
hap <- read.csv("../parents.csv", as.is=TRUE)

  
# max position on each chromosome
L <- sapply(map, max)

# no. lines
nlines <- 150

# generations (uniform 6-18)
ngen <- sample(6:18, nlines, repl=TRUE)

# crosses
crosses <- matrix(nrow=nlines, ncol=8)
for(i in 1:nlines) crosses[i,] <- sample(8)

# sexes
sex <- sample(rep(0:1, nlines/2))

# simulate meioses
gen <- vector("list", length(map))
for(i in seq(along=map)) {
  gen[[i]] <- vector("list", nlines)

  for(j in 1:nlines) {
    # generation G1
    g1 <- vector("list", 4)
    for(k in 1:4) g1[[k]] <- create.parent(L[i], crosses[j,k*2 + (-1:0)])

    # generation G2
    mom <- cross(g1[[1]], g1[[2]], m=10, xchr=(names(map)[i]=="X"),
                 male=FALSE)
    dad <- cross(g1[[3]], g1[[4]], m=10, xchr=(names(map)[i]=="X"),
                 male=TRUE)

    for(k in 1:ngen[j]) {
      kid1 <- cross(mom, dad, m=10, xchr=(names(map)[i]=="X"),
                    male=FALSE)
      kid2 <- cross(mom, dad, m=10, xchr=(names(map)[i]=="X"),
                    male=TRUE)
      mom <- kid1
      dad <- kid2
    }
    if(sex[j]==0) gen[[i]][[j]] <- kid1
    else gen[[i]][[j]] <- kid2
  }
}

# turn into discrete genotypes
margen <- vector("list", length(map))
for(i in seq(along=map)) {
  margen[[i]] <- matrix("", ncol=nlines, nrow=length(map[[i]]))
  thishap <- hap[match(names(map[[i]]), hap[,1]),-1]

  for(j in 1:nlines) {
    mav <- gen[[i]][[j]]$mat
    pav <- gen[[i]][[j]]$pat
    
    marow <- parow <- seq(along=map[[i]])

    for(k in 2:ncol(mav)) 
      marow[map[[i]] > mav[1,k-1] & map[[i]] <= mav[1,k]] <- mav[2,k]
    for(k in 2:ncol(pav))
      parow[map[[i]] > pav[1,k-1] & map[[i]] <= pav[1,k]] <- pav[2,k]
    
    machar <- pachar <- rep("", length(map[[i]]))

    for(k in 1:8) {
      wh <- marow==k
      if(any(wh)) machar[wh] <- thishap[wh,k]

      wh <- parow==k
      if(any(wh)) pachar[wh] <- thishap[wh,k]
    }
    margen[[i]][,j] <- paste(machar, pachar, sep="")
  }
}


allgeno <- margen[[1]]
for(i in 2:3) allgeno <- rbind(allgeno, margen[[i]])

# reorder
reord <- rbind(c("CA", "AC"), c("GA", "AG"),
               c("GC", "CG"), c("TA", "AT"),
               c("TC", "CT"), c("TG", "GT"))

for(k in 1:nrow(reord)) 
  allgeno[allgeno==reord[k,1]] <- reord[k,2]

# add 0.3% random missing data
allgeno[sample(c(FALSE, TRUE), length(allgeno), repl=TRUE, prob=c(0.997, 0.003))] <- "-"

colnames(allgeno) <- sprintf("preCC%03d", 1:150)
write.table(cbind("marker"=hap[,1], allgeno), file="../genotypes.csv", sep=",", row.names=FALSE,
            col.names=TRUE, quote=FALSE)

# cross info / sex
crossinfo <- matrix(as.character(crosses), nrow=nlines)
for(i in 1:8) crossinfo <- gsub(i, LETTERS[i], crossinfo)
crossinfo <- apply(crossinfo, 1, paste, collapse="")

write.table(cbind(id=colnames(allgeno), cross=crossinfo, "no. gen"=ngen, sex=c("Female", "Male")[sex+1]),
            file="../phenotypes.csv", sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)

library(qtl)

set.seed(99274828)

# read and reorganize map
map <- read.csv("../map.csv")
map$pos <- map$pos/2e6
chr <- c(1:19, "X")
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
fsnp <- vector("list", 20)
lo <- hi <- rep("", nrow(hap))
for(i in 1:20) {
  fsnp[[i]] <- matrix(0, ncol=8, nrow=length(map[[i]]))
  mn <- rownames(fsnp[[i]]) <- names(map[[i]])
  wh <- match(mn, hap[,1])
  if(any(is.na(wh))) stop("mismatch in marker names on chr ", i)
  for(j in seq(along=wh)) {
    g <- unlist(hap[wh[j],-1])
    tab <- table(g)
    hi[wh[j]] <- gcommon <- names(tab)[tab==max(tab)][1]
    if(length(tab)==1) lo[wh[j]] <- gcommon
    else lo[wh[j]] <- names(tab)[names(tab) != gcommon]
    fsnp[[i]][j,g==gcommon] <- 1
  }
}
  
# simulate CC data
cc <- sim.cross(map, n.ind=250, type="ri8sib", m=10, founderGeno=fsnp, random.cross=TRUE)

# grab the "true" genotypes
gnum <- t(cc$geno[[1]]$truegeno)
for(i in 2:20)
  gnum <- rbind(gnum, t(cc$geno[[i]]$truegeno))

# convert back to A/T/C/G calls
ccg <- matrix("", ncol=ncol(gnum), nrow=nrow(gnum))
for(i in 1:ncol(ccg)) {
  for(k in 1:8) {
    wh <- gnum[,i]==k
    if(any(wh)) ccg[wh,i] <- hap[wh,k+1]
  }
}


# [a more complicated way to do it]
#gbin <- pull.geno(cc) %% 2
#ccgback <- matrix("", ncol=nrow(gbin), nrow=ncol(gbin))
#for(i in 1:nrow(gbin)) {
#  firstpar <- cc$cross[i,1]
#  ishi <- (gbin[i,]==1 & hap[,firstpar+1]==hi) | (gbin[i,]==0 & hap[,firstpar+1]==lo)
#  ccgback[ishi,i] <- hi[ishi]
#  ccgback[!ishi,i] <- lo[!ishi]
#}


# add 0.1% random missing data
ccg[sample(c(FALSE, TRUE), length(ccg), repl=TRUE, prob=c(0.999, 0.001))] <- "-"

# ids
colnames(ccg) <- id <- sprintf("CC%03d", 1:250)

# write genotype data
write.table(cbind("marker"=markernames(cc), ccg), file="../genotypes.csv", sep=",", row.names=FALSE,
            col.names=TRUE, quote=FALSE)

# cross info
crossinfo <- matrix(as.character(cc$cross), nrow=250)
rownames(crossinfo) <- id
for(i in 1:8)
  crossinfo <- gsub(i, LETTERS[i], crossinfo)
write.table(crossinfo, file="../crossinfo.csv", sep=",", row.names=TRUE, col.names=FALSE, quote=FALSE)

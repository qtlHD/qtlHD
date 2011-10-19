library(qtl)

# load original haplotypes
hap <- read.csv("orig_haplotypes.csv", as.is=TRUE)
cn <- scan("orig_haplotypes.csv", what=character(), nlines=1, sep=",")
colnames(hap) <- cn

# write map file
write.table(hap[,1:3], file="../map.csv", sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)

# write haplotype file
write.table(hap[,c(1,4:7)], file="../haplotypes.csv", sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE) 

# create map
map <- split(hap[,3], hap[,2])
chr <- c(1:9, 10:19, "X")
map <- map[chr]
for(i in chr) {
  names(map[[i]]) <- hap[hap[,2]==i,1]
  class(map[[i]]) <- "A"
}
class(map[["X"]]) <- "X"
class(map) <- c("map", "list")

# sex-specific version
mapss <- map
for(i in seq(along=mapss))  {
  mapss[[i]] <- rbind(mapss[[i]]*1.2, mapss[[i]]*0.8)
  class(mapss[[i]]) <- "A"
}
mapss[["X"]] <- rbind(map[["X"]], rep(0, length(map[["X"]])))
class(mapss[["X"]]) <- "X"


# simulate cross
set.seed(72922275)
x <- sim.cross(mapss, n.ind=250, type="4way", model=c(2, 30, 0.2, 0.5, 0.8))

# include some males
x$pheno$sex <- sample(0:1, nind(x), repl=TRUE)
g <- x$geno[["X"]]$data[x$pheno$sex==1,]
g <- g+2
x$geno[["X"]]$data[x$pheno$sex==1,] <- g

# re-estimate map
mapss.est <- est.map(x)

# make tabular
mapss.est <- pull.map(replacemap(x, mapss.est), as.table=TRUE)

# write sex-specific map
write.table(data.frame(marker=rownames(mapss.est), chr=mapss.est[,1], round(mapss.est[,2:3], 4)),
            file="../map_sexsp.csv", sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)

# genotype data: markers by individuals 
g <- t(pull.geno(x))

# mother and father alleles
ma <- pa <- g
for(i in 1:nrow(g)) {
  ma[i,ma[i,]==1 | ma[i,]==3] <- hap[i,4]
  ma[i,ma[i,]==2 | ma[i,]==4] <- hap[i,5]

  pa[i,pa[i,]==1 | pa[i,]==2] <- hap[i,6]
  pa[i,pa[i,]==3 | pa[i,]==4] <- hap[i,7]
}


# function for sorting then combining alleles 
sortpaste <- function(a, b, sep="|") 
{
  z <- rep("", length(a))
  for(i in seq(along=a)) {
    if(a[i] > b[i]) z[i] <- paste(b[i], a[i], sep=sep)
    else z[i] <- paste(a[i], b[i], sep=sep)
  }
  z
}

# genotypes as strings like "192|194"
gstr <- matrix("", ncol=ncol(g), nrow=nrow(g))
for(i in 1:nrow(g)) 
  gstr[i,] <- sortpaste(ma[i,], pa[i,], sep="|")


# simulate missing data
gstr[sample(c(FALSE, TRUE), length(gstr), repl=TRUE, prob=c(0.95, 0.05))] <- "-"


# individual IDs
x$pheno$id <- 1:nind(x)

# genotype data with IDs
gstr <- rbind(id=1:nind(x), gstr)
rownames(gstr)[-1] <- markernames(x)

# write phenotype data
write.table(round(x$pheno, 4), file="../phenotypes.csv", quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)

# write genotype data
write.table(gstr, file="../genotypes_mar-by-ind.csv", quote=FALSE, sep=",", row.names=TRUE, col.names=FALSE)

# write genotype data, transposed
write.table(t(gstr), file="../genotypes_ind-by-mar.csv", quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)


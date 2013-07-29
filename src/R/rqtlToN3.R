# CrossObject to N3 format
CrossTotripleStore <- function(cross, dataset = "default", file = ""){
  cat("@prefix : <http://www.rqtl.org/ns/#> .\n", file=file)
  cat("@prefix individual: <http://www.rqtl.org/ns/individual#> .\n", file=file, append=TRUE)
  cat("@prefix marker: <http://www.rqtl.org/ns/marker#> .\n",     file=file, append=TRUE)
  cat("@prefix location: <http://www.rqtl.org/ns/location#> .\n",   file=file, append=TRUE)
  cat("@prefix phenotype: <http://www.rqtl.org/ns/phenotype#> .\n",  file=file, append=TRUE)
  map <- pull.map(cross)
  for(chr in names(map)){
    for(mar in names(map[[chr]])){
      mname <- sub("/","_",mar)
      cat(paste0("marker:", mname," location:chr ", chr," .\n"), file=file, append=TRUE)
      cat(paste0("marker:", mname," location:cm ", map[[chr]][mar]," .\n"), file=file, append=TRUE)
      geno <- pull.geno(cross)[,mar]
      for(ind in 1:length(geno)){
        cat(paste0("individual:", ind," marker:", mname," \"", geno[ind],"\" .\n"), file=file, append=TRUE)
      }
    }
  }
  for(phe in names(pull.pheno(multitrait))){
    pheno <- pull.pheno(cross)[,phe]
    for(ind in 1:length(pheno)){
      cat(paste0("individual:", ind," phenotype:", phe," \"", pheno[ind],"\" .\n"), file=file, append=TRUE)
    }
  }
}

# CrossObject to N3->RDF and upload to 4 store tripleStore
testCrossTotripleStore <- function(dataset = "default", file="default", toRDF = FALSE, to4Store = FALSE, host = "localhost", port = 9000){
  library(qtl)
  data(multitrait)
  file <- paste0(file,".n3")
  CrossTotripleStore(multitrait, dataset="multitrait", file=file)
  if(file != ""){
    rdfname <- paste0(dataset, ".rdf")
    cat("For RDF and 4Store upload execute:\n")
    cat(paste0("cd ", getwd(), "\n"))
    if(toRDF) cat(paste0("rdfcat -in N3 -out RDF/XML -n ", file, " > ", rdfname,"\n"))
    if(to4Store) cat(paste0("curl -T '", rdfname, "' http://",host,":",port,"/data/", rdfname,"\n"))
  }
}

testCrossTotripleStore(toRDF = TRUE, to4Store = TRUE);

#./sparql-query http://localhost:8080/sparql/
#PREFIX : <http://www.rqtl.org/ns/#>
#PREFIX individual: <http://www.rqtl.org/ns/individual#>
#PREFIX marker: <http://www.rqtl.org/ns/marker#>
#PREFIX location: <http://www.rqtl.org/ns/location#>
#PREFIX phenotype: <http://www.rqtl.org/ns/phenotype#>
#SELECT * WHERE { ?s ?p ?o };


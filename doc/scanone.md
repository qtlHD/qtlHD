# Scanone

## Introduction

Scanone is a routine that takes as input a vector of markers, genotypes of individuals at
those markers, and phenotypes, and outputs LOD scores at
marker locations. A score beyond a certain LOD threshold is a significant 
quantitative trait locus (QTL).

The markers represent genomic locations on chromosomes, with individuals having
an observed genotype, e.g., the marker/individual matrix:

            M1  M2  M3
    Ind1    AA  AB  BB
    Ind2    BB  AB  BB

The observed genotypes AA, AB and BB get indexed by number, i.e. AA=1, AB=2 and
BB=3. At M2 we have 2 times genotype 3. With more genotypes the number simply 
expands:

            M1  M2  M3
    Ind1     1   2   3
    Ind2     3   2   3

In the next step we expand the markers with intermediate steps. If markers M1-3
are on one chromosome, we introduce pseudo-markers at fixed intervals. So we get

            M1  PM11 M2  PM21 PM22 M3
    Ind1     1   ?    2    ?    ?   3
    Ind2     3   ?    2    ?    ?   3
  
Of course we do not know the true genotypes at the pseudo marker locations, but
we can infer that PM11 for individual 1 is likely to be 1 or 2. To calculate
the probabilities of the inferred genotypes we create a 3D matrix with the
probabilities on the z-axis. The (inferred) genotypeprobabilities are
calculated using an HMM.

Using the genotype probabilities we can use models (where the model is we have
a QTL at a marker location, or we have no QTL, i.e.  the LOD score compares the
probability of obtaining the test data if the two markers are linked to the
probability of obtaining the test data if the two markers are not linked.) to
map phenotypes against genotype. The outcome is a LOD score (logarithm of odds)
at every marker position. Additional parameters are covariates and weights.
Covariates (such as sex, laboratory, or greenhouse, but also other phenotypes
and even genotypes) are additional inputs that are used at every marker
location - i.e.  they are added to the marker/inferred genotype matrix before
the linear regression.  Weights are used to give weight to certain grouped
individuals also before regression (this is rarely used in practise, it could
be used to correct for population size, for example).

## Input files

The input files reflect above inputs. The main exercise of importing data is to
translate a genotype encoding scheme (AA, AB, BB, etc) into observed genotypes. 

When you check the listeria.csv file in ./test/data, it contains a row
of marker names, followed by a row of chromosome numbers for each marker, next 
the marker locations on the chromosomes. The rest of the file consists
of individuals, each row containing a name, the phenotype (a double) followed by
the genotypes at every marker location, in this case A, B, H and - for missing
data.



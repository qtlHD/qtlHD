/**
 * Module for Haley-Knott scanone (single QTL mapping by HK regression).
 */

module qtl.core.scaneone_hk;

/*
 * Performs genome scan using the Haley-Knott regression method
 * (regressing phenotypes on conditional genotype probabilities; the
 * multipoint genotype probabilities have already been calculated in
 * calc.genoprob)
 * 
 * The old R/qtl interface read:
 *
 *   n_ind        Number of individuals
 *   n_pos        Number of marker positions
 *   n_gen        Number of different genotypes
 *   Genoprob     Array of conditional genotype probabilities
 *   Addcov       Matrix of additive covariates: Addcov[cov][ind]
 *   n_addcov     Number of columns of Addcov
 *   Intcov       Number of interactive covariates: Intcov[cov][ind]
 *   n_intcov     Number of columns of Intcov
 *   pheno        Phenotype data, as a vector
 *   nphe         Number of phenotypes
 *   weights      Vector of positive weights, of length n_ind
 *   ind_noqtl    Indicators (0/1) of which individuals should be excluded 
 *                from QTL effects.  
 *   Result       Result matrix of size [n_pos x (nphe)] containing the
 *                LOD scores for each phenotype
 *
 * The new (tentative) interface is:
 *
 * In:
 *
 *   markers          List of markers and positions, including refs to
 *   -> genotype[][]     Matrix of genotypes (one or more)
 *   phenotypes       List of phenotypes
 *   individuals      List of individuals with refs to
 *   -> weights
 *   -> exclude
 *   covariates       List of covariates 
 *
 * Out:
 *  
 *   mappedqtls       The LOD scores for each phenotype
 *
 **********************************************************************/

MappedQTLs scanone_hk(in Markers markers,
                      in Phenotypes phenotypes, 
                      in Individuals individuals,
                      in Covariates covariates)
{
}


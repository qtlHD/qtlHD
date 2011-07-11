/**
 * Module for Haley-Knott scanone (single QTL mapping by HK regression).
 *
 * NOTE: this module is under development
 */

module qtl.core.scanone_hk;

import std.container;
import qtl.core.primitives;
import std.stdio;
import std.algorithm;
import qtl.core.genotype;
import qtl.plugins.input.read_csv;


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
 * Notes: 
 *
 *   This version skips weights, covariates
 *
 **********************************************************************/

double[] scanone_hk(Ms,Ps,Is,Gs)(in Ms markers, in Ps phenotypes, in Is individuals, in Gs genotypes) 
{
  auto sorted_markers = markers.sorted();
  immutable n_mar = sorted_markers.length;
  immutable n_phe = phenotypes.length;
  immutable n_ind = individuals.length;
  immutable n_gen = genotypes.length;
  immutable n_rss = n_phe;
  immutable n_intcov = 0; // later
  immutable n_addcov = 0; // later
  double[] rss;
  rss.reserve(n_rss);
  double[] tmp_pheno;  
  tmp_pheno.reserve(n_ind * n_phe);
  // design matrix for full model (genotypes + cov)
  immutable ncolx = n_gen + (n_gen-1)*n_intcov+n_addcov;
  // immutable rank = ncolx;
  immutable lwork = 3*ncolx + max(n_ind, n_phe);
  auto dwork = new double[][](200,200);
  // dwork.reserve((2*n_ind+1)*ncolx+lwork+(n_ind+ncolx)*n_phe);
  foreach(ref i; dwork)
    i[] = 0;
  auto singular = dwork;
  auto work = singular[ncolx];
  //auto x = work + lwork;
  //auto x_bk = x + n_ind*ncolx;
  //auto yfit = x_bk + n_ind*ncolx;
  //auto coef = yfit + n_ind*n_phe;

  // for(j=0; j<n_ind; j++) 
  //  for(k=0; k<n_phe; k++)
  //    pheno[j+k*n_ind] *= weights[j];

  

  return null;
}

unittest {
  writeln("Unit test " ~ __FILE__);
}


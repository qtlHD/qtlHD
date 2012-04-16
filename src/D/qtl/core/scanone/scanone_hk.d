/**
 * Module for Haley-Knott scanone (single QTL mapping by HK regression).
 *
 * NOTE: this module is under development
 */

module qtl.core.scanone.scanone_hk;

import std.container;
import std.stdio;
import std.algorithm;
import std.math;
import std.string;
import std.conv;
import std.exception;

import qtl.core.primitives;
import qtl.core.scanone.linreg;

// fill up X matrix for scanone, one position
void create_scanone_Xmatrix(in Probability[][] genoprobs, in double[][] addcovar,
                            in double [][] intcovar, in double[] weights,
                            out double[] Xmatrix)
{
  // genoprobs is [individuals][genotypes]
  
  auto n_ind = genoprobs.length;
  auto n_gen = genoprobs[0].length;

  size_t n_addcovar=0, n_intcovar=0;

  if(addcovar.length > 0) n_addcovar = addcovar[0].length;
  if(intcovar.length>0) n_intcovar = intcovar[0].length;

  // output matrix should be n_ind x (n_gen + n_addcovar + (n_gen-1)*n_intcovar)

  size_t i, j, k1, k2, s;

  for(i=0; i<n_ind; i++) {
    // genotype probability columns
    for(j=0; j<n_gen; j++)
      Xmatrix[i+j*n_ind] = genoprobs[i][j]*weights[i];

    // addcovar columns
    for(j=0; j<n_addcovar; j++)
      Xmatrix[i+(j+n_gen)*n_ind] = addcovar[i][j]*weights[i];

    // intcovar columns
    for(k1=0,s=0; k1<n_gen-1; k1++)
      for(k2=0; k2<n_intcovar; k2++,s++) 
        Xmatrix[i+(n_gen+n_addcovar+s)*n_ind] = genoprobs[i][k1]*intcovar[i][k2]*weights[i];
  }
}

// scanone, one position
double[] scanone_hk_onelocus(in Probability[][] genoprobs, in double[][] pheno,
                             in double[][] addcovar, in double[][] intcovar,
                             double[] weights)
{
  auto n_phe = pheno[0].length;
  auto rss = new double[](n_phe);

  return(rss);
}

unittest {
  // See also ./test/scanone/test_scanone.d
  writeln("Unit test " ~ __FILE__);
}


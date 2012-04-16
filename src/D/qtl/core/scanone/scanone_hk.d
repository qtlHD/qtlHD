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
double[] create_scanone_Xmatrix(in Probability[][] genoprobs, in double[][] addcovar,
                                in double [][] intcovar, in double[] weights)
{
  // genoprobs is [individuals][genotypes]
  // addcovar, intcovar are [individuals][covariates]
  
  auto n_ind = genoprobs.length;
  auto n_gen = genoprobs[0].length;

  size_t n_addcovar=0, n_intcovar=0;

  if(addcovar.length > 0) n_addcovar = addcovar[0].length;
  if(intcovar.length>0) n_intcovar = intcovar[0].length;

  // X matrix
  auto ncolx = n_gen + n_addcovar + (n_gen-1)*n_intcovar;
  auto Xmatrix = new double[](n_ind * ncolx);

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

  return Xmatrix;
}

// scanone, one position
//    returns vector of residual sums of squares
double[] scanone_hk_onelocus(in Probability[][] genoprobs, in double[][] pheno,
                             in double[][] addcovar, in double[][] intcovar,
                             double[] weights, double tol=1e-8)
{
  // genoprobs is [individuals][genotypes]
  // addcovar, intcovar are [individuals][covariates]
  // phenotypes is [individuals][phenotypes]
  
  auto n_ind = genoprobs.length;
  auto n_gen = genoprobs[0].length;

  auto n_phe = pheno[0].length;

  size_t n_addcovar=0, n_intcovar=0;

  if(addcovar.length > 0) n_addcovar = addcovar[0].length;
  if(intcovar.length>0) n_intcovar = intcovar[0].length;

  // create simple double[] array for phenotypes
  auto Ymatrix = new double[](n_ind * n_phe);
  foreach(i; 0..n_ind)
    foreach(j; 0..n_phe)
      Ymatrix[i+j*n_ind] = pheno[i][j];
  
  // create simple double[] array with X matrix
  auto ncolx = n_gen + n_addcovar + (n_gen-1)*n_intcovar;
  auto Xmatrix = create_scanone_Xmatrix(genoprobs, addcovar, intcovar, weights);

  auto rss = calc_linreg_rss(Xmatrix, n_ind, ncolx, Ymatrix, n_phe, tol);

  return rss;
}

// scanone, many positions
//    returns matrix of residual sums of squares [position][phenotype]
double[][] scanone_hk(in Probability[][][] genoprobs, in double[][] pheno,
                      in double[][] addcovar, in double[][] intcovar,
                      double[] weights, double tol=1e-8)
{
  // genoprobs is [positions][individuals][genotypes]
  
  auto rss = new double[][](genoprobs.length, pheno[0].length);

  // if weights has length 0, fill with 1's
  if(weights.length==0) {
    weights = new double[](pheno.length);
    foreach(i; 0..pheno.length)
      weights[i] = 1.0;
  }
 
  // scan one position at a time
  foreach(i; 0..genoprobs.length)
    rss[i] = scanone_hk_onelocus(genoprobs[i], pheno, addcovar, intcovar,
                                 weights, tol);

  return rss;
}



unittest {
  // See ./test/scanone/test_scanone.d
  writeln("Unit test " ~ __FILE__);
}


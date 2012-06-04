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
import math.lapack.linreg;

// fill up X matrix for scanone, one position
double[] create_scanone_Xmatrix(in Probability[][] genoprobs, in double[][] addcovar,
                                in double [][] intcovar, in double[] weights)
{
  if(genoprobs.length != weights.length)
    throw new Exception("Mismatch in no. individuals in genoprobs and weights");
  if(addcovar.length > 0 && addcovar[0].length > 0 && genoprobs.length != addcovar.length)
    throw new Exception("Mismatch in no. individuals in genoprobs and addcovar");
  if(intcovar.length > 0 && intcovar[0].length > 0 && genoprobs.length != intcovar.length)
    throw new Exception("Mismatch in no. individuals in genoprobs and intcovar");
  // genoprobs is [individuals][genotypes]
  // addcovar, intcovar are [individuals][covariates]

  auto n_ind = genoprobs.length;
  auto n_gen = genoprobs[0].length;

  size_t n_addcovar=0, n_intcovar=0;

  if(addcovar.length > 0) n_addcovar = addcovar[0].length;
  if(intcovar.length > 0) n_intcovar = intcovar[0].length;

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
double[] scanone_hk_onelocus(in Probability[][] genoprobs, in Phenotype!(double)[][] pheno,
                             in double[][] addcovar, in double[][] intcovar,
                             in double[] weights, in double tol=1e-8)
{
  if(genoprobs.length != pheno.length)
    throw new Exception("Mismatch in no. individuals in genoprobs and pheno");
  if(genoprobs.length != weights.length)
    throw new Exception("Mismatch in no. individuals in genoprobs and weights");
  if(addcovar.length > 0 && addcovar[0].length > 0 && genoprobs.length != addcovar.length)
    throw new Exception("Mismatch in no. individuals in genoprobs and addcovar");
  if(intcovar.length > 0 && intcovar[0].length > 0 && genoprobs.length != intcovar.length)
    throw new Exception("Mismatch in no. individuals in genoprobs and intcovar");

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
      Ymatrix[i+j*n_ind] = pheno[i][j].value;

  // create simple double[] array with X matrix
  auto ncolx = n_gen + n_addcovar + (n_gen-1)*n_intcovar;
  auto Xmatrix = create_scanone_Xmatrix(genoprobs, addcovar, intcovar, weights);

  auto rss = calc_linreg_rss(Xmatrix, n_ind, ncolx, Ymatrix, n_phe, tol);

  return rss;
}

// scanone, many positions
//    returns matrix of residual sums of squares [position][phenotype]
double[][] scanone_hk(in Probability[][][] genoprobs, in Phenotype!(double)[][] pheno,
                      in double[][] addcovar, in double[][] intcovar,
                      in double[] weights, in double tol=1e-8)
{
  if(genoprobs.length == 0)
    throw new Exception("genoprobs is empty");
  if(genoprobs[0].length != pheno.length)
    throw new Exception("Mismatch in no. individuals in genoprobs and pheno");
  if(weights.length > 0 && genoprobs[0].length != weights.length)
    throw new Exception("Mismatch in no. individuals in genoprobs and weights");
  if(addcovar.length > 0 && addcovar[0].length > 0 && genoprobs[0].length != addcovar.length)
    throw new Exception("Mismatch in no. individuals in genoprobs and addcovar");
  if(intcovar.length > 0 && intcovar[0].length > 0 && genoprobs[0].length != intcovar.length)
    throw new Exception("Mismatch in no. individuals in genoprobs and intcovar");
  // genoprobs is [positions][individuals][genotypes]

  auto rss = new double[][](genoprobs.length, pheno[0].length);

  // if weights has length 0, fill with 1's
  auto local_weights = new double[](pheno.length);
  if(weights.length==0) {
    foreach(i; 0..pheno.length)
      local_weights[i] = 1.0;
  }
  else {
    local_weights = weights.dup;
  }

  // scan one position at a time
  foreach(i; 0..genoprobs.length)
    rss[i] = scanone_hk_onelocus(genoprobs[i], pheno, addcovar, intcovar,
                                 local_weights, tol);

  return rss;
}


// one phenotype version
double[] scanone_hk(in Probability[][][] genoprobs, in Phenotype!(double)[] pheno,
                    in double[][] addcovar, in double[][] intcovar,
                    in double[] weights, in double tol=1e-8)
{
  auto pheno_2d = new Phenotype!double[][](pheno.length, 1);
  foreach(i; 0..pheno.length)
    pheno_2d[i][0].value = pheno[i].value;

  auto rss_2d = scanone_hk(genoprobs, pheno_2d, addcovar, intcovar, weights, tol);
  auto rss_1d = new double[](rss_2d.length);
  foreach(i; 0..rss_2d.length)
    rss_1d[i] = rss_2d[i][0];

  return rss_1d;
}


// scanone for null model
double[] scanone_hk_null(in Phenotype!(double)[][] pheno, in double[][] addcovar,
                         in double[] weights, in double tol=1e-8)
{
  if(weights.length > 0 && pheno.length != weights.length)
    throw new Exception("Mismatch in no. individuals in pheno and weights");
  if(addcovar.length > 0 && addcovar.length > 0 && pheno.length != addcovar.length)
    throw new Exception("Mismatch in no. individuals in pheno and addcovar");

  auto vector_of_ones = new Probability[][](pheno.length,1);
  foreach(i; 0..pheno.length)
    vector_of_ones[i][0] = 1.0;

  // if weights has length 0, fill with 1's
  auto local_weights = new double[](pheno.length);
  if(weights.length==0) {
    foreach(i; 0..pheno.length)
      local_weights[i] = 1.0;
  }
  else {
    local_weights = weights.dup;
  }

  auto intcovar = new double[][](0,0);

  return scanone_hk_onelocus(vector_of_ones, pheno,
                             addcovar, intcovar, local_weights, tol);
}

// one phenotype version
double scanone_hk_null(in Phenotype!(double)[] pheno, in double[][] addcovar, 
                       in double[] weights, in double tol=1e-8)
{
  auto pheno_2d = new Phenotype!double[][](pheno.length, 1);
  foreach(i; 0..pheno.length)
    pheno_2d[i][0].value = pheno[i].value;

  auto rss = scanone_hk_null(pheno_2d, addcovar, weights, tol);

  return rss[0];
}


// calculate lod scores
double[][] rss_to_lod(in double[][] rss_alt, in double[] rss_null, in size_t n_ind)
{
  if(rss_alt[0].length != rss_null.length)
    throw new Exception("Mismatch in lengths between rss_alt and rss_null");

  auto lod = new double[][](rss_alt.length, rss_alt[0].length);

  foreach(i; 0..rss_alt.length)
    foreach(j; 0..rss_alt[i].length)
      lod[i][j] = rss_to_lod(rss_alt[i][j], rss_null[j], n_ind);

  return lod;
}

// one phenotype version
double[] rss_to_lod(in double[] rss_alt, in double rss_null, in size_t n_ind)
{
  auto lod = new double[](rss_alt.length);

  foreach(i; 0..rss_alt.length)
    lod[i] = rss_to_lod(rss_alt[i], rss_null, n_ind);

  return lod;
}

// single value version
double rss_to_lod(in double rss_alt, in double rss_null, in size_t n_ind)
{
  return cast(double)n_ind/2.0 * log10(rss_null / rss_alt);
}


unittest {
  // See ./test/scanone/test_scanone.d
  writeln("Unit test " ~ __FILE__);
}


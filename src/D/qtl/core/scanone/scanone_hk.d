/**
 * Module for Haley-Knott scanone (single QTL mapping by HK regression).
 *
 * NOTE: this module is under development
 */

module qtl.core.scanone.scanone_hk;

import std.container;
import qtl.core.primitives;
import std.stdio;
import std.algorithm;
import std.math;
import std.string;
import std.conv;
import std.exception;
import qtl.core.deprecate.genotype_enum;
import qtl.plugins.csv.read_csv;

// import core.stdc.stdlib;  // for malloc
// import core.stdc.string;  // for memcpy

static immutable bool VERBOSE = true;

void message(string s) {
  static if (VERBOSE)
    writeln("DEBUG: " ~ s);
  else
    debug(1) { write("DEBUG: "); writeln(s); }
}

version (Windows) {
  import qtl.core.libs.libload;
  import std.loader;
}else{
  pragma(lib, "blas");
  pragma(lib, "lapack");
}

immutable TOL = 1e-12;  // tolerance for linear regression

// Prototypes for the raw Fortran interface to BLAS
extern(C) {

alias float f_float;
alias double f_double;
alias cfloat f_cfloat;
alias cdouble f_cdouble;
alias int f_int;

version (Windows) {
  void function(char *trans, f_int *m, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *work, f_int *lwork, f_int *info) dgels_;
  void function (f_int *m, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *s, f_double *rcond, f_int *rank, f_double *work, f_int *lwork, f_int *info)dgelss_;
  void function (char *transa, char *transb, f_int *m, f_int *n, f_int *k, f_double *alpha, f_double *A, f_int *lda, f_double *B, f_int *ldb, f_double *beta, f_double *C, f_int *ldc) dgemm_;
  void function (char *uplo, f_int *n, f_double *a, f_int *lda, f_int *info) dpotrf_;
  void function (char *uplo, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_int *info) dpotrs_;

  static this(){
    HXModule blaslib = load_library("Rblas");
    load_function(dgemm_)(blaslib,"dgemm_");
    if (VERBOSE) message("mapped Rblas.dll");

    HXModule lapacklib = load_library("Rlapack");
    load_function(dgels_)(lapacklib,"dgels_");
    load_function(dgelss_)(lapacklib,"dgelss_");
    load_function(dpotrf_)(lapacklib,"dpotrf_");
    load_function(dpotrs_)(lapacklib,"dpotrs_");
    if (VERBOSE) message("mapped Rlapack.dll");
  }
} 
else {

  void dgels_(char *trans, f_int *m, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *work, f_int *lwork, f_int *info);

  void dgelss_(f_int *m, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *s, f_double *rcond, f_int *rank, f_double *work, f_int *lwork, f_int *info);

  // void dgemm_(char *transa, char *transb, f_int *m, f_int *n, f_int *k, f_double *alpha, f_double *A, f_int *lda, f_double *B, f_int *ldb, f_double *beta, f_double *C, f_int *ldc, f_int transa_len, f_int transb_len);
  void dgemm_(char *transa, char *transb, f_int *m, f_int *n, f_int *k, f_double *alpha, f_double *A, f_int *lda, f_double *B, f_int *ldb, f_double *beta, f_double *C, f_int *ldc);

  // void dpotrf_(char *uplo, f_int *n, f_double *a, f_int *lda, f_int *info, f_int uplo_len);
  void dpotrf_(char *uplo, f_int *n, f_double *a, f_int *lda, f_int *info);

  // void dpotrs_(char *uplo, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_int *info, f_int uplo_len);
  void dpotrs_(char *uplo, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_int *info);

  }
}

/**********************************************************************
 *
 * lapackutil.c
 *
 * copyright (c) 2006, Hao Wu
 *
 * last modified Feb, 2006 
 * first written Jan, 2006 
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU General Public License,
 *     version 3, as published by the Free Software Foundation.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but without any warranty; without even the implied warranty of
 *     merchantability or fitness for a particular purpose.  See the GNU
 *     General Public License, version 3, for more details.
 * 
 *     A copy of the GNU General Public License, version 3, is available
 *     at http://www.r-project.org/Licenses/GPL-3
 *
 * C functions for the R/qtl package
 *
 * These are some wrapper functions for several LAPACK routines.
 *
 * Contains: mydgelss, mydgemm, mydpotrf, mydpotrs
 *
 **********************************************************************/

/* DGELSS function */
private void mydgelss (int n_ind, int ncolx0, int n_phe, double x0[], double x0_bk[],
               double *pheno, double tmppheno[], double singular[], double tol, 
               int rank, double *work, int *lwork)
{
  bool is_singular=false;

  /* 
     SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,
     $                  INFO )
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
   
      M      = rows in A
      N      = cols in A
      NRHS   = number of right hand sides
      A      = (input/output) DOUBLE PRECISION array, dimension (LDA,N)
      LDA    = The leading dimension of the array A.  LDA >= max(1,M).
      B      = (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
               the matrix B of right hand side vectors, stored
               B is M-by-NRHS 
      LDB    = The leading dimension of the array B. LDB >= MAX(1,M,N).
      WORK   = (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
      LWORK  = INTEGER The dimension of the array WORK.
      INFO    (output) INTEGER  0:  successful exit


         ncolx0    n_pheno
         ngen+cov
         1 .. N  | NRHS   |
                 
      1  
 nind ..
      M


F77_NAME(dgels)(const char* trans, const int* m, const int* n,
                const int* nrhs, double* a, const int* lda,
                double* b, const int* ldb,
                double* work, const int* lwork, int* info);

  */
  writeln(n_ind); writeln(ncolx0); writeln(n_phe); writeln(x0[0..10]); writeln(n_ind); writeln(*work); writeln(*lwork); 

  auto m = n_ind;   // ind
  auto n = ncolx0;  // gen + cov
  auto nrhs = n_phe;
  auto lda = n_ind;
  double *a = x0.ptr; // (gen + cov + phe) x ind
  double *b = tmppheno.ptr;
  auto ldb = n_ind;

  auto info = 0;
  assert(lda >= max(1,m));
  writeln("x0:",x0[0]);
  assert(!isnan(x0[0]), to!string(x0[0]));
  message("dgels_");
  dgels_(cast(char *)toStringz("N"), &m, &n, &nrhs, a, &lda, b, &ldb,
      work, lwork, &info);

  /* if there's problem like singular, use dgelss */
  /* note that x0 will contain the result for QR decomposition. 
  If any diagonal element of R is zero, then input x0 is rank deficient */
  for(auto i=0; i<ncolx0; i++)  {
    if(abs(x0[n_ind*i+i]) < TOL) {
      is_singular = true;
      break;
    }
  }

  
  if(is_singular) { /* switch to dgelss if input x0 is not of full rank */
    /* note that tmppheno and x0 have been destroyed already,
    we need to make another copy of them */
    /*mexPrintf("Warning - Design matrix is signular \n"); */
    /* note that at this stage both x0 and tmppheno might be destroyed,
    we need to make a copy of them */
    message("singular");

    // memcpy(x0, x0_bk, *n_ind*(*ncolx0)*double.sizeof);
    for (auto idx=0; idx<n_ind*(ncolx0); idx++)
      x0[idx] = x0_bk[idx];
    // memcpy(tmppheno, pheno, *n_ind*(*n_phe)*double.sizeof);
    for (auto idx=0; idx<n_ind*(n_phe); idx++)
      tmppheno[idx] = pheno[idx];
    info = 0;
    dgelss_(&n_ind, &ncolx0, &n_phe, a, &n_ind, tmppheno.ptr, &n_ind, 
      singular.ptr, &tol, &rank, work, lwork, &info);
  }
}


/* DGEMM */
private void mydgemm(int *n_phe, int *n_ind, double *alpha, double *tmppheno, double *beta, double *rss_det) 
{
  dgemm_(cast(char *)toStringz("T"),cast(char *)toStringz("N"), n_phe, n_phe, n_ind, alpha, tmppheno, n_ind, tmppheno, n_ind, beta, rss_det, n_phe);
}

/* DPOTRF */
private void mydpotrf(int *n_phe, double *rss_det, int *info) 
{
  dpotrf_(cast(char *)toStringz("U"), n_phe, rss_det, n_phe, info);
}

/*DPOTRS */
private void mydpotrs(char *uplo, int *n, int *nrhs, double *A, 
                      int *lda, double *B, int *ldb, int *info)
{
  dpotrs_(uplo, n, nrhs, A, lda, B, ldb, info);
}

/* end of lapackutil.c */

private double [] matmult(in double a[], int nrowa, int ncola, in double b[], int ncolb)

{
  int i, j, k;
  double result[] = new double[](nrowa*ncola);

  for(i=0; i<nrowa; i++) {
    for(j=0; j<ncolb; j++) {
      /* clear the content of result */
      // result[j*nrowa+i] = 0.0;
      /*result[i*ncolb+j] = 0.0;*/
      for(k=0; k<ncola; k++)
        result[j*nrowa+i] += a[k*nrowa+i]*b[j*ncola+k];
    }
  }
  return result;
}

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
 *   Genoprob     Matrix of conditional genotype probabilities
 *   addcov       Matrix of additive covariates: addcov[cov][ind] (nyi)
 *   n_addcov     Number of columns of addcov (nyi)
 *   intcov       Number of interactive covariates: intcov[cov][ind] (nyi)
 *   n_intcov     Number of columns of intcov (nyi)
 *   pheno        Phenotype data, as a vector -- matix?
 *   n_phe        Number of phenotypes
 *   weights      Vector of positive weights, of length n_ind
 *   ind_noqtl    Indicators (0/1) of which individuals should be excluded 
 *                from QTL effects.  
 *   Result       Result matrix of size [n_pos x (n_phe)] containing the
 *                LOD scores for each phenotype
 *
 * The new (tentative) interface is:
 *
 * Notes: 
 *
 *   This version skips weights, covariates
 *
 * Check: http://www.netlib.org/lapack/lug/node71.html
 *
 **********************************************************************/

double[] scanone_hk(Is,Ps,Ms)(in Is individuals, in Ps phenotypes, in Ms markers, in GenoProbs genoprob) 
{
  // inputs
  message("Starting scanone_HK");
  // immutable double **addcov, intcov;

  // calculated
  immutable n_pos = markers.length;
  int n_phe = phenotypes.length;
  int n_ind = individuals.length;
  immutable n_gen = genoprob[0].length; // get columns of matrix

  // unused now
  int[] ind_noqtl;
  ind_noqtl.length = n_ind; 
  assert(ind_noqtl[0] == 0);
  double[] weights;
  weights.length = n_ind;
  weights[] = 1.0;
  immutable n_intcov = 0;
  immutable n_addcov = 0;
  // initialize
  auto pheno_size = n_ind * n_phe; // matrix of ind x phenotypes
  auto pheno_memsize = pheno_size * double.sizeof;
  auto pheno = new double[pheno_size];
  // genoprob = new double[][][](n_gen,n_pos,n_ind);
  assert(!isnan(genoprob[0][0][0].value));

  // local
  int  i, j, k, k2, s, nrss, lwork, ind_idx;
  double alpha=1.0;
  auto beta=0.0;
  double tol=TOL;
  double dtmp;

  /* number of rss's, currently multivar is not used so it's always 0 */
  nrss = n_phe;

  /* allocate memory */
  auto rss = new double[nrss];  

  auto tmppheno = new double[pheno_size];
  /* number of columns in design matrix X for full model */
  int ncolx = n_gen + (n_gen-1)*n_intcov+n_addcov; 
  // 1..n_gen ; 1..n_gen * n_intcov ; n_addcov
  /*ncol0 = n_addcov+1;*/
  auto rank = ncolx;

  /* allocate space and set things up*/
  /*  x = (double *)R_alloc(n_ind*ncol, double.sizeof);
  coef = (double *)R_alloc(ncol, double.sizeof);
  resid = (double *)R_alloc(n_ind, double.sizeof);
  qty = (double *)R_alloc(n_ind, double.sizeof);
  jpvt = (int *)R_alloc(ncol, sizeof(int));
  qraux = (double *)R_alloc(ncol, double.sizeof);
  work = (double *)R_alloc(2 * ncol, double.sizeof); */
  lwork = 3*ncolx + max(n_ind, n_phe);
  // 120 n_ind, 133 ncolx (n_geno+cov), 519 lwork, 58200 dwork_size
  auto dwork_size = (2*n_ind+1)*ncolx+lwork+(n_ind+ncolx)*n_phe;
  writeln("!!",n_ind,',',n_phe,',',ncolx,',',lwork,',',dwork_size);
  auto dwork = new double[dwork_size];
  dwork[] = 0.0;
  assert(dwork[0]==0.0);
  assert(dwork[dwork_size-10]==0.0);

  /* split the memory block 
    singular (lwork) | x (ind*ncolx) | x_bk (ind*ncolx) | yfit (ind*ncolx) | coef (n_ind*phe)
  */
  auto singular = dwork;
  auto work = singular[ncolx..$];
  assert(work[0]==0.0);
  auto x = work[lwork..$];
  auto x_bk = x[n_ind*ncolx..$];
  assert(x_bk[0]==0.0);
  auto yfit = x_bk[n_ind*ncolx..$];
  auto coef = yfit[n_ind*n_phe..$];
  x.length = n_ind * ncolx;
  assert(x[0]==0.0);
  x_bk.length = n_ind * ncolx; 
  yfit.length = n_ind * ncolx; 
  coef.length = n_ind * n_phe;
  assert(work.length == 62799, to!string(work.length));

  /* NULL model is now done in R ********************
     (only do it once!)
  for(j=0; j<n_ind; j++) {
    x[j] = 1.0;
    for(k=0; k<n_addcov; k++) 
      x[j+(k+1)*n_ind] = addcov[k][j];
  }
  F77_CALL(dqrls)(x, &n_ind, &ncol0, pheno, &ny, &tol, coef, resid,
                  qty, &k, jpvt, qraux, work);
  rss0 = 0.0;
  for(j=0; j<n_ind; j++)  rss0 += (resid[j]*resid[j]);
  Null model is now done in R ********************/
  // disabling weights, for now

  /*
  for(j=0; j<n_ind; j++) 
    for(k=0; k<n_phe; k++)
      pheno[j+k*n_ind] *= weights[j];  note: weights are really square-root of weights 
  */
  for(i=0; i<n_pos; i++) { /* loop over positions (markers) */
    // for(k=0; k<n_ind * ncolx; k++) x[k] = 0.0;
    x[] = 0.0;
    assert(x[0]==0.0);

    /* fill up X matrix */
    for(j=0; j<n_ind; j++) {
      if(!ind_noqtl[j]) {
        for(k=0; k<n_gen; k++)
          x[j+k*n_ind] = genoprob[k][i][j].value*weights[j];
      }
      // for(k=0; k<n_addcov; k++)
      //  x[j+(k+n_gen)*n_ind] = addcov[k][j]*weights[j];
      // if(!ind_noqtl[j]) {
      //   for(k=0,s=0; k<n_gen-1; k++)
      //     for(k2=0; k2<n_intcov; k2++,s++) 
      //       x[j+(n_gen+n_addcov+s)*n_ind] = genoprob[k][i][j]*intcov[k2][j]*weights[j];
      //}
    }
    /* linear regression of phenotype on QTL genotype probabilities */
    /*    F77_CALL(dqrls)(x, &n_ind, &ncol, pheno, &ny, &tol, coef, resid,
                    qty, &k, jpvt, qraux, work);
    */
    /* make a copy of x matrix, we may need it */
    // memcpy(x_bk.ptr, x.ptr, n_ind*ncolx*double.sizeof);
    
    /* make a copy of phenotypes. I'm doing this because 
       dgelss will destroy the input rhs array */
    // note using memcpy can be dangerous... FIXME
    for(auto idx=0; idx < pheno.length; idx++)  
      tmppheno[idx] = pheno[idx];

    /* linear regression of phenotype on QTL genotype probabilities */
    mydgelss(n_ind, ncolx, n_phe, x, x_bk, pheno.ptr, tmppheno, singular,
      tol, rank, work.ptr, &lwork);
    /* calculate residual sum of squares */
    if(n_phe == 1) {
      /* only one phenotype, this is easier */
      /* if the design matrix is full rank */
      if(rank == ncolx) {
        for (k=rank, rss[0]=0.0; k<n_ind; k++)
          rss[0] += tmppheno[k]*tmppheno[k];
      }
      else {
        /* the desigm matrix is not full rank, this is trouble */
        /* calculate the fitted value */
        yfit = matmult(x_bk, n_ind, ncolx, tmppheno, 1);
        /* calculate rss */
        for (k=0, rss[0]=0.0; k<n_ind; k++)
          rss[0] += (pheno[k]-yfit[k]) * (pheno[k]-yfit[k]);
      }
    }
    else { /* multiple phenotypes */
        if(rank == ncolx) {
          for(k=0; k<nrss; k++) {
            ind_idx = k*n_ind;
            for(j=rank, rss[k]=0.0; j<n_ind; j++) {
              dtmp = tmppheno[ind_idx+j];
              rss[k] += dtmp * dtmp;
            }
          }
        }
        else { /* not full rank, this is troubler */
          /* note that the result tmppheno has dimension n_ind x n_phe,
          the first ncolx rows contains the estimates. */
          for (k=0; k<n_phe; k++) 
            for(auto idx=0; idx < ncolx; idx++) 
              coef[k*ncolx+idx] = tmppheno[k*n_ind+idx];
          /* calculate yfit */
          yfit = matmult(x_bk, n_ind, ncolx, coef, n_phe);
          /* calculate residual, put the result in tmppheno */
          for (k=0; k<n_ind*n_phe; k++)
            tmppheno[k] = pheno[k] - yfit[k];
          /* calculate rss */
          for(k=0; k<nrss; k++) {
            ind_idx = k*n_ind;
            for(j=0, rss[k]=0.0; j<n_ind; j++) {
              dtmp = tmppheno[ind_idx+j];
              rss[k] += dtmp * dtmp;
            }
          }
          
        }
        
    }
    /* make the result */
    /* log10 likelihood */
    //@for(k=0; k<nrss; k++) 
    //@   Result[k][i] = (double)n_ind/2.0*log10(rss[k]);

  } /* end loop over positions */
  return null;
}

unittest {
  // See also ./test/scanone/test_scanone.d
  writeln("Unit test " ~ __FILE__);
}


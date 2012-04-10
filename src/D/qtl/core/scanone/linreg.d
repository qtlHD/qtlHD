/*
 * linreg: linear regression utility functions
 **/

// was using dqrls, as that is what lm() uses, but that's Linpack
// The Lapack function is dgelss; see R-2.15.0/src/modules/lapack/dlapack1.f
// ...looks like I should use dgelsd rather than dgelss
// dgelsd uses SVD; dgelsy uses QR
// dgels uses QR, but X matrix must be full rank
//
// Lapack benchmarks at http://www.netlib.org/lapack/lug/node71.html
// dgelsy indistinguishable from dgels
// dgelsd 3-5x slower; dgelss 7-34x slower


module qtl.core.scanone.linreg;

import std.algorithm;
import std.math;
import std.range;
import std.c.stdio;
import std.stdio;
import std.string;
import std.conv;

import std.c.stdlib;

extern(C) {
  alias float f_float;
  alias double f_double;
  alias int f_int;
}

version(Windows){
  private import std.loader;
  private import qtl.core.util.windows;
  private import qtl.plugins.renv.libload;

  extern (C) void function(f_int *m, f_int *n, f_int *nrhs, f_double *A, f_int *lda,
                           f_double *B, f_int *ldb, f_int *jpvt, f_double *rcond,
                           f_int *rank, f_double *work, f_int *lwork, f_int *info) dgelsy_;

  static this(){
    HXModule lib = load_library("Rblas");
    load_function(dqrls_)(lib,"dgelsy_");
    writeln("Loaded BLAS functionality");
  }

}else{
  pragma(lib, "blas");
  pragma(lib, "lapack");

  // The C interface from the SciD library:
  extern (C) void dgelsy_(f_int *m, f_int *n, f_int *nrhs, f_double *A, f_int *lda,
                          f_double *B, f_int *ldb, f_int *jpvt, f_double *rcond,
                          f_int *rank, f_double *work, f_int *lwork, f_int *info);
}

// The D interface is a D-ified call which calls the C interface dgelsy_
// see 
void gelsy(f_int m,         // number of rows in A
           f_int n,         // number of columns in A
           f_int nrhs,      // number of right-hand sides (no. columns in B)
           f_double *A,     // [m x n] covariate matrix
           f_int lda,       // leading dimension of A [== m]
           f_double *B,     // [m x nrhs] outcome matrix
           f_int ldb,       // leading dimension of B [== m]
           f_int *jpvt,     // n-vector to keep track of reordering of columns of A
           f_double rcond,  // used to determine the effective rank of A (condition number < 1/rcond)
           f_int *rank,      // on output, the rank of A
           f_double *work,  // [lwork] vector of workspace
           f_int lwork,    // dimension of work (should be >= max(mn+3n+1, 2*mn+nrhs
           f_int *info)     // on output, =0 indicates success; =-i indicates ith argument had illegal value
{
  // see R-2.15.0/src/appl/dgelsy.f

  dgelsy_(&m, &n, &nrhs, A, &lda, B, &ldb, jpvt, &rcond, rank, work, &lwork, info);

  if(*info!=0) throw new Exception("dgelsy_: illegal value in argument " ~ to!string(*info));
}


// fit linear regression model and return residual sum of squares
double[] calc_linreg_rss(double x[], int nrow, int ncolx, double y[], int ncoly, double tol=1e-7)
{
  int nrhs=ncoly, lda=nrow, ldb=nrow, info, rank;
  auto lwork = max(nrow*ncolx + 3*ncolx + 1, 2*nrow*ncolx+nrhs);
  auto work = new double[lwork];
  auto jpvt = new int[ncolx];

  foreach(i; 0..ncolx) jpvt[i] = i+1;  // keeps track of pivoted columns

  double rcond = tol;

  rank = info = -999;
  gelsy(nrow, ncolx, nrhs, x.ptr, lda, y.ptr, ldb, jpvt.ptr, rcond, &rank, work.ptr, lwork, &info);

  auto rss = new double[ncoly];
  foreach(i; 0..ncoly) rss[i]=0.0; // fill with 0's

  if(rank == ncolx) { // X is of full rank
    auto row_index = 0;

    // in each column of y:
    //  first rank values = estimated coefficients
    //  sum of squares of the rest gives RSS (residual sum of squares)
    foreach(i; 0..ncoly) {  
      foreach(j; rank..nrow)
        rss[i] += y[row_index+j]^^2;
      row_index += nrow;
    }
  }
  else {
  }

  return rss;
}


unittest {
  writeln("Unit test " ~ __FILE__);
  writeln("  --X matrix with full rank");

  double[] x = [ 8, 5, 4, 3, 3,
                 9, 2, 2, 2, 9,
                 5, 1, 8, 1, 4];

  double[] y = [42, 12, 32, 10, 33];


  double rss_R = 0.42328733709665872231;

  int nrow = cast(int)y.length;
  int ncolx = cast(int)x.length / nrow;
  int ncoly = 1;

  auto rss = calc_linreg_rss(x, nrow, ncolx, y, ncoly);

  assert(abs(rss[0] - rss_R) < 1e-12);
}

/* this part isn't ready
unittest {
  writeln("Unit test " ~ __FILE__);
  writeln("  --X matrix with full rank");

  double[] x = [ 8, 5, 4, 3, 3,
                 9, 2, 2, 2, 9,
                 11,9, 0, 5, 2,
                 5, 1, 8, 1, 4
                 6, 5, 3, 2, 1];

  double[] y = [42, 12, 32, 10, 33];


  double rss_R = 0.42328733709665872231;

  int nrow = cast(int)y.length;
  int ncolx = cast(int)x.length / nrow;
  int ncoly = 1;

  auto rss = calc_linreg_rss(x, nrow, ncolx, y, ncoly);

  assert(abs(rss[0] - rss_R) < 1e-12);
}
*/
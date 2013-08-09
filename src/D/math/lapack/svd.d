/*
 * SVG: singular value decomposition
 *
 * dgesvd calculates singular value decomposition
 */

module math.lapack.svd;

import std.algorithm;
import std.math;
import std.range;
import std.c.stdio;
import std.stdio;
import std.string;
import std.conv;
import std.c.stdlib;

import math.matrix;

extern(C) {
  alias float f_float;
  alias double f_double;
  alias int f_int;
  alias char f_char;
}

version(Windows){
  private import std.loader;
  private import arch.windows;
  private import qtl.plugins.renv.libload;

  extern (C) void function(f_char *jobu, f_char *jobvt, f_int *m, f_int *n, f_double *A,
                           f_int *lda, f_double *S, f_double *U, f_int *ldu, f_double *VT,
                           f_int *ldvt, f_double *work, f_int *lwork, f_int *info) dgesvd_;

  static this(){
    HXModule lib_lapack = load_library("Rlapack");
    load_function(dgesvd_)(lib_lapack,"dgesvd_");
    writeln("Loaded Rlapack functionality");
  }

}else{
  pragma(lib, "lapack");

  // Lapack routine for SVD
  extern (C) void dgesvd_(f_char *jobu, f_char *jobvt, f_int *m, f_int *n, f_double *A,
                          f_int *lda, f_double *S, f_double *U, f_int *ldu, f_double *VT,
                          f_int *ldvt, f_double *work, f_int *lwork, f_int *info);
}

// The D interface is a D-ified call which calls the C interface dgesvd_
void gesvd(f_char jobu,    // whether to consider left singular vectors (='N' for no)
           f_char jobvt,   // whether to consider right singular vectors (='N' for no)
           f_int m,        // number of rows
           f_int n,        // number of columns
           f_double *A,     // the key matrix [rows x columns]
           f_int lda,      // leading dimension of A (= nrow)
           f_double *S,     // to contain singular values, length = min(m, n)
           f_double *U,     // to contain left singular vectors
           f_int ldu,      // leading dimension of U
           f_double *VT,    // to contain right singular vectors
           f_int ldvt,     // leading dimension of VT
           f_double *work,  // workspace of dimension lwork
           f_int lwork,    // > max { 3*min(m,n) + max(m,n) , 5*min(m,n) }
           f_int *info)     // on output, =0 indicates success; =-i indicates ith argument had illegal value; =+i if converge problem
{
  // see R-3.0.1/src/module/lapack/dlapack.f

  dgesvd_(&jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt, work, &lwork, info);

  if(*info<0) throw new Exception("dgesvd_: illegal value in argument " ~ to!string(*info));
}


// calculate singular values for a matrix
double[] calc_singular_values(double x[], size_t nrow, size_t ncol)
{
  int lda=cast(int)nrow, info, rank;

  int lwork = cast(int)max(3*min(nrow,ncol) + max(nrow,ncol), 5*min(nrow,ncol));
  auto work = new double[lwork];

  auto result = new double[min(nrow, ncol)];

  int ldu=1, ldvt=1;
  auto U = new double[1];
  auto VT = new double[1];

  gesvd('N', 'N', cast(int)nrow, cast(int)ncol, x.ptr, lda, result.ptr, U.ptr, ldu, VT.ptr, ldvt, work.ptr, lwork, &info);

  return result;
}

// same, for doubly-indexed array X[row][col]
double[] calc_singular_values(double [][]matrix)
{
  auto vector_tuple = matrix_as_vector(matrix);

  return calc_singular_values(vector_tuple[0], vector_tuple[1], vector_tuple[2]);
}

// matrix rank
int matrix_rank(double x[], size_t nrow, size_t ncol, double tol=1e-16)
{
  auto sv = calc_singular_values(x, nrow, ncol);

  tol *= sv[0]*max(nrow, ncol);
  int rank=0;
  foreach(s; sv)
    if(s > tol) rank++;

  return(rank);
}

int matrix_rank(double[][] matrix, double tol=1e-16)
{
  auto vector_tuple = matrix_as_vector(matrix);

  return matrix_rank(vector_tuple[0], vector_tuple[1], vector_tuple[2], tol);
}


unittest {
  writeln("Unit test " ~ __FILE__);

  auto x = [ [4.75, 3.12, 3.37, 3.81, 3.95],
             [4.89, 3.81, 2.79, 2.95, 3.63],
             [2.91, 2.40, 2.89, 3.39, 3.74],
             [2.63, 3.70, 3.05, 3.62, 2.72] ];

  // singular values calculated in R: print(svd(x)$d, digits=15)
  auto sv = [15.346690975369022, 1.721557287148784, 1.288596129368399, 0.437077371833305];

  auto xv = matrix_as_vector(x);

  auto sv1 = calc_singular_values(x);
  auto sv2 = calc_singular_values(xv[0], xv[1], xv[2]);

  foreach(i, val; sv)
    writeln(val, " ", sv1[i], " ", sv2[i]);

  assert(sv1.length == sv.length);
  assert(sv2.length == sv.length);
  foreach(i, svv; sv) {
    assert(abs(sv1[i] - svv) < 1e-12);
    assert(abs(sv2[i] - svv) < 1e-12);
  }

  assert(matrix_rank(x) == 4);
}

unittest {
  writeln("A less-than-full rank matrix");

  auto x = [ [1.0, 1.0, 1.0, 1.0],
             [1.0, 1.0, 1.0, 1.0],
             [1.0, 1.0, 1.0, 1.0],
             [1.0, 1.0, 1.0, 1.0],
             [1.0, 2.0, 1.0, 0.0],
             [1.0, 2.0, 1.0, 0.0],
             [1.0, 2.0, 1.0, 0.0],
             [1.0, 2.0, 2.0, 0.0],
             [1.0, 2.0, 2.0, 0.0],
             [1.0, 2.0, 2.0, 0.0] ];

  assert(matrix_rank(x) == 3);

  auto xv = matrix_as_vector(x);
  assert(matrix_rank(xv[0], xv[1], xv[2]) == 3);
}

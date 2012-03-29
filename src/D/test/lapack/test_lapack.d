/**

The following example uses this function to calculate the matrix-vector product

  / 3 1 3 \   / -1 \ 
  | 1 5 9 | * | -1 |. 
  \ 2 6 5 /   \  1 / 

D2 version based on the C version at http://www.seehuhn.de/pages/linear

Compile and run with

  dmd -oftest_lapack test_lapack.d

Or use gdb

  dmd -gc -debug test_lapack.d ; gdb ./test_lapack 

Should produce

  3.0  1.0  3.0
  1.0  5.0  9.0
  2.0  6.0  5.0
 -1.0
  3.0
 -3.0


*/

module test.lapack.test_lapack;

import std.algorithm;
import std.range;
import std.c.stdio;
import std.stdio;
import std.string;

pragma(lib, "blas");
pragma(lib, "lapack");

double m[9] = [
  3, 1, 3,
  1, 5, 9,
  2, 6, 5
];

double m2[9] = [
  3, 1, 2,
  1, 5, 6,
  3, 9, 5
];


double x[3] = [
  -1, -1, 1
];

double y[3] = [
  0, 0, 0
];

extern(C) {
  alias float f_float;
  alias double f_double;
  alias int f_int;
}

// The C interface from the SciD library:
extern (C) void dgemv_(char *trans, f_int *m, f_int *n, f_double *alpha, f_double *A, f_int *lda, f_double *x, f_int *incx, f_double *beta, f_double *y, f_int *incy, f_int trans_len);

// The D interface is a D-ified call which calls the C interface dgemv_
void gemv(char trans, f_int m, f_int n, f_double alpha, f_double *A, f_int lda,
f_double *x, f_int incx, f_double beta, f_double *y, f_int incy) {
    dgemv_(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy, 1);
}

int main()
{
  int i, j;

  for (i=0; i<3; ++i) {
    for (j=0; j<3; ++j) printf("%5.1f", m[i*3+j]);
    putchar('\n');
  }

  gemv('N', 3, 3, 1.0, m2.ptr, 3, x.ptr, 1, 0.0, y.ptr, 1);

  for (i=0; i<3; ++i)  printf("%5.1f\n", y[i]);

  writeln("
Above output is supposed to look like

  3.0  1.0  3.0
  1.0  5.0  9.0
  2.0  6.0  5.0
 -1.0
  3.0
 -3.0
");
  return 0;
}

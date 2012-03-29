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

void dgels_(char *trans, f_int *m, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *work, f_int *lwork, f_int *info);

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

extern (C) enum CBLAS_ORDER
    {CblasRowMajor=101, CblasColMajor=102 };
extern (C) enum CBLAS_TRANSPOSE
    {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113, AtlasConj=114};

const char CblasRowMajor = 101;
const char CblasNoTrans = 111;

extern(C) {

alias float f_float;
alias double f_double;
alias int f_int;
}

// From scid:
extern (C) void dgemv_(char *trans, f_int *m, f_int *n, f_double *alpha, f_double *A, f_int *lda, f_double *x, f_int *incx, f_double *beta, f_double *y, f_int *incy, f_int trans_len);


extern (C) void cblas_dgemv(int order, int TransA , 
  int M , int N , double alpha , f_double *A , int lda ,
  f_double *X , int incX , double beta , f_double *Y , int incY);

extern (C) void gemv(char *trans, f_int m, f_int n, f_double alpha, f_double *A, f_int lda,
f_double *x, f_int incx, f_double beta, f_double *y, f_int incy) {
    dgemv_(trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy, 1);
}

void dblas_dgemv(int order, int TransA , 
  int M , int N , double alpha , double A[] , int lda ,
  double X[] , int incX , double beta , double Y[] , int incY)
{
   writeln("order=",order);
   writeln("TransA=",TransA);
   writeln("M=",M);
   writeln("N=",N);
   writeln("X=",X);
   writeln("A=",A);
   writeln("Y=",Y);
   writeln("lda=",lda);
   writeln("incX=",incX);
   writeln("beta=",beta);
   writeln("incY=",incY);
   double *x = X.ptr;
   double *a = A.ptr;
   double *y = Y.ptr;
   writeln("calling into cblas_dgemv");
   // dgemv_(order,TransA,M,N,alpha,a,lda,x,incX,beta,y,incY);
   writeln("returned from cblas_dgemv");
   writeln(Y);
}

int main()
{
  int i, j;

  for (i=0; i<3; ++i) {
    for (j=0; j<3; ++j) printf("%5.1f", m[i*3+j]);
    putchar('\n');
  }

  // gemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, m, 3,
	//      x, 1, 0.0, y, 1);
  gemv(cast(char *)toStringz("N"), 3, 3, 1.0, m2.ptr, 3,
	     x.ptr, 1, 0.0, y.ptr, 1);

  for (i=0; i<3; ++i)  printf("%5.1f\n", y[i]);

  return 0;
}

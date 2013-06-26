/**
 * matrix functions for multiple QTL modelling and scanning
 **/

module qtl.core.mqm.matrix;

import std.c.stdlib;
import std.stdio;
import std.conv;

/*
 * Produces a N-dim vector
 */
pure T[] newvector(T)(size_t length, T value = T.init){
  T[] x;
  x.length = length;
  for(size_t j=0; j < length; j++){ x ~= value; }
  return x;
}

/*
 * Create a copy of a vector
 */
pure T[] copyvector(T)(in T[] c){
  T[] x;
  x.length = c.length;
  for(size_t j=0; j < c.length;j++){
    x[j]=c[j];
  }
  return x;
}

/*
 * Produces a R*C-dim matrix
 */
pure T[][] newmatrix(T)(size_t nrow, size_t ncol, T init = T.init){
  T[][] x;
  x.length=nrow;
  for(size_t i=0;i < nrow;i++){
    x[i] = newvector!T(ncol,init);
  }
  return x;
}

pure T[][] translate(T)(in T[][] i){
  if(i.length == 0) throw new Exception("Matrix needs to be at least of dim 1,1");
  T[][] m = newmatrix!T(i[0].length, i.length);
  for(size_t r=0;r<i.length;r++){
    for(size_t c=0;c<i[0].length;c++){
      m[c][r] = i[r][c];
    }
  }
  return m;
}



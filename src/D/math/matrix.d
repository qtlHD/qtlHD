/*
 * Matrix-related stuff
 */

module math.matrix;

import std.algorithm;
import std.math;
import std.range;
import std.stdio;
import std.string;
import std.conv;
import std.container;
import std.typecons;

import math.lapack.svd;

// convert doubly-index matrix to single-index thing
Tuple!(double[],size_t,size_t) matrix_as_vector(double [][]matrix)
{
  auto nrow = matrix.length;
  auto ncol = matrix[0].length;

  auto vector = new double[nrow*ncol];

  foreach(i, row; matrix)
    foreach(j, el; row)
      vector[j*nrow + i] = el;

  return tuple(vector, nrow, ncol);
}

unittest {
  writeln("Unit test " ~ __FILE__);

  auto x = [ [4.75, 3.12, 3.37, 3.81, 3.95],
             [4.89, 3.81, 2.79, 2.95, 3.63],
             [2.91, 2.40, 2.89, 3.39, 3.74],
             [2.63, 3.70, 3.05, 3.62, 2.72] ];

  auto xv = matrix_as_vector(x);

  assert(xv[1] == 4);
  assert(xv[2] == 5);
  assert(xv[0][5] == 3.81);
}

// select subset of columns of matrix
T [][] subset_columns(T)(T [][] matrix, size_t[] columns)
{
  auto nrow = matrix.length;
  auto ncol = matrix[0].length;
  foreach(col; columns)
    if(col < 0 || col >= ncol)
      throw new Exception("Invalid column " ~ to!string(col));

  auto result = new T[][](nrow, columns.length);
  foreach(i, row; matrix) {
    foreach(j, col; columns) {
      result[i][j] = row[col];
    }
  }

  return(result);
}

unittest {
  auto x = [ [4.75, 3.12, 3.37, 3.81, 3.95],
             [4.89, 3.81, 2.79, 2.95, 3.63],
             [2.91, 2.40, 2.89, 3.39, 3.74],
             [2.63, 3.70, 3.05, 3.62, 2.72] ];

  auto y = subset_columns(x, [0, 1, 4]);
  assert(y.length == 4);
  assert(y[0].length == 3);
  foreach(i, row; x) {
    assert(row[0] == y[i][0]);
    assert(row[1] == y[i][1]);
    assert(row[4] == y[i][2]);
  }
}

// select subset of columns of matrix
T [][] omit_column(T)(T [][] matrix, size_t column)
{
  auto nrow = matrix.length;
  auto ncol = matrix[0].length;
  if(column < 0 || column >= ncol)
    throw new Exception("column is outside allowable range");

  size_t[] columns;
  foreach(i; 0..ncol)
    if(i != column) columns ~= i;

  return subset_columns(matrix, columns);
}

unittest {
  auto x = [ [4.75, 3.12, 3.37, 3.81, 3.95],
             [4.89, 3.81, 2.79, 2.95, 3.63],
             [2.91, 2.40, 2.89, 3.39, 3.74],
             [2.63, 3.70, 3.05, 3.62, 2.72] ];

  auto y = omit_column(x, 3);
  assert(y.length == 4);
  assert(y[0].length == 4);
  foreach(i, row; x) {
    assert(row[0] == y[i][0]);
    assert(row[1] == y[i][1]);
    assert(row[2] == y[i][2]);
    assert(row[4] == y[i][3]);
  }
}

// for matrix of less-than-full rank, find a dependent column to omit
size_t find_dependent_column(double[][] matrix)
{
  auto ncol = matrix[0].length;
  auto rank = matrix_rank(matrix);

  if(rank == ncol) return(-2); // full rank

  foreach_reverse(column; 0..ncol) {
    auto newrank = matrix_rank(omit_column(matrix, column));
    if(newrank == rank) return(column);
  }

  return(-1); // can't find one
}

// for matrix of less-than-full rank, omit linearly dependent columns
double[][] omit_dependent_columns(double [][] matrix)
{
  auto ncol = matrix[0].length;
  auto rank = matrix_rank(matrix);

  if(rank == ncol) return(matrix);

  auto col_to_omit = find_dependent_column(matrix);
  if(col_to_omit == -1)
    throw new Exception("Trouble finding dependent column to omit");
  auto result = omit_column(matrix, col_to_omit);
  rank = matrix_rank(result);

  while(rank < result[0].length) {
    col_to_omit = find_dependent_column(result);
    if(col_to_omit == -1)
      throw new Exception("Trouble finding dependent column to omit");
    result = omit_column(result, col_to_omit);
    rank = matrix_rank(result);
  }

  return(result);
}

unittest {
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

  auto y = omit_dependent_columns(x);
  assert(y[0].length == 3);
}

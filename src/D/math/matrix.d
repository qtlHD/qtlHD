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
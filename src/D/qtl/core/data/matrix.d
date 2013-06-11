/**
 * Generic matrix functions 
 */

module qtl.core.data.matrix;

import std.conv;
import std.stdio;
import std.string;
import std.array;
import std.typecons;
import std.algorithm;

/**
 * Return those rows with index that match the test function on row elements (non-lazy)
 * FIXME: write test and use line operator
 */

auto filter_matrix_by_row_with_index(T)(T[][] matrix, bool function (T) test ) {
  uint i=0;
  return map!( (row) { 
    foreach(item; row) { 
      if (!test(item)) return null; 
      i += 1;
    }
    return tuple(i,row);
  })(matrix).remove(null).array();
}

/**
 * Return the rows matching the test function on row elements (non-lazy)
 * FIXME: write test and use line operator
 */

T[] filter_matrix_by_row(T)(T[][] matrix, bool function (T) test ) {
  return map!( (row) { 
    foreach(item; row) { if (!test(item)) return null; }
    return row;
  })(matrix).remove(null).array();
}

/**
 * Return a list of booleans for rows matching the test function for every item in the row
 * (non-lazy)
 */

bool[] filter_matrix_by_row_2bool(T)(T[][] matrix, bool function (T) test ) {
  return map!( (row) { return reduce!( (count,item) => count+test(item) )(0,row) > 0 ; } )(matrix).array();

}

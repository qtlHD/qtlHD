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
 * Return the rows matching the test function on row elements
 */

T[] filter_matrix_by_row(T)(T[][] matrix, bool function (T) test ) {
  return map!( (row) { 
    foreach(item; row) {
      if (!test(item))
        return null;
    }
    return row;
  })(matrix).array();
}

/**
 * Return a list of booleans for rows matching the test function for every item in the row
 */

bool[] filter_matrix_by_row_2bool(T)(T[][] matrix, bool function (T) test ) {
  return map!( (row) { return reduce!( (count,item) => count+test(item) )(0,row) > 0 ; } )(matrix).array();

}

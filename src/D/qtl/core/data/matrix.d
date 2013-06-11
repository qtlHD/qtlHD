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
 * Return a list of booleans for rows matching the test function for every item in the row
 */

bool[] filter_matrix_by_row_2bool(T)(T[][] matrix, bool delegate (T) test ) {
  return map!( (ind_p) { return reduce!( (count,p) => count+test(p) )(0,ind_p) > 0 ; } )(matrix).array();

}

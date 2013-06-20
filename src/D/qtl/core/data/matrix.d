/**
 * Generic matrix functions 
 *
 * These functions are considered standard algorithms for D and many of them
 * return lazy ranges. For laziness in D see http://ddili.org/ders/d.en/ranges.html
 * and Andrei's http://www.informit.com/articles/printerfriendly.aspx?p=1407357
 */

module qtl.core.data.matrix;

import std.conv;
import std.stdio;
import std.string;
import std.array;
import std.typecons;
import std.algorithm;

/** 
 * Test each row of a matrix and return the test result, the index and the row
 * element in a tuple. This is the full implementation that is lazy and can apply 
 * to almost all situations.
 */

auto test_matrix_by_row(T)(T[][] matrix, bool function (T) test ) {
  uint i=0;
  return map!( (row) { 
    foreach(col; row) { 
      if (!test(col)) return tuple(false,i,row); 
      i += 1;
    }
    return tuple(true,i,row);
  })(matrix);
}

/**
 * Return those rows with index that match the test function on row elements (non-lazy). 
 *
 * For each matching row a tuple is returned containing the row index, and a pointer to
 * the row array.
 */

Tuple!(uint, T[])[] filter_matrix_by_row_with_index(T)(T[][] matrix, bool function (T) test ) {
  auto range = filter!"a[0]"(
    test_matrix_by_row!double(matrix,test)
  );
  return map!( result => tuple(result[1],result[2]) )(range).array();
}

/**
 * Return the rows matching the test function on row elements (non-lazy)
 */

auto filter_matrix_by_row(T)(T[][] matrix, bool function (T) test ) {
  return remove!"a==[]"(map!( (row) { 
    foreach(item; row) { if (!test(item)) return null; }
    return row;
  })(matrix).array());
}

unittest {
  auto matrix = new double[][](3,3);
  matrix.length = 3;
  matrix[0] = [1.0,1.0,1.0];  
  matrix[1] = [1.0,0.0,1.0];  
  matrix[2] = [1.0,1.0,1.0];  
  auto rows1 = test_matrix_by_row!double(matrix, (item) => item==1.0);
  assert(rows1.length == 3,to!string(rows1.length));
  assert(rows1[0][0] == true,to!string(rows1[0][0]));
  assert(rows1[1][0] == false,to!string(rows1[1][0]));
  assert(rows1[2][0] == true,to!string(rows1[0][0]));
  auto rows = filter_matrix_by_row!double(matrix, (item) => item==1.0);
  assert(rows.length == 2,to!string(rows.length));
  auto tuples = filter_matrix_by_row_with_index!double(matrix, (item) => item==1.0);
  assert(tuples.length == 2,to!string(tuples.length));
}

/**
 * Return a list of booleans for rows matching the test function for every item in the row
 * (non-lazy)
 */

bool[] filter_matrix_by_row_2bool(T)(T[][] matrix, bool function (T) test ) {
  return map!( (row) {
    return reduce!( (count,item) => count+test(item) )(0,row) == row.length ; 
  } )(matrix).array();
}

unittest {
  auto matrix = new double[][](3,3);
  matrix.length = 3;
  matrix[0] = [1.0,1.0,1.0];  
  matrix[1] = [1.0,0.0,1.0];  
  matrix[2] = [1.0,1.0,1.0];  
  auto rows = filter_matrix_by_row_2bool!double(matrix, item => item==1.0);
  assert(rows == [true,false,true],to!string(rows));
}



/**
 * Generic matrix functions 
 *
 * These functions are considered standard algorithms for qtlHD and most of them
 * return lazy ranges. Some come with non-lazy alternatives to make type checking a 
 * little easier. 
 *
 * For laziness in D see
 *
 * http://ddili.org/ders/d.en/ranges.html and Andrei's
 * http://www.informit.com/articles/printerfriendly.aspx?p=1407357 and
 * https://www.semitwist.com/articles/article/view/combine-coroutines-and-input-ranges-for-dead-simple-d-iteration
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

auto test_matrix_by_row(T)(T[][] matrix, bool function (T[]) test_row ) {
  uint i=0;
  return map!( (row) { 
    i += 1;
    return tuple(test_row(row),i-1,row); 
  })(matrix);
}

auto test_matrix_by_row_element(T)(T[][] matrix, bool function (T) test_element ) {
  // return test_matrix_by_row!T(matrix, row => true );
  // return test_matrix_by_row!T(matrix, (row) { return true; } );
  /* This should work, but does not
  return test_matrix_by_row!T(matrix, (row) { 
    foreach(element ; row) {
      if (!test2(element)) return false;
    }
    return true;
  });
  */
  uint i=0;
  return map!( (row) { 
    foreach(col; row) { 
      if (!test_element(col)) return tuple(false,i,row); 
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

auto filter_matrix_by_row_with_index(T)(T[][] matrix, bool function (T) test ) {
  // filter on all elements that returned false
  auto range = filter!"a[0]"( test_matrix_by_row_element!double(matrix,test) );
  return map!( result => tuple(result[1],result[2]) )(range);
}

Tuple!(uint, T[])[] non_lazy_filter_matrix_by_row_with_index(T)(T[][] matrix, bool function (T) test ) {
  return array(filter_matrix_by_row_with_index!T(matrix,test));
}

/**
 * Return the rows matching the test function on row elements (non-lazy)
 */

T[][] non_lazy_filter_matrix_by_row(T)(T[][] matrix, bool function (T) test ) {
  // the first element of each tuple is the bool
  auto range = filter!"a[0]"( test_matrix_by_row_element(matrix,test) );
  // the third element is the row
  return map!( result => result[2] )(range).array();
}

unittest {
  auto matrix = new double[][](3,3);
  matrix.length = 3;
  matrix[0] = [1.0,1.0,1.0];  
  matrix[1] = [1.0,0.0,1.0];  
  matrix[2] = [1.0,1.0,1.0];  
  // Test the matrix by row element
  auto rows1 = test_matrix_by_row_element!double(matrix, (item) => item==1.0);
  assert(rows1.length == 3,to!string(rows1.length));
  assert(rows1[0][0] == true,to!string(rows1[0][0]));
  assert(rows1[1][0] == false,to!string(rows1[1][0]));
  assert(rows1[2][0] == true,to!string(rows1[0][0]));
  auto rows = non_lazy_filter_matrix_by_row!double(matrix, (item) => item==1.0);
  assert(rows.length == 2,to!string(rows.length));
  assert(rows[1] == [1.0,1.0,1.0]);
  auto tuples = non_lazy_filter_matrix_by_row_with_index!double(matrix, (item) => item==1.0);
  assert(tuples.length == 2,to!string(tuples.length));
  assert(tuples[1][1] == [1.0,1.0,1.0]);
  // Test the matrix by row 
  auto rows2 = test_matrix_by_row!double(matrix, r => true );
  assert(rows2.length == 3,to!string(rows2.length));
  auto rows3 = test_matrix_by_row!double(matrix, (r) { 
    return (reduce!"a+(b==1.0)"(0,r) == r.length); }
  );
  assert(rows3.length == 3,to!string(rows3.length));
  assert(rows3[0][0] == true,to!string(rows3[0][0]));
  assert(rows3[1][0] == false,to!string(rows3[1][0]));
  assert(rows3[2][0] == true,to!string(rows3[0][0]));
}

/**
 * Return a list of booleans for rows matching the test function for every item in the row
 * (non-lazy)
 */

bool[] filter_matrix_by_row_2bool(T)(T[][] matrix, bool function (T) test ) {
  // the first element of each tuple is the bool
  return map!"a[0]"( test_matrix_by_row_element(matrix,test)).array();
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



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
import std.range;

/** 
 * Test each row of a matrix and return the test result, the index and the row
 * element in a tuple. This is the full implementation that is lazy and can apply 
 * to almost all situations.
 */

auto test_matrix_by_row(T)(T[][] matrix, bool function (T[]) test_row ) {
  uint i=0;
  return map!( (row) { 
    return tuple(test_row(row),i++,row); 
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
    i += 1;
    foreach(element; row) { 
      // breaks early when element test is false (no need to test the rest)
      if (!test_element(element)) return tuple(false,i-1,row); 
    }
    return tuple(true,i-1,row);
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
  // writeln(matrix);
  // Test the matrix by row element
  auto rows1 = test_matrix_by_row_element!double(matrix, (item) => item==1.0).array();
  assert(rows1.length == 3,to!string(rows1.length));
  writeln(rows1);
  assert(rows1[0][0] == true,to!string(rows1[0][0]));
  assert(rows1[0][1] == 0,to!string(rows1[0][1]));
  assert(rows1[1][0] == false,to!string(rows1[1][0]));
  assert(rows1[2][0] == true,to!string(rows1[2][0]));
  assert(rows1[2][1] == 2,to!string(rows1[2][1]));
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
  ).array();
  assert(rows3.length == 3,to!string(rows3.length));
  assert(rows3[0][0] == true,to!string(rows3[0][0]));
  assert(rows3[1][0] == false,to!string(rows3[1][0]));
  assert(rows3[2][0] == true,to!string(rows3[2][0]));
  assert(rows3[2][1] == 2,to!string(rows3[2][1]));
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

/** 
 * Test each column of a matrix and return the test result, the index and the column
 * element in a tuple. This is the full implementation that is lazy and can apply 
 * to almost all situations.
 *
 * To create a lazy return value we loop with map. The column is returned as a 
 * concrete array (if we make that lazy, it can get kinda crazy!)
 */

auto gen_t_matrix(T)(T[][] matrix) {
  assert(matrix.length > 0);
  uint i = 0;
  return map!( column => transversal(matrix,i++).array() )(matrix[0]);
}

auto test_matrix_by_col(T)(T[][] matrix, bool function (T[]) test_col ) {
  uint i = 0;
  return map!( (col) { 
    return tuple(test_col(col),i++,col); 
  })(gen_t_matrix(matrix));
}

auto test_matrix_by_col_element(T)(T[][] matrix, bool function (T) test_element ) {
  uint i=0;
  return map!( (col) { 
    i += 1;
    foreach(element; col) { 
      // breaks early when element test is false (no need to test the rest)
      if (!test_element(element)) return tuple(false,i-1,col); 
    }
    return tuple(true,i-1,col);
  })(gen_t_matrix(matrix));
}

T[][] non_lazy_filter_matrix_by_col(T)(T[][] matrix, bool function (T) test ) {
  // the first element of each tuple is the bool - at this point the result 
  // is concrete, so as to not invoke the lambdas again. Maybe this can be
  // improved...
  auto range = filter!"a[0]"( test_matrix_by_col_element(matrix,test).array() );
  // the third element is the col, and now make concrete...
  return map!( result => result[2] )(range).array();
}

auto filter_matrix_by_col_with_index(T)(T[][] matrix, bool function (T) test ) {
  // filter on all elements that returned false
  auto range = filter!"a[0]"( test_matrix_by_col_element!double(matrix,test).array() );
  return map!( result => tuple(result[1],result[2]) )(range);
}

Tuple!(uint, T[])[] non_lazy_filter_matrix_by_col_with_index(T)(T[][] matrix, bool function (T) test ) {
  return array(filter_matrix_by_col_with_index!T(matrix,test));
}


unittest {
  writeln("Start column-wise");
  auto matrix = new double[][](3,3);
  matrix.length = 3;
  matrix[0] = [1.0,1.0,1.0];  
  matrix[1] = [1.0,0.0,1.0];  
  matrix[2] = [1.0,2.0,1.0];  
  // writeln(matrix);
  // Test the matrix by col element (all items in column should be 1.0)
  auto cols1 = test_matrix_by_col_element!double(matrix, (item) => item==1.0).array();
  // lazy version returns information on 3 columns
  assert(cols1.length == 3,to!string(cols1.length));
  assert(cols1[0][0] == true,to!string(cols1[0][0]));
  assert(cols1[0][1] == 0,to!string(cols1[0][1]));
  assert(cols1[1][0] == false,to!string(cols1[1][0]));
  assert(cols1[2][0] == true,to!string(cols1[2][0]));
  assert(cols1[2][1] == 2,to!string(cols1[2][1]));
  // Now filter down to 2 columns
  auto cols = non_lazy_filter_matrix_by_col!double(matrix, item => item==1.0);
  writeln("##",cols);
  assert(cols.length == 2,to!string(cols.length));
  assert(cols[1] == [1.0,1.0,1.0]);
  auto tuples = non_lazy_filter_matrix_by_col_with_index!double(matrix, (item) => item==1.0);
  assert(tuples.length == 2,to!string(tuples.length));
  assert(tuples[1][1] == [1.0,1.0,1.0]);
  // Test the matrix by col 
  auto cols2 = test_matrix_by_col!double(matrix, r => true ).array();
  assert(cols2.length == 3,to!string(cols2.length));
  auto cols3 = test_matrix_by_col!double(matrix, (r) { 
    return (reduce!"a+(b==1.0)"(0,r) == r.length); }
  ).array();
  assert(cols3.length == 3,to!string(cols3.length));
  writeln(cols3);
  assert(cols3[0][0] == true,to!string(cols3[0][0]));
  assert(cols3[1][0] == false,to!string(cols3[1][0]));
  assert(cols3[2][0] == true,to!string(cols3[2][0]));
  assert(cols3[2][1] == 2,to!string(cols3[2][1]));
}

/**
 * pull out a column of a matrix indexed as matrix[row][col]
 **/
T[] get_column(T)(in size_t column_index, T[][] matrix)
{
  if(column_index < 0 || column_index >= matrix[0].length)
    throw new Exception("column_index outside of allowable range");

  return map!( row => row[column_index] )(matrix).array();
}


/**
 * transpose of a matrix
 */
T[][] transpose_matrix(T)(in T[][] matrix)
{
  auto ret = new T[][](matrix[0].length, matrix.length);

  foreach(i, row; matrix) {
    foreach(j, element; row)
      ret[j][i] = element;
  }

  return ret;
}

import qtl.core.phenotype;

unittest {
  writeln("Test get_column and transpose_matrix");

  Phenotype pheno[][];
  auto pdbl = [ ["0.0", "0.0", "0.0", "0.0", "0.0"],
                ["0.0", "0.0", "0.0", "0.0", "0.0"],
                [  "-", "0.0", "0.0",   "-", "0.0"],
                ["0.0", "0.0", "0.0", "0.0", "0.0"],
                ["0.0", "0.0", "0.0", "0.0", "0.0"],
                ["0.0", "0.0", "0.0", "0.0", "0.0"],
                [  "-",   "-",   "-",   "-",   "-"],
                ["0.0", "0.0", "0.0", "0.0", "0.0"],
                [  "-",   "-", "0.0",   "-",   "-"],
                ["0.0", "0.0", "0.0", "0.0", "0.0"],
                ["0.0",   "-",   "-", "0.0",   "-"],
                ["0.0",   "-",   "-", "0.0",   "-"],
                ["0.0", "0.0", "0.0", "0.0", "0.0"]];

  foreach(line; pdbl) {
    Phenotype[] ps = std.array.array(map!((a) {return
    set_phenotype(a);})(line));
    pheno ~= ps;
  }

  writeln("Test get_column");
  auto pcol = get_column(1, pheno);
  foreach(i, row; pheno)
    assert(isSame(row[1], pcol[i]));

  writeln("Test transpose of phenotype matrix");
  auto tpheno = transpose_matrix(pheno);

  foreach(i, row; tpheno)
    foreach(j, element; row)
      assert(isSame(element, pheno[j][i]));
}

module qtl.core.matrices;

//import core.memory;
import std.c.stdlib;
import std.stdio;
import std.conv;

/*
 * Produces an array with start,start+1,start+2 .. start+length-1
 */
pure int[] doRange(int start,int length){
  int array[];
  for(int i = 0; i < (length-1); i++){
   array ~= start+i;
  }
  return array;
}

/*
 * Produces an array of size length, filled with value
 */
pure T[] doArray(T)(int length,T value){
  T array[];
  for(int i = 0; i < (length-1); i++){
   array ~= value;
  }
  return array;
}

/*
 * Produces a R*C-dim matrix
 */
T** newmatrix(T)(size_t rows, size_t cols) {
  T** m;
  m = cast(T**)calloc(rows, (T*).sizeof);
  if(m is null){
    writeln("Not enough memory for new matrix");
  }
  for(size_t i=0; i<rows; i++) {
    m[i]= newvector!T(cols);
  }
  return m;
}

/*
 * Prints a matrix to stdout
 */
void printmatrix(T)(T** m, int rows, int cols) {
  for (int r=0; r<rows; r++) {
    for (int c=0; c<cols; c++) {
      writeln("%s",to!string(m[r][c]));
    }
  }
}

/*
 * Produces a N-dim vector
 */
T* newvector(T)(size_t dim) {
  T* v;
  v = cast(T*)calloc(dim, T.sizeof);
  if(v is null){
    writeln("Not enough memory for new vector of dimension %d",(dim+1));
  }
  return v;
}

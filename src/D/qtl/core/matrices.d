module qtl.core.matrices;

import core.memory;
import std.stdio;
import std.conv;

T** newmatrix(T)(int rows, int cols) {
  T** m;
  m = cast(T**)GC.malloc(rows, T.sizeof);
  if(m is null){
    writeln("Not enough memory for new matrix");
  }
  for(int i=0; i<rows; i++) {
    m[i]= newvector!T(cols);
  }
  return m;
}

void printmatrix(T)(T** m, int rows, int cols) {
  for (int r=0; r<rows; r++) {
    for (int c=0; c<cols; c++) {
      writeln("%f",m[r][c]);
    }
    writeln("col done");
  }
}

T* newvector(T)(int dim) {
  T* v;
  v = cast(T*)GC.malloc(dim, T.sizeof);
  if(v is null){
    writeln("Not enough memory for new vector of dimension %d",(dim+1));
  }
  return v;
}

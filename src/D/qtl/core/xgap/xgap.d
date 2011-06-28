/**
 * XGap types 
 **/

module qtl.core.xgap.xgap;

import std.file;
import std.conv;
 
enum MatrixType : uint { 
  EMPTY = 0, 
  INTMATRIX = 1, 
  DOUBLEMATRIX = 2, 
  FIXEDCHARMATRIX = 3, 
  VARCHARMATRIX = 4
};

immutable byte[2] b_footprint = [ 0, 5 ];
immutable byte[3] b_version = [ 0, 0, 1 ];

/*
 * Helper function to go from sizes in bytes to KBs
 */
double toKb(in string filename){
  return cast(double) getSize(filename)/1024;
}

/*
 * Helper function to go from ubyte[] to int
 */
int byteToInt(ubyte[] bits, bool little_endian = true ){
  return *cast(int*)bits;
}

/*
 * Helper function to go from ubyte[] to double
 */  
double byteToDouble(ubyte[] bits){
  return *cast(double*)bits;
}

/*
 * Helper function to go from ubyte[] to string
 */  
string byteToString(ubyte[] bits){
  char[] r;
  foreach(ubyte b;bits){
    r ~= cast(char)b;
  }
  r ~= '\0';
  return to!string(r);
}

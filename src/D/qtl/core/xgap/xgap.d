/**
 * XGap types 
 **/

module qtl.core.xgap.xgap;

import std.file;
import std.conv;

//Aliasses to define Footprint and Version
alias byte[8] Footprint;
alias byte[4] Version;
//Aliasses to more easily refer to the different matrices
alias XgapMatrixData!int     integermatrix;
alias XgapMatrixData!double  doublematrix;
alias XgapMatrixData!string  stringmatrix;
//Current Footprint and Version
immutable Footprint xgap_footprint = [ 0, 'X', 'G', 'A', 'P', 'B', 0, 1 ];
immutable Version   xgap_version   = [ 0, 0, 1, 'A' ];

enum MatrixType : uint { 
  EMPTY = 0, 
  INTMATRIX = 1, 
  DOUBLEMATRIX = 2, 
  FIXEDCHARMATRIX = 3, 
  VARCHARMATRIX = 4
};

enum MatrixClass : uint { 
  EMPTY = 0,
  PHENOTYPE = 1, 
  GENOTYPE = 2, 
  MAP = 3,
  ANNOTATION = 4
};

/*
 * XgapBinary file header
 */
struct XgapFileHeader{
  Footprint   magicn = xgap_footprint;
  Version     fileversion = xgap_version;
  byte[5]     p0 = [0,0,0,0,0];
  int         nmatrices;
  byte[4]     p1 = [0,0,0,0];
}

interface Container {

}

//XgapBinary matrix structure holding a lengths and templated data
class XgapMatrixData(T) : Container{
  int[]         lengths;
  T[][]         data;
}

//XgapBinary matrix structure holding a lengths and templated data
struct XgapMatrixHeader{
  Footprint     magicn = xgap_footprint;
  
  MatrixType    type;
  MatrixClass   mclass;
  
  int           size;      //former skip
  int           nrow;

  int           ncol;
  int           start;     //start of the matrix in the file
}

class XgapMatrix {
  XgapMatrixHeader  header;
  Container data;
}

/*
 * Helper function to go from sizes in bytes to KBs
 */
double toKb(in string filename){
  return cast(double) getSize(filename)/1024;
}

/*
 * Helper function to go from ubyte[] to int
 */
int byteToInt(ubyte[] bits){
  return convbyte!int(bits);
}

/*
 * Helper function to go from ubyte[] to double
 */  
double byteToDouble(ubyte[] bits){
  return convbyte!double(bits);
}

T convbyte(T)(ubyte[] bits){
  return *cast(T*)bits;
}

/*
 * Helper function to go from ubyte[] to string
 */  
string byteToString(ubyte[] bits){
  char[] r;
  foreach(ubyte b;bits){
    r ~= convbyte!char([b]);
  }
  r ~= '\0';
  return to!string(r);
}

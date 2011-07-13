/**
 * XGap types 
 **/

module qtl.core.xgap.xgap;

import std.file;
import std.conv;

import qtl.core.primitives;

//Aliasses to define MagicNumber and Version
alias byte[8] MagicNumber;
alias byte[4] Version;
//Aliasses to more easily refer to the different matrices
alias XgapMatrixData!int     IntegerMatrix;
alias XgapMatrixData!double  DoubleMatrix;
alias XgapMatrixData!string  StringMatrix;
//Current MagicNumber and Version
immutable MagicNumber xgap_magicnumber = [ 0, 'X', 'G', 'A', 'P', 'B', 0, 1 ];
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

//XgapBinary file header, see doc/input/XGap.md
struct XgapFileHeader{
  MagicNumber   magicn = xgap_magicnumber;
  Version     fileversion = xgap_version;
  int         nmatrices;
}

//Container interface
interface Container {
}

//Xgap matrix data structure holding a lengths and templated data
class XgapMatrixData(T) : Container{
  int[]         lengths;
  T[][]         data;
}

//Xgap matrix header structure, see doc/input/XGap.md
struct XgapMatrixHeader{
  MagicNumber     magicn = xgap_magicnumber;
  
  MatrixType    type;
  MatrixClass   mclass;
  
  int           size;      //former skip
  int           nrow;

  int           ncol;
  byte[4]       pad16b;
}

//A matrix is defined as a header with data
class XgapMatrix {
  XgapMatrixHeader  header;
  Container data;
}

//Helper function to go from marker[] to a 3 column format
string[][] getMarkerInfoMatrix(Marker[] markers){
  string[][] return_matrix;
  foreach(Marker m;markers){
    string[] row;
    row ~= m.name;
    row ~= m.chromosome.name;
    row ~= to!string(m.position);
    return_matrix ~= row;
  }
  return return_matrix;
}

//Helper function to go from sizes in bytes to KBs
 double toKb(in string filename){
  return cast(double) getSize(filename)/1024;
}

//Helper function to go from ubyte[] to int
int byteToInt(ubyte[] bits){
  return convbyte!int(bits);
}

//Helper function to go from ubyte[] to double 
double byteToDouble(ubyte[] bits){
  return convbyte!double(bits);
}

//Helper function to go from ubyte[] to any class
T convbyte(T)(ubyte[] bits){
  return *cast(T*)bits;
}

//Helper function to go from ubyte[] to string
string byteToString(ubyte[] bits){
  char[] r;
  foreach(ubyte b;bits){
    r ~= cast(char)(b);
  }
  //r ~= '\0';
  return to!string(r);
}

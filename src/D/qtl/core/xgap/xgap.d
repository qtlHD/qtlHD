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

/*
 * Different types of data types we allow for storage
 */
enum MatrixType : uint { 
  EMPTY = 0, 
  INTMATRIX = 1, 
  DOUBLEMATRIX = 2, 
  FIXEDCHARMATRIX = 3, 
  VARCHARMATRIX = 4
};

/*
 * We tag our data as being of a certain 'Class'
 * This determines how data can be used by the algorithms
 */
enum MatrixClass : uint { 
  EMPTY = 0,
  PHENOTYPE = 1, 
  GENOTYPE = 2, 
  MAP = 3,
  ANNOTATION = 4
};

/* XgapBinary file header
 * This is the first structure we encounter in an XGAP binary file
 * We start by printing a magicnumber followed 
 * by the version and the number of matrices stored
 * for more information see: doc/input/XGap.md
 */
struct XgapFileHeader{
  MagicNumber   magicn = xgap_magicnumber;
  Version     fileversion = xgap_version;
  int         nmatrices;
}

//Container interface
interface Container {
}

/*
 * Xgap matrix data structure 
 * This object stores the length of all the elements inside.
 * Also this objects contains the templated data.
 */
class XgapMatrixData(T) : Container{
  int[]         lengths;
  T[][]         data;
}

/*
 * Xgap matrix header structure
 * Each data field (representing a matrix) inside the xgap binary file needs a header.
 * This header holds: 
 *  - A magic number (for file integrity checks)
 *  - Matrix type and class
 *  - size: The number of elements inside the matrix
 *  - nrow: Number of matrix rows
 *  - ncol: Number of matrix columns
 *  - pad16b: Padding for the header
 * For more information see: doc/input/XGap.md
 */
struct XgapMatrixHeader{
  MagicNumber   magicn = xgap_magicnumber;
  
  MatrixType    type;
  MatrixClass   mclass;
  
  int           size;      //former skip
  int           nrow;

  int           ncol;
  byte[4]       pad16b;
}

struct XgapMatrixNames{
  int[]         rowlengths;
  string[]      rownames;
  int[]         collengths;
  string[]      colnames;
}

//A matrix is defined as a header with data
class XgapMatrix {
  XgapMatrixHeader  header;
  XgapMatrixNames   names;
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

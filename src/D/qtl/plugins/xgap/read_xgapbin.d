/**
 * Read binary XGap format
 *
 **/
 
module qtl.core.xgap.read_xgapbin;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.deprecate.genotype_enum;
import qtl.core.xgap.xgap;

import qtl.plugins.csvr.read_csvr;
import qtl.plugins.xgap.write_xgapbin;

import std.stdio;
import std.typecons;
import std.conv;
import std.string;
import std.path;
import std.file;

class XbinReader {
  XgapFileHeader      header;         //Binary file header
  XgapMatrixHeader[]  headers;        //Matrix headers
  XgapMatrixNames[]   names;          //Row and column names for the matrices
  private File    f;              //File pointer
  ubyte[]         inputbuffer;    //Buffered file content
  bool            correct;        //Is the file correct?

  //Checks if the buffer is equal to the magic number
  bool checkMagicNumber(in ubyte[] buffer){
    if(xgap_magicnumber == buffer) return true;
    return false;
  }
  
  //Checks to see if the file is loaded into memory.
  //We can check this by checking the magic numbers at the beginning and ending
  bool checkBuffer(ubyte[] buffer){
    ubyte[] startprint = buffer[0..MagicNumber.sizeof];
    ubyte[] endprint = buffer[buffer.length-MagicNumber.sizeof..$];
    if(checkMagicNumber(startprint) && checkMagicNumber(endprint)) return (correct = true);
    return (correct=false);
  }
  
  T[] toType(T)(ubyte[] buffer){
    T[] returnbuffer;
    foreach(int i, byte b ; buffer){
      returnbuffer ~= to!T(b);
    }
    return returnbuffer;
  }
  
  //CInternal function to load and cast data to the types specified by the header of the matrix
  T[][] loadData(T)(XgapMatrixHeader h, int[] lengths, int start){
    T[][] data;
    int skip=start;
    for(int r=0;r<h.nrow;r++){
      T[] row;
      for(int c=0;c<h.ncol;c++){
        int s = lengths[(r*h.ncol)+c];
        if(h.type == MatrixType.FIXEDCHARMATRIX || h.type == MatrixType.VARCHARMATRIX){
          row ~= byteToString(inputbuffer[skip..skip+s]);
        }else{
          row ~= convbyte!T(inputbuffer[skip..skip+s]);
        }
        skip += s;
      }
      data ~= row;
    }
    return data;
  }
  
  XgapMatrixNames getNames(int start, XgapMatrixHeader header){
    XgapMatrixNames data;
    for(int r=0;r<header.nrow;r++){
      data.rowlengths ~= byteToInt(inputbuffer[(start+(r*int.sizeof))..(start + int.sizeof + (r*int.sizeof))]);
    }
    auto skip = start + (header.nrow*int.sizeof);
    for(int r=0;r<header.nrow;r++){
      if(data.rowlengths[r] != 0){
        data.rownames ~= byteToString(inputbuffer[skip..skip+data.rowlengths[r]]);
      }else{
        data.rownames ~= "";
      }
      skip += data.rowlengths[r];
    }
    for(int c=0;c<header.ncol;c++){
      data.collengths ~= byteToInt(inputbuffer[(skip+(c*int.sizeof))..(skip + int.sizeof + (c*int.sizeof))]);
    }
    skip = skip + (header.ncol*int.sizeof);
    for(int c=0;c<header.ncol;c++){
      if(data.collengths[c] != 0){
        data.colnames ~= byteToString(inputbuffer[skip..skip+data.collengths[c]]);
      }else{
        data.colnames ~= "";
      }
      skip += data.collengths[c];
    }
    data.size = cast(int)(skip - start);
    return data;
  }
  
 /*
  * Loads the matrixid'th XgapMatrix from the file
  */
  XgapMatrix load(int matrixid = 0){
    if(!(matrixid >= 0 && matrixid < headers.length)){
      throw new Exception("No such matrix");
    }
    int skip = XgapFileHeader.sizeof;
    for(int i=0;i < matrixid;i++){
      skip += headers[i].size;
      skip += names[i].size;
    }
    XgapMatrix returnmatrix = new XgapMatrix();
    XgapMatrixHeader h = headers[matrixid];
    returnmatrix.header = h;
    skip += XgapMatrixHeader.sizeof;
  //*  XgapMatrixNames names;
    returnmatrix.names = getNames(skip,h);
    skip += returnmatrix.names.size;

    int start; //In matrix location of start fo the data

    int[] lengths;
    for(int x=0;x<(h.nrow*h.ncol);x++){
      if(h.type == MatrixType.VARCHARMATRIX){
        lengths ~= byteToInt(inputbuffer[(skip+(x*int.sizeof))..(skip + int.sizeof + (x*int.sizeof))]);
        start = cast(int) (skip + ((h.nrow*h.ncol)*int.sizeof));      
      }else{
        lengths ~= byteToInt(inputbuffer[skip..(skip+int.sizeof)]);
        start = cast(int) (skip + int.sizeof);
      }
    }
    switch(h.type){
      case MatrixType.INTMATRIX:
        IntegerMatrix c = new IntegerMatrix();
        c.lengths = lengths;
        c.data = loadData!int(h, lengths, start);
        returnmatrix.data = c;
      break;        
      case MatrixType.DOUBLEMATRIX:
        DoubleMatrix c = new DoubleMatrix();
        c.lengths = lengths;
        c.data = loadData!double(h, lengths, start);
        returnmatrix.data = c;
      break;        
      case MatrixType.FIXEDCHARMATRIX:
        StringMatrix c = new StringMatrix();
        c.lengths = lengths;
        c.data = loadData!string(h, lengths, start);
        returnmatrix.data = c;
      break;        
      case MatrixType.VARCHARMATRIX:
        StringMatrix c = new StringMatrix();
        c.lengths = lengths;
        c.data = loadData!string(h, lengths, start);
        returnmatrix.data = c;
      break;
      default:
        throw new Exception("Trying to load unsupported matrix type format");
      break;             
    }
    return returnmatrix;
  }
  
  //Get the version number of the file
  Version getVersion(){
    return header.fileversion;
  }
  
  //Return the number of matrices inside the file
  int getNumberOfMatrices(){
    return header.nmatrices;
  }
 
 /*
  *Parses the matrix header, reading type, dimensions and element sizes
  */
  int parseMatrixHeader(ubyte[] buffer, int matrix, int start){
    XgapMatrixHeader header = convbyte!(XgapMatrixHeader)(inputbuffer[start..start+XgapMatrixHeader.sizeof]);
    headers ~= header;
    XgapMatrixNames matrixnames = getNames(cast(int)(start+XgapMatrixHeader.sizeof),header);
    names ~= matrixnames;
    writefln("      Matrix %d, type=%d, rows: %d, columns: %d, %d %d", matrix, header.type, header.nrow, header.ncol, header.size, matrixnames.size);
    return header.size + matrixnames.size;
  }
  
 /*
  *Parses the file header, versions and # of matrices
  */
  bool parseFileHeader(bool verbose){
    writeln("    File OK? " ~ to!string(checkBuffer(inputbuffer)));
    header = convbyte!(XgapFileHeader)(inputbuffer[0..XgapFileHeader.sizeof]);
    writeln("    Version: " ~ to!string(header.fileversion));
    writeln("    Matrices: " ~ to!string(header.nmatrices));
    assert(correct,"Footprint failed");
    return correct;
  }
  
  //Read the entire file to memory and parse the matrixheaders
  this(in string filename,in bool verbose = false){
    assert(getSize(filename) < uint.max);
    inputbuffer = new ubyte[cast(uint)getSize(filename)];
    auto f = new File(filename,"rb");
    f.rawRead(inputbuffer);
    //Header
    parseFileHeader(verbose);

    //Loop through the matrices
    int skip = XgapFileHeader.sizeof;
    for(int m=0; m<getNumberOfMatrices(); m++){
      skip += parseMatrixHeader(inputbuffer, m, skip);
    }
  }
}

unittest{
  writeln("Unit test " ~ __FILE__);
  writeln("  - writing XBIN ");
  alias std.path.buildPath buildPath;
  auto infn = to!string(dirName(__FILE__) ~ sep ~ buildPath("..","..","..","..","..","test","data","input","multitrait.csvr"));
  auto outfn = to!string(dirName(__FILE__) ~ sep ~ buildPath("..","..","..","..","..","test","data","input","multitrait.xbin"));
  writeln("  - reading CSVR " ~ infn ~" to " ~ outfn);
  auto indata = new CSVrReader!RIL(infn);
  auto result = new BinaryWriter!(CSVrReader!RIL,RIL)(indata,outfn);
  writefln("Size (txt to xbin): (%.2f Kb to %.2f Kb)", toKb(infn), toKb(outfn));
  
  writeln("  - reading XBIN " ~ outfn);
  auto data = new XbinReader(outfn);
  assert(data.correct == true);
  assert(data.header.nmatrices == 3);
  assert(data.header.fileversion == [0,0,1, 'A']);
  XgapMatrix m1 = data.load(0);
  assert(set_phenotype!double(to!string((cast(DoubleMatrix)(m1.data)).data[0][1])) == indata.phenotypes[0][1]);
  assert(set_phenotype!double(to!string((cast(DoubleMatrix)(m1.data)).data[5][10])) == indata.phenotypes[5][10]);
  assert(set_phenotype!double(to!string((cast(DoubleMatrix)(m1.data)).data[10][5])) == indata.phenotypes[10][5]);
  writeln("   - Reloaded phenotypes from xgap binary");
  XgapMatrix m2 = data.load(1);
  assert(set_genotype!RIL((cast(StringMatrix)(m2.data)).data[3][1]) == indata.genotypes[3][1]);
  assert(set_genotype!RIL((cast(StringMatrix)(m2.data)).data[7][15]) == indata.genotypes[7][15]);
  assert(set_genotype!RIL((cast(StringMatrix)(m2.data)).data[20][3]) == indata.genotypes[20][3]);
  writeln("   - Reloaded genotypes from xgap binary");
  XgapMatrix m3 = data.load(2);
  assert((cast(StringMatrix)(m3.data)).data[1][0] == getMarkerInfoMatrix(indata.markers)[1][0]);
  assert((cast(StringMatrix)(m3.data)).data[10][2] == getMarkerInfoMatrix(indata.markers)[10][2]);
  assert((cast(StringMatrix)(m3.data)).data[20][1] == getMarkerInfoMatrix(indata.markers)[20][1]);  
  writeln("   - Reloaded markers from xgap binary");
}


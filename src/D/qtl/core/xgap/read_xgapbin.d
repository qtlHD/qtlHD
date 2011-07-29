/**
 * Write binary XGap format
 **/
 
module qtl.core.xgap.read_xgapbin;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.core.xgap.xgap;

import qtl.plugins.input.read_csvr;
import qtl.plugins.output.write_xgapbin;

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.file;

class XbinReader {
  XgapFileHeader      header;         //Binary file header
  XgapMatrixHeader[]  headers;        //Matrix headers
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
    }
    XgapMatrix returnmatrix = new XgapMatrix();
    XgapMatrixHeader h = headers[matrixid];
    returnmatrix.header = h;
    skip += XgapMatrixHeader.sizeof;

    int[] lengths;
    int start; //In matrix location of start fo the data
    for(int x=0;x<(h.nrow*h.ncol);x++){
      if(h.type == MatrixType.VARCHARMATRIX){
        lengths ~= byteToInt(inputbuffer[(skip+(x*int.sizeof))..(skip + int.sizeof + (x*int.sizeof))]);
        start = (skip + ((h.nrow*h.ncol)*int.sizeof));      
      }else{
        lengths ~= byteToInt(inputbuffer[skip..(skip+int.sizeof)]);
        start = (skip + int.sizeof);
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
    writefln("      Matrix %d, type=%d, rows: %d, columns: %d", matrix, header.type, header.nrow, header.ncol);
    headers ~= header;
    return header.size;
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
  alias std.path.join join;
  auto infn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","multitrait.csvr");
  auto outfn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","multitrait.xbin");
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


/**
 * Write binary XGap format
 **/
 
module qtl.core.xgap.read_xgapbin;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.core.xgap.xgap;

import qtl.plugins.input.read_interface;
import qtl.plugins.input.read_csvr;
import qtl.plugins.output.write_xgapbin;

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.file;

class XbinReader(XType) : GenericReader!XType{
  XgapBinHeader   header;         //Binary file header
  MatrixHeader[]  headers;        //Matrix headers
  int[]           skips;          //Stores the matrix skips, to get to the start of a matrix
  private File    f;              //File pointer
  ubyte[]         inputbuffer;    //Buffered file content
  bool            correct;        //Is the file correct?
  
  Matrix[]        matrices;       //MatrixHeaders
  
  bool checkFootprint(in ubyte[] buffer){
    if(xgap_footprint == buffer) return true;
    return false;
  }
  
  bool checkBuffer(ubyte[] buffer){
    ubyte[] startprint = buffer[0..Footprint.sizeof];
    ubyte[] endprint = buffer[buffer.length-Footprint.sizeof..$];
    if(checkFootprint(startprint) && checkFootprint(endprint)) return (correct = true);
    return (correct=false);
  }
  
  T[] toType(T)(ubyte[] buffer){
    T[] returnbuffer;
    foreach(int i, byte b ; buffer){
      returnbuffer ~= to!T(b);
    }
    return returnbuffer;
  }
  
 /*
  * Loads the ith matrix of MatrixClasstype from the file
  */
  void load(MatrixClass, int i = 0){
  
  }
  
  Version getVersion(){
    return header.fileversion;
  }
  
  int getNumberOfMatrices(){
    return header.nmatrices;
  }
 
 /*
  *Parses the matrix header, reading type, dimensions and element sizes
  */
  int parseMatrixHeader(ubyte[] buffer, int matrix, int start){
    MatrixHeader header = *cast(MatrixHeader*) inputbuffer[start..start+MatrixHeader.sizeof];
    writefln("      Matrix %d, type=%d, rows: %d, columns: %d", matrix, header.type, header.nrow, header.ncol);
    headers ~= header;
    return header.size;
  }
  
 /*
  *Parses the file header, versions and # of matrices
  */
  bool parseFileHeader(bool verbose){
    writeln("    File OK? " ~ to!string(checkBuffer(inputbuffer)));
    header = *cast(XgapBinHeader*) inputbuffer[0..XgapBinHeader.sizeof];
    writeln("    Version: " ~ to!string(header.fileversion));
    writeln("    Matrices: " ~ to!string(header.nmatrices));
    assert(correct,"Footprint failed");
    return correct;
  }
  
  this(in string filename,in bool verbose = false){
    assert(getSize(filename) < uint.max);
    inputbuffer = new ubyte[cast(uint)getSize(filename)];
    auto f = new File(filename,"rb");
    f.rawRead(inputbuffer);
    //Header
    parseFileHeader(verbose);

    //Loop through the matrices
    int skip = XgapBinHeader.sizeof;
    for(int m=0; m<getNumberOfMatrices(); m++){
      skips ~= skip;
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
  auto data = new XbinReader!RIL(outfn);
  assert(data.correct == true);
  assert(data.header.nmatrices == 3);
  assert(data.header.fileversion == [0,0,1, 'A']);
  //data.loadPhenotypes(data.matrices[0]);
  //data.loadGenotypes(data.matrices[1]);
  //data.loadMarkers(data.matrices[2]);
}


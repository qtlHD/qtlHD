/**
 * Write binary XGap format
 **/
 
module qtl.core.read_xgapbin;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.core.xgap;

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.file;

class XbinReader(XType){
  private File f;
  bool correct;
  int[] fileversion;
  int nmatrices;
  
  bool checkFootprint(in ubyte[] buffer){
    if(b_footprint == buffer) return true;
    return false;
  }
  
  bool checkBuffer(ubyte[] buffer){
    ubyte[] startprint = buffer[0..2];
    ubyte[] endprint = buffer[buffer.length-2..$];
    if(checkFootprint(startprint) && checkFootprint(endprint)){
      correct = true;
    }else{
      correct = false;
    }
    return correct;
  }
  
  int byteToInt(ubyte[] bits, bool little_endian = true ){
    int result = 0;
    if(little_endian){
      for(int n = bits.length-1; n >= 0; n--)
        result = (result << 8) +bits[ n ];
    }else{
      for(int n = 0; n < bits.length; n++)
        result = (result << 8) +bits[ n ];
    }
    return result;
  }
  
  double byteToDouble(ubyte[] bits){
    return cast(double)bits;
  }
  
  T[] toType(T)(ubyte[] buffer){
    T[] returnbuffer;
    foreach(int i, byte b ; buffer){
      returnbuffer ~= to!T(b);
    }
    return returnbuffer;
  }
  
  int[] getVersion(ubyte[] buffer){
    fileversion = toType!int(buffer[2..5]);
    return fileversion;
  }
  
  int getNumberOfMatrices(ubyte[] buffer){
    nmatrices = byteToInt(buffer[5..9]);
    return nmatrices;
  }
  
  int getMatrix(ubyte[] buffer, int matrix, int skip){
    int type = byteToInt(buffer[skip..skip+4]);
    int nrow = byteToInt(buffer[(skip+4)..(skip+8)]);
    int ncol = byteToInt(buffer[(skip+8)..(skip+12)]);
    int elem=0;
    int matrix_skip;
    writef("      Matrix %d, type=%d, rows: %d, columns: %d", matrix, type, nrow, ncol);
    if(type == MatrixType.VARCHARMATRIX){
      for(int x=0;x<(nrow*ncol);x++){
        elem += byteToInt(buffer[(skip+12+(x*4))..(skip+16+(x*4))]);
      }
      matrix_skip = elem+(4*((nrow*ncol)-1));
      writefln(", isize: %.2f, skip: %d", cast(double)elem/(nrow*ncol), matrix_skip);
    }else{
      elem = byteToInt(buffer[(skip+12)..(skip+16)]);
      matrix_skip = (nrow*ncol*elem);
      writefln(", isize: %d, skip: %d", elem, matrix_skip);
    }
    return 16+matrix_skip;
  }
  
  bool parseHeader(ubyte[] inputbuffer, bool verbose){
    writeln("    File OK? " ~ to!string(checkBuffer(inputbuffer)));
    writeln("    Version: " ~ to!string(getVersion(inputbuffer)));
    writeln("    Matrices: " ~ to!string(getNumberOfMatrices(inputbuffer)));
    assert(correct,"Footprint failed");
    assert(fileversion[0]==0 && fileversion[1]==0 && fileversion[2]==1,"Version incorrect");
    return correct;
  }
  
  this(in string filename,in bool verbose = false){
    assert(getSize(filename) < uint.max);
    ubyte[] inputbuffer = new ubyte[cast(uint)getSize(filename)];
    auto f = new File(filename,"rb");
    writeln("    Prereading");
    f.rawRead(inputbuffer);
    //Header
    writeln("    Read: " ~ to!string(getSize(filename)) ~ " bytes");
    parseHeader(inputbuffer,verbose);
    //Loop through the matrices
    int skip = 11;
    for(int m=0; m<getNumberOfMatrices(inputbuffer); m++){
      skip += getMatrix(inputbuffer, m, skip);
      assert(checkFootprint(inputbuffer[(skip)..(skip+2)]),"File corrupted ? No footprint at:" ~to!string(skip));
      skip += 2; //Skip the footprint
    }
  }
}

unittest{
  writeln("Unit test " ~ __FILE__);
  alias std.path.join join;
  auto infn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","test","data","input","multitrait.xbin");
  writeln("  - reading XBIN " ~ infn);
  auto data = new XbinReader!RIL(infn);
}


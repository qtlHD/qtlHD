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
  XgapFileHeader      header;         //Binary file header
  XgapMatrixHeader[]  headers;        //Matrix headers
  private File        f;              //File pointer
  ubyte[]             inputbuffer;    //Buffered file content
  bool                correct;        //Is the file correct?
  
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
  
  T[][] loadData(T)(XgapMatrixHeader h, int[] lengths, int start){
    T[][] data;
    int skip=start;
    for(int r=0;r<h.nrow;r++){
      T[] row;
      for(int c=0;c<h.ncol;c++){
        int s = lengths[(r*h.ncol)+c];
        row ~= convbyte!T(inputbuffer[skip..skip+s]);
        skip += s;
      }
      data ~= row;
    }
    return data;
  }
  
 /*
  * Loads the ith matrix of MatrixClasstype from the file
  */
  void load(MatrixClass c, int i = 0){
    int skip = XgapFileHeader.sizeof;
    foreach(XgapMatrixHeader h;headers){
      if(h.mclass==c){
        if(i==0){
        //extract matrix
        skip += XgapMatrixHeader.sizeof;
        int[] lengths;
        int start;
        if(h.type == MatrixType.VARCHARMATRIX){
          for(int x=0;x<(h.nrow*h.ncol);x++){
            lengths ~= byteToInt(inputbuffer[(skip+(x*int.sizeof))..(skip + int.sizeof + (x*int.sizeof))]);
          }
          start = (skip + int.sizeof + ((h.nrow*h.ncol)*int.sizeof));
        }else{
          for(int x=0;x<(h.nrow*h.ncol);x++){
            lengths ~= byteToInt(inputbuffer[skip..(skip+int.sizeof)]);
          }
          start = (skip + int.sizeof);
        }
        /* Here we need to convert data based on mclass and type */
        /* this should be solved by a single call, but i dun see how with all the templates going on */
        switch(h.mclass){
          case MatrixClass.EMPTY:
          break;
          case MatrixClass.PHENOTYPE:
          switch(h.type){
            case MatrixType.INTMATRIX:
              int[][]  phenotypes = loadData!int(h, lengths, start);
            break;        
            case MatrixType.DOUBLEMATRIX:
              double[][]  phenotypes = loadData!double(h, lengths, start);
            break;        
            case MatrixType.FIXEDCHARMATRIX:
              string[][]  phenotypes = loadData!string(h, lengths, start);
            break;        
            case MatrixType.VARCHARMATRIX:
              string[][]  phenotypes = loadData!string(h, lengths, start);
            break;
            default:
              throw new Exception("Unsupported format for PHENOTYPES");
            break;             
          }
          break;
          case MatrixClass.GENOTYPE:
          switch(h.type){
            case MatrixType.INTMATRIX:
              int[][] phenotypes = loadData!int(h, lengths, start);
            break;        
            case MatrixType.DOUBLEMATRIX:
              double[][]  phenotypes = loadData!double(h, lengths, start);
            break;        
            case MatrixType.FIXEDCHARMATRIX:
              string[][]  phenotypes = loadData!string(h, lengths, start);
            break;        
            default:
              throw new Exception("Unsupported format for GENOTYPES");
            break;   
          }
          break;
          case MatrixClass.MAP:
            switch(h.type){
              case MatrixType.VARCHARMATRIX:
                string[][]  map = loadData!string(h, lengths, start);
              break;       
              default:
                throw new Exception("Unsupported format for MAP");
              break;
            }
          break;
          case MatrixClass.ANNOTATION:

          break;
        }
        }else{
          i--;
        }
      }
      skip += h.size;
    }
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
  auto data = new XbinReader!RIL(outfn);
  assert(data.correct == true);
  assert(data.header.nmatrices == 3);
  assert(data.header.fileversion == [0,0,1, 'A']);
  data.load(MatrixClass.PHENOTYPE,0);
  data.load(MatrixClass.GENOTYPE,0);
  data.load(MatrixClass.MAP,0);
  //data.loadGenotypes(data.matrices[1]);
  //data.loadMarkers(data.matrices[2]);
}


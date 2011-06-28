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
  private File f;
  bool correct;       //Is the file correct?
  int[] fileversion;  //Version of the file
  int nmatrices;      //Number of matrices
  int[] skips;        //Stores the matrix skips, to get to the start of a matrix
  Matrix[] matrices;
  ubyte[] inputbuffer;
  
  bool checkFootprint(in ubyte[] buffer){
    if(xgap_footprint == buffer) return true;
    return false;
  }
  
  bool checkBuffer(ubyte[] buffer){
    ubyte[] startprint = buffer[0..8];
    ubyte[] endprint = buffer[buffer.length-8..$];
    if(checkFootprint(startprint) && checkFootprint(endprint)){
      correct = true;
    }else{
      correct = false;
    }
    return correct;
  }
  
  T[] toType(T)(ubyte[] buffer){
    T[] returnbuffer;
    foreach(int i, byte b ; buffer){
      returnbuffer ~= to!T(b);
    }
    return returnbuffer;
  }
  
  int[] getVersion(ubyte[] buffer){
    fileversion = toType!int(buffer[8..12]);
    return fileversion;
  }
  
  int getNumberOfMatrices(ubyte[] buffer){
    nmatrices = byteToInt(buffer[12..16]);
    return nmatrices;
  }
  
  void loadPhenotypes(Matrix m){
    for(int c=0;c<m.ncol;c++){
      Phenotype!double[] ps;
      for(int r=0;r<m.nrow;r++){
        string pheno = to!string(byteToDouble(inputbuffer[m.skip+(r*c)*m.lengths[0]..m.skip+m.lengths[0]+(r*c)*m.lengths[0]]));
        //writefln("%d %d %s",r,c,pheno);
        ps ~= set_phenotype!double(pheno);
      }
      phenotypes ~= ps;
    }
  }

  void loadGenotypes(Matrix m){
    for(int r=0;r<m.nrow;r++){
      Genotype!XType[] gs;
        for(int c=0;c<m.ncol;c++){  
        string geno = to!string(*cast(char*)inputbuffer[m.skip+(r*c)*m.lengths[0]..m.skip+m.lengths[0]+(r*c)*m.lengths[0]]);
        //writefln("%d %d %s",r,c,geno);
        gs ~= set_genotype!XType(geno);
      }
      genotypes ~= gs;
    }
  }
  
  void loadMarkers(Matrix m){
    int sofar=0;
    for(int r=0;r<m.nrow;r++){
      for(int c=0;c<m.ncol;c++){
        int s = m.lengths[(r*m.ncol)+c];
        string pheno = byteToString(inputbuffer[m.skip+sofar..m.skip+sofar+s]);
        //writefln("%d %d %s",r,c,pheno);
        sofar += s;
      }
    }
  }

  
 /*
  *Parses the matrix header, reading type, dimensions and element sizes
  */
  int parseMatrix(ubyte[] buffer, int matrix, int skip){
    int type = byteToInt(buffer[skip..skip+4]);
    int nrow = byteToInt(buffer[(skip+4)..(skip+8)]);
    int ncol = byteToInt(buffer[(skip+8)..(skip+12)]);
    
    int[] lengths;
    int elem=0;
    int matrix_skip;
    int length_skip=0;
    writef("      Matrix %d, type=%d, rows: %d, columns: %d", matrix, type, nrow, ncol);
    if(cast(MatrixType)type == MatrixType.VARCHARMATRIX){
      for(int x=0;x<(nrow*ncol);x++){
        int l = byteToInt(buffer[(skip+12+(x*4))..(skip+16+(x*4))]);
        lengths ~= l;
        elem += l;
      }
      matrix_skip = elem+(4*((nrow*ncol)-1));
      length_skip = (4*((nrow*ncol)-1));
      writefln(", isize: %.2f, skip: %d", cast(double)elem/(nrow*ncol), matrix_skip);
    }else{
      int l = byteToInt(buffer[(skip+12)..(skip+16)]);
      lengths ~= l;
      elem = l;
      matrix_skip = (nrow*ncol*elem);
      writefln(", isize: %d, skip: %d", elem, matrix_skip);
    }
    matrices ~= Matrix(cast(MatrixType)type,skip+16+length_skip,nrow,ncol,lengths);
    return 16+matrix_skip;
  }
  
 /*
  *Parses the file header, versions and # of matrices
  */
  bool parseFileHeader(bool verbose){
    writeln("    File OK? " ~ to!string(checkBuffer(inputbuffer)));
    writeln("    Version: " ~ to!string(getVersion(inputbuffer)));
    writeln("    Matrices: " ~ to!string(getNumberOfMatrices(inputbuffer)));
    assert(correct,"Footprint failed");
    assert(fileversion[0]==0 && fileversion[1]==0 && fileversion[2]==1,"Version incorrect");
    return correct;
  }
  
  this(in string filename,in bool verbose = false){
    assert(getSize(filename) < uint.max);
    inputbuffer = new ubyte[cast(uint)getSize(filename)];
    auto f = new File(filename,"rb");
    f.rawRead(inputbuffer);
    //Header
    writeln("    Read: " ~ to!string(getSize(filename)) ~ " bytes");
    parseFileHeader(verbose);

    //Loop through the matrices
    int skip = 24;
    for(int m=0; m<getNumberOfMatrices(inputbuffer); m++){
      skips ~= skip;
      skip += parseMatrix(inputbuffer, m, skip);
      assert(checkFootprint(inputbuffer[(skip)..(skip+8)]),"File corrupted ? No footprint at:" ~to!string(skip));
      skip += 8; //Skip the footprint
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
  assert(data.nmatrices == 3);
  assert(data.fileversion == [0,0,1]);
  data.loadPhenotypes(data.matrices[0]);
  data.loadGenotypes(data.matrices[1]);
  data.loadMarkers(data.matrices[2]);
}


/**
 * Write binary XGap files
 **/

module qtl.plugins.output.write_xgapbin;

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

import qtl.plugins.input.read_csvr;
import qtl.plugins.input.read_csv;

/** 
 * Convert a simple CSVR file containing marker names, chromosome nrs, position, 
 * phenotype and genotype - such as the multitrait.csvr file used in R/qtl.
 *
 * The file is parsed once on class instantiation. Elements can be queried.
 */
class BinaryWriter(Reader, XType) {
  private File f;
  Reader data;
 
  void myWrite(T)(T[] x,File f, MatrixType t = MatrixType.EMPTY, bool bin = true){
    if(bin){
      if(t == MatrixType.FIXEDCHARMATRIX || t == MatrixType.VARCHARMATRIX){
        foreach(T c;x){
          f.write(c);
        }
      }else{
        f.rawWrite(x);
      }
    }else{
      f.writeln(x);
    }
  }
  
  void write_matrix(T)(T[][] towrite, File outfile, MatrixType t = MatrixType.EMPTY){
    uint[1] type = [t];
    uint[] sizes =[towrite.length,towrite[0].length];
    switch(t){
      case MatrixType.EMPTY:
        sizes ~= cast(uint)0;
        break;
      case MatrixType.INTMATRIX:
        sizes ~= int.sizeof;
        break;        
      case MatrixType.DOUBLEMATRIX:
        sizes ~= double.sizeof;
        break;        
      case MatrixType.FIXEDCHARMATRIX:
        string s = to!string(towrite[0][0]);
        sizes ~= char.sizeof * s.length;
        break;        
      case MatrixType.VARCHARMATRIX:
        for(int r=0;r < towrite.length;r++){
          for(int c=0;c < towrite[r].length;c++){
            sizes ~= char.sizeof * to!string(towrite[r][c]).length;
          }
        }
        break;        
      default:
        break;
    }
    myWrite(type,outfile);
    myWrite(sizes,outfile);
    foreach(T[] e;towrite){
      myWrite(e,outfile,t);
    }
  }
  
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
  
  void write_binary(in string filename){
    f = File(filename,"wb");
    myWrite(b_footprint,f);
    myWrite(b_version,f);
    uint[1] nmatrix = [ 3 ];
    myWrite(nmatrix,f);
    write_matrix!(Phenotype!double)(data.phenotypes, f, MatrixType.DOUBLEMATRIX);
    myWrite(b_footprint,f);
    write_matrix!(Genotype!XType)(data.genotypes,f, MatrixType.FIXEDCHARMATRIX);
    myWrite(b_footprint,f);
    auto markermatrix = getMarkerInfoMatrix(data.markers);
    write_matrix!(string)(markermatrix,f, MatrixType.VARCHARMATRIX);
    myWrite(b_footprint,f);
    f.close();
  }

  this(Reader indata, in string outfilename){
    data = indata;
    write_binary(outfilename);
  }
}

unittest {
  writeln("Unit test " ~ __FILE__);
  alias std.path.join join;
  auto infn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","multitrait.csvr");
  auto outfn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","multitrait.xbin");
  writeln("  - reading CSVR " ~ infn ~" to " ~ outfn);
  auto data = new CsvrReader!RIL(infn);
  auto result = new BinaryWriter!(CsvrReader!RIL,RIL)(data,outfn);
  writefln("Size (txt to xbin): (%.2f Kb to %.2f Kb)", toKb(infn), toKb(outfn));
  
  auto infn1 = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","listeria.csv");
  auto outfn1 = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","listeria.xbin");
  writeln("  - reading CSVR " ~ infn1 ~" to " ~ outfn1);
  auto data1 = new ReadSimpleCSV!F2(infn1);
  auto result1 = new BinaryWriter!(ReadSimpleCSV!F2,F2)(data1,outfn1);
  writefln("Size (txt to xbin): (%.2f Kb to %.2f Kb)", toKb(infn1), toKb(outfn1));
}



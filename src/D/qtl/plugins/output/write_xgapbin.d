/**
 * Write binary XGap files
 **/

module qtl.plugins.output.write_xgapbin;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.core.xgap.xgap;

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
  
  void write_matrix(T)(T[][] towrite, File outfile, MatrixType t = MatrixType.EMPTY, MatrixClass mclass = MatrixClass.EMPTY){
    XgapMatrixHeader header = XgapMatrixHeader(xgap_footprint, t, mclass);
    header.nrow = towrite.length;
    header.ncol = towrite[0].length;
    int datasize = 0;
    int[] elementsizes;
    switch(t){
      case MatrixType.EMPTY:
        elementsizes ~= cast(uint)0;
        datasize = (0)*(header.nrow * header.ncol);
        break;
      case MatrixType.INTMATRIX:
        elementsizes ~= int.sizeof;
        datasize = (int.sizeof)*(header.nrow * header.ncol);
        break;        
      case MatrixType.DOUBLEMATRIX:
        elementsizes ~= double.sizeof;
        datasize = (double.sizeof)*(header.nrow * header.ncol);
        break;        
      case MatrixType.FIXEDCHARMATRIX:
        string s = to!string(towrite[0][0]);
        elementsizes ~= char.sizeof * s.length;
        datasize = (char.sizeof * s.length)*(header.nrow * header.ncol);
        break;        
      case MatrixType.VARCHARMATRIX:
        for(int r=0;r < towrite.length;r++){
          for(int c=0;c < towrite[r].length;c++){
            elementsizes ~= char.sizeof * to!string(towrite[r][c]).length;
            datasize += char.sizeof * to!string(towrite[r][c]).length;
          }
        }
        break;        
      default:
        break;
    }
    header.size = XgapMatrixHeader.sizeof + (elementsizes.length*4) + datasize;
    myWrite([header],outfile);
    myWrite(elementsizes,outfile);
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
    XgapFileHeader h = XgapFileHeader(xgap_footprint,xgap_version,[0,0,0,0,0],3,[0,0,0,0]);
    myWrite([h],f);
    write_matrix!(Phenotype!double)(data.phenotypes, f, MatrixType.DOUBLEMATRIX,MatrixClass.PHENOTYPE);
    write_matrix!(Genotype!XType)(data.genotypes,f, MatrixType.FIXEDCHARMATRIX,MatrixClass.GENOTYPE);
    write_matrix!(string)(getMarkerInfoMatrix(data.markers),f, MatrixType.VARCHARMATRIX,MatrixClass.MAP);
    myWrite(xgap_footprint,f);
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
  auto data = new CSVrReader!RIL(infn);
  auto result = new BinaryWriter!(CSVrReader!RIL,RIL)(data,outfn);
  writefln("Size (txt to xbin): (%.2f Kb to %.2f Kb)", toKb(infn), toKb(outfn));
  
  auto infn1 = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","listeria.csv");
  auto outfn1 = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","listeria.xbin");
  writeln("  - reading CSVR " ~ infn1 ~" to " ~ outfn1);
  auto data1 = new ReadSimpleCSV!F2(infn1);
  auto result1 = new BinaryWriter!(ReadSimpleCSV!F2,F2)(data1,outfn1);
  writefln("Size (txt to xbin): (%.2f Kb to %.2f Kb)", toKb(infn1), toKb(outfn1));
}



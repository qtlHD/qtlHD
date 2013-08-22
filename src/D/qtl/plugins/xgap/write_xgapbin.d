/**
 * Write binary XGap files
 **/

module qtl.plugins.xgap.write_xgapbin;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.core.xgap.xgap;
import example.genotype_examples;

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.file;

import qtl.plugins.csvr.read_csvr;
import qtl.plugins.csv.read_csv;

/** 
 * Convert a simple CSVR file containing marker names, chromosome nrs, position, 
 * phenotype and genotype - such as the multitrait.csvr file used in R/qtl.
 *
 * The file is parsed once on class instantiation. Elements can be queried.
 */
class BinaryWriter(Reader, XType) {
  private File f;
  Reader data;
 
  void myWrite(T)(T[] x,File f, MatrixType t = MatrixType.EMPTY, MatrixClass mclass = MatrixClass.EMPTY, bool bin = true){
    if(bin){
      if(t == MatrixType.FIXEDCHARMATRIX || t == MatrixType.VARCHARMATRIX){
        foreach(T c;x){
          f.write(c);
        }
      }else{
        if(mclass == MatrixClass.GENOTYPE){
          foreach(T c;x){
            f.rawWrite("1");
          }
        }else{
          f.rawWrite(x);
        }
      }
    }else{
      f.writeln(x);
    }
  }
  
  void write_matrix(T)(T[][] towrite, File outfile, MatrixType t = MatrixType.EMPTY, MatrixClass mclass = MatrixClass.EMPTY){
    XgapMatrixHeader header = XgapMatrixHeader(xgap_magicnumber, t, mclass);
    header.nrow = cast(int) towrite.length;
    header.ncol = cast(int) towrite[0].length;
    int datasize = 0;
    int[] elementsizes;
    switch(t){
      case MatrixType.EMPTY:
        elementsizes ~= cast(uint)0;
        datasize = (0)*(header.nrow * header.ncol);
        break;
      case MatrixType.INTMATRIX:
        elementsizes ~= int.sizeof;
        datasize = cast(int) (int.sizeof)*(header.nrow * header.ncol);
        break;        
      case MatrixType.DOUBLEMATRIX:
        elementsizes ~= double.sizeof;
        datasize = cast(int) (double.sizeof)*(header.nrow * header.ncol);
        break;        
      case MatrixType.FIXEDCHARMATRIX:
        string s = to!string(towrite[0][0]);
        elementsizes ~= cast(int) (char.sizeof * s.length);
        datasize = cast(int) (char.sizeof * s.length)*(header.nrow * header.ncol);
        break;        
      case MatrixType.VARCHARMATRIX:
        for(int r=0;r < towrite.length;r++){
          for(int c=0;c < towrite[r].length;c++){
            elementsizes ~= cast(int) (char.sizeof * to!string(towrite[r][c]).length);
            datasize += char.sizeof * to!string(towrite[r][c]).length;
          }
        }
        break;        
      default:
        break;
    }
    header.size = cast(int) (XgapMatrixHeader.sizeof + (elementsizes.length*4) + datasize);
    myWrite([header],outfile);
    //Here we need to write rownames and colnames, Now STUBS
    int[] rowlength;
    string[] rownames;
    for(auto r=0;r<header.nrow;r++){
      rowlength ~= 1;
      rownames ~= "A";
    }
    myWrite(rowlength,outfile);
    myWrite(rownames,outfile,MatrixType.VARCHARMATRIX);
    
    int[] collength;
    string[] colnames;
    for(auto c=0;c<header.ncol;c++){
      collength ~= 2;
      colnames ~= "BA";
    }

    myWrite(collength,outfile);
    myWrite(colnames,outfile,MatrixType.VARCHARMATRIX);
  
    myWrite(elementsizes,outfile);
    foreach(T[] e;towrite){
//      writeln(e);
      myWrite(e,outfile,t);
    }
  }
  
  void write_binary(in string filename){
    f = File(filename,"wb");
    XgapFileHeader header = XgapFileHeader(xgap_magicnumber,xgap_version,3);
    myWrite([header],f);
    write_matrix!(Phenotype)(data.phenotypes, f, MatrixType.DOUBLEMATRIX,MatrixClass.PHENOTYPE);
    write_matrix!(GenotypeSymbolMapper)(data.genotypecombinator,f, MatrixType.VARCHARMATRIX,MatrixClass.GENOTYPE);
    write_matrix!(string)(getMarkerInfoMatrix(data.markers),f, MatrixType.VARCHARMATRIX,MatrixClass.MAP);
    myWrite(xgap_magicnumber,f);
    f.close();
  }

  this(Reader indata, in string outfilename){
    data = indata;
    write_binary(outfilename);
  }
}

unittest {
  writeln("Unit test " ~ __FILE__);
  alias std.path.buildPath buildPath;
  auto infn = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","test","data","input","multitrait.csvr"));
  auto outfn = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","test","data","input","multitrait.xbin"));
  writeln("  - reading CSVR " ~ infn ~" to " ~ outfn);
  auto data = new CSVrReader!(RISELF,ObservedRISELF)(infn);
  auto result = new BinaryWriter!(CSVrReader!(RISELF,ObservedRISELF),RISELF)(data,outfn);
  writefln("Size (txt to xbin): (%.2f Kb to %.2f Kb)", toKb(infn), toKb(outfn));

  auto infn1 = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","test","data","input","listeria.csv"));
  auto outfn1 = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","test","data","input","listeria.xbin"));
  writeln("  - reading CSVR " ~ infn1 ~" to " ~ outfn1);
  auto data1 = new ReadSimpleCSV!(F2,ObservedF2)(infn1);
  auto result1 = new BinaryWriter!(ReadSimpleCSV!(F2,ObservedF2),F2)(data1,outfn1);
  writefln("Size (txt to xbin): (%.2f Kb to %.2f Kb)", toKb(infn1), toKb(outfn1));
}



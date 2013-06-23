/**
 * XGap converter - translate between XGap types and qtlHD primitives
 **/

module qtl.core.xgap.converter;

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.file;

import qtl.core.primitives;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.core.xgap.xgap;
import qtl.core.xgap.read_xgapbin;

class XbinConverter{
  
  this(){
  
  }

  GenotypeCombinator[][] toGenotypes(T)(XgapMatrix m){
    GenotypeCombinator[][] all_genotypes;
    for(int r=0;r<m.header.nrow;r++){
      GenotypeCombinator[] genotype;
      for(int c=0;c<m.header.ncol;c++){
        //genotype ~= T.decode((cast(StringMatrix)(m.data)).data[r][c]);
      }
      all_genotypes ~= genotype;
    }
    return all_genotypes;
  }
  
  Phenotype[][] toPhenotype(T)(XgapMatrix m){
    Phenotype[][] all_phenotypes;
    for(int r=0;r<m.header.nrow;r++){
      Phenotype[] phenotype;
      for(int c=0;c<m.header.ncol;c++){
        phenotype ~= set_phenotype(to!string((cast(DoubleMatrix)(m.data)).data[r][c]));
      }
      all_phenotypes ~= phenotype;
    }
    return all_phenotypes;
  }
  
  Marker[] asMarkers(XgapMatrix m){
    if(m.header.type != MatrixType.VARCHARMATRIX){
      throw new Exception("Marker matrices should be of type VARCHAR");
    }else{
      Marker[] markers;
      for(int r=0;r<m.header.nrow;r++){
        string mname = (cast(StringMatrix)m.data).data[r][0];
        string mchr = (cast(StringMatrix)m.data).data[r][1];
        double mloc = to!double((cast(StringMatrix)m.data).data[r][2]);
        markers ~= new Marker(mloc,mname);
      }
      return markers;
    }
  }
}

unittest{
  writeln("Unit test " ~ __FILE__);
  alias std.path.buildPath buildPath;
  auto outfn = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","test","data","input","multitrait.xbin"));
  auto data = new XbinReader(outfn);
  auto convertor = new XbinConverter();
  auto phenotypes = convertor.toPhenotype(data.load(0));
  auto genotypes = convertor.toGenotypes!ObservedRISELF(data.load(1));
  auto markers = convertor.asMarkers(data.load(2));
}

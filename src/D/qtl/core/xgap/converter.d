/**
 * XGap converter
 **/

module qtl.core.xgap.converter;

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.file;

import qtl.core.primitives;
import qtl.core.phenotype;
import qtl.core.deprecate.genotype_enum;
import qtl.core.xgap.xgap;
import qtl.core.xgap.read_xgapbin;

class XbinConverter{
  
  this(){
  
  }
  
  Genotype!T[][] toGenotypes(T)(XgapMatrix m){
    Genotype!T[][] all_genotypes;
    for(int r=0;r<m.header.nrow;r++){
      Genotype!T[] genotype;
      for(int c=0;c<m.header.ncol;c++){
        genotype ~= set_genotype!T((cast(StringMatrix)(m.data)).data[r][c]);
      }
      all_genotypes ~= genotype;
    }
    return all_genotypes;
  }
  
  Phenotype!T[][] toPhenotype(T)(XgapMatrix m){
    Phenotype!T[][] all_phenotypes;
    for(int r=0;r<m.header.nrow;r++){
      Phenotype!T[] phenotype;
      for(int c=0;c<m.header.ncol;c++){
        phenotype ~= set_phenotype!T(to!string((cast(DoubleMatrix)(m.data)).data[r][c]));
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
  auto outfn = to!string(dirName(__FILE__) ~ sep ~ buildPath("..","..","..","..","..","test","data","input","multitrait.xbin"));
  auto data = new XbinReader(outfn);
  auto convertor = new XbinConverter();
  auto phenotypes = convertor.toPhenotype!double(data.load(0));
  auto genotypes = convertor.toGenotypes!RIL(data.load(1));
  auto markers = convertor.asMarkers(data.load(2));
}

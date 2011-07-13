/**
 * XGap converter
 **/

module qtl.core.xgap.converter;

import std.file;
import std.conv;
import std.stdio;

import qtl.core.primitives;

class XbinConverter(XType){
  
  this(){
  
  }
  
  Genotype!T[][] toGenotypes(T)(XgapMatrix m){
    Genotype!T[][] all_genotypes;
    for(int r=0;r<m.nrows;r++){
      Genotype!T[] genotype;
      for(int c=0;c<m.ncols;c++){
        genotype ~= set_genotype!T((cast(StringMatrix)(m.data)).data[r][c]);
      }
      all_genotypes ~= genotype;
    }
    return all_genotypes;
  }
  
  Phenotype!T[][] toPhenotype(T)(XgapMatrix m){
    Phenotype!T[][] all_phenotypes;
    for(int r=0;r<m.nrows;r++){
      Phenotype!T[] phenotype;
      for(int c=0;c<m.ncols;c++){
        phenotype ~= set_phenotype!T(to!string((cast(DoubleMatrix)(m.data)).data[r][c]));
      }
      all_phenotypes ~= phenotype;
    }
    return all_phenotypes;
  }
  
  Marker[] toMarkers(T)(XgapMatrix m){
    if(m.header.type != MatrixType.VARCHARMATRIX){
      throw new Exception("Marker matrices should be of type VARCHAR");
    }else{
    
    }
  }
}

unittest{
  writeln("Unit test " ~ __FILE__);
}
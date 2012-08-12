/**
 * Phenotype module
 */

module qtl.core.phenotype;

import qtl.core.primitives;
import std.conv;
import std.stdio;
import std.string;

Phenotype!T set_phenotype(T)(in string s) {
  // writeln(s);
  Phenotype!T p;
  if(s == "NA" || s == "-"){
    p.value = to!T(PHENOTYPE_NA);
  }else{
    if(s.indexOf(".") != -1){  // FIXME: this should only be for floats
      p.value = to!T(s);
    }else{
      p.value = to!T(s~".0");
    }
  }
  return p;
}

/**
 * Check whether a phenotype is missing
 */
bool isNA(T)(Phenotype!T phe) { 
  return(phe.value == PHENOTYPE_NA);
}

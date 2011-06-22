/**
 * Phenotype module
 */

module qtl.core.phenotype;

import qtl.core.primitives;
import std.conv;
import std.stdio;

immutable PHENOTYPE_NA = double.max;

Phenotype!T set_phenotype(T)(in string s) {
  // writeln(s);
  Phenotype!T p;
  if(s == "NA" || s == "-") 
    p.value = PHENOTYPE_NA;
  else
    p.value = to!T(s);
  return p;
}


module qtl.core.phenotype;

import qtl.core.primitives;
import std.conv;
import std.stdio;

immutable PHENO_NA = double.max;

Phenotype!T set_phenotype(T)(in string s) {
  writeln(s);
  Phenotype!T p;
  if (s == "NA") 
    p.value = PHENO_NA;
  else
    p.value = to!T(s);
  return p;
}

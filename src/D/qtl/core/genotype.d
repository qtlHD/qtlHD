/**
 * Genotype module
 */

module qtl.core.genotype;

import qtl.core.primitives;
import std.conv;
import std.stdio;

immutable GENOTYPE_NA = -1;

enum RIL { NA, A, B, C, H };
enum BC  { NA, A, B, C, H };

Genotype!T set_genotype(T)(in string s) {
  // writeln(s);
  Genotype!T g;
  switch(s) {
    case "-":
      g.value = T.NA;
      break;  
    case "A":
      g.value = T.A;
      break;
    case "B":
      g.value = T.B;
      break;
    case "C":
      g.value = T.C;
      break;
    case "H":
      g.value = T.H;
      break;
  }
  return g;
}

/**
 * Genotype module
 */

module qtl.core.genotype;

import qtl.core.primitives;
import std.conv;
import std.stdio;

immutable GENOTYPE_NA = -1;

enum RIL { NA, A, B };
enum F2  { NA, A, H, B, HorB, HorA }; 
enum BC  { NA, A, H };

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
    case "H":
      g.value = T.H;
      break;
    case "B":
      g.value = T.B;
      break;
    case "C":
      g.value = T.HorB;
      break;
    case "D":
      g.value = T.HorA;
      break;
  }
  return g;
}


unittest {
  writeln("Unit test " ~ __FILE__);
  Genotype!F2[] gs;
  gs ~= set_genotype!F2("A");
  // assert(gs[0] == F2.A);
}

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
  gs ~= set_genotype!F2("-");
  gs ~= set_genotype!F2("A");
  gs ~= set_genotype!F2("B");
  gs ~= set_genotype!F2("H");
  gs ~= set_genotype!F2("C");
  gs ~= set_genotype!F2("D");
  assert(gs[0].value == F2.NA);
  assert(gs[1].value == F2.A);
  assert(gs[2].value == F2.B);
  assert(gs[3].value == F2.H);
  assert(gs[4].value == F2.HorB);
  assert(gs[5].value == F2.HorA);
}

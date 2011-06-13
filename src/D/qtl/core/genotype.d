/**
 * Genotype module
 */

module qtl.core.genotype;

import qtl.core.primitives;
import std.conv;
import std.stdio;

immutable GENOTYPE_NA = -1;

enum RIL { NA, A, B };
enum F2  { NA, A, H, B, C, D }; // C = not A = H or B; D = not B = A or H
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
      g.value = T.C;
      break;
    case "D":
      g.value = T.D;
      break;
  }
  return g;
}

/**
 * Genotype module
 */

module qtl.core.genotype;

import qtl.core.primitives;
import std.conv;
import std.stdio;

enum RIL { NA = GENOTYPE_NA, A, B };
enum F2  { NA = GENOTYPE_NA, A, H, B, HorB, HorA }; 
enum F2pk {AA, AB, BA, BB}; // pk = "phase known"
enum BC  { NA = GENOTYPE_NA, A, H };

Genotype!T set_genotype(T)(in string s) {
  Genotype!T g;
  static if (T.stringof == "RIL") {
    switch(s) {
      case "-":
        g.value = T.NA;
        break;  
      case "A":
        g.value = T.A;
        break;
      case "AA":
        g.value = T.A;
        break;        
      case "B":
        g.value = T.B;
        break;
      case "BB":
        g.value = T.A;
        break;
      default:
        throw new Exception("Unknown genotype " ~ s ~ " for " ~ T.stringof);
    }
  }
  static if (T.stringof == "BC") {
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
      default:
        throw new Exception("Unknown genotype " ~ s ~ " for " ~ T.stringof);
    }
  }
  static if (T.stringof == "F2") {
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
      default:
        throw new Exception("Unknown genotype " ~ s ~ " for " ~ T.stringof);
    }
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

unittest {
  Genotype!RIL[] ril;
  ril ~= set_genotype!RIL("-");
  ril ~= set_genotype!RIL("A");
  ril ~= set_genotype!RIL("B");
  // ril ~= set_genotype!RIL("H"); raises exception
  // ril ~= set_genotype!RIL("C"); raises exception
  assert(ril[0].value == RIL.NA);
  assert(ril[1].value == RIL.A);
  assert(ril[2].value == RIL.B);
}

import std.algorithm;

unittest {
  Genotype!BC[] bc;
  bc = [ set_genotype!BC("-"), set_genotype!BC("A"), set_genotype!BC("H") ];
  auto list = map!((x) { return x.value; })(bc);
  assert(list[0] == BC.NA);
  assert(list[1] == BC.A);
  assert(list[2] == BC.H);
}

unittest {
  // This should compile
  auto map = new FullMap!F2();
  /*
  foreach ( c ; map.chromosome_map ) {
    auto markers = c.markers;
    foreach ( m ; markers.list ) {
    }
  }
  */
}

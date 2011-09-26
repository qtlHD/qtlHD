/**
 * Genotype module
 */

module qtl.core.genotype;

import std.conv;
import std.stdio;
import qtl.core.primitives;

/**

  Discussing storing of genotypes:

  The number of possible true, or real, genotypes at a marker location is
  limited, as they depend on the number of founders (K^2 types). Unfortunately,
  due to the scoring technology, there are a lot more 'observed' genotypes -
  i.e.  combinations of possible true genotypes,

  Inside each dataset, however, there is a limited number of stored
  combinations. So, rather than setting all types in advance we only create the
  actual types used in the dataset. The genotype matrix contains pointers to a
  bit array. The bit array references combinations of supported types.
  Supported types are Tuples of 'alleles', where alleles are numbers referring
  to the founders. I.e.

  founders:       [1]   [2]   [3]   [4]
  
  actual types:   [1,1] [1,2] [2,2] [2,3] [2,4] etc.  (TrueGenotype)

  combinations: (GenotypeCombinator[])

   -1        ->   NA
    0        ->   [ 1      0     0     0    0 ]   i.e. [1,1]
    1        ->   [ 0      1     1     0    0 ]   i.e. [1,2] or [2,2] 
    2        ->   [ 0      0     1     0    0 ]   i.e. [2,2]

  (GenotypeCombinatorIndex - the left column)

  Types are defined in primitives.d.

  Every marker column has its own combinations defined. So combination #1 may
  refer to different genotypes, between different marker columns. 

  Above bit matrix may suggest a form of binary storage/packing/unpacking. In
  reality we store it as a list of references to TrueGenotype. This is because
  the bit vectors would get really wide with many true genotypes (K^2 sized).
  So in above example:

   Index          GenotypeCombinator
   -1        ->   NA
    0        ->   [ 0 ]    i.e. [1,1]
    1        ->   [ 1,2 ]  i.e. [1,2] or [2,2] 
    2        ->   [ 2 ]    i.e. [2,2]

  Note that types can be reused - i.e. there are no duplicates defined.

  Each marker has a list of GenotypeCombinators. The actual 
  GenotypeCombinators are also maintained in a separate list - so as
  to share observed types, rather than duplicating them for each 
  marker. For this we define ObservedGenotypes, which maintains
  this list.

  Sex should be queried at the chromosome level.

  Note that outbreeding can be supported by creating artificial founders
  consisting of singular genotypes (similar to a RIL), with one maternal and
  one paternal genotype as founders.

  With SNP markers, two representations can be supported using this genotyping
  system. In the first, each SNP is treated as a marker - therefore there are
  always two genotypes involved. This is the simplest approach. The second
  approach works at the allele level. SNPs are combined into ORFs, or genes.
  Multiple allele combinations can be supported with this genotyping
  representation.

 */

/**
 * FounderIndex is an index into Founders. This may be turned into a Founder
 * ref, later.
 */

alias uint FounderIndex;   

/**
 * A true genotype consists of a Tuple of genotypes/alleles - one from each
 * (founder) parent. These genotypes can, potentially, be directional 
 */

class TrueGenotype {
  Tuple!(FounderIndex,FounderIndex) founders;
  this(FounderIndex founder1,FounderIndex founder2) {
    founders = Tuple!(FounderIndex,FounderIndex)(founder1, founder2);
  }
  this(TrueGenotype g) { founders = g.founders; }
  auto homozygous()   { return founders[0] == founders[1]; };
  auto heterozygous() { return !homozygous(); };
  int opCmp(Object other) {
    auto founders2 = (cast(TrueGenotype)other).founders;
    if (founders[0] > founders2[0]) return 1;
    if (founders[0] < founders2[0]) return -1;
    if (founders[1] > founders2[1]) return 1;
    if (founders[1] < founders2[1]) return -1;
    return 0;
  }
  bool opEquals(Object other) {
    auto founders2 = (cast(TrueGenotype)other).founders;
    return founders[0] == founders2[0] && founders[1] == founders2[1];
  }
  string toString() { 
    uint f0 = founders[0]; uint f1 = founders[1];
    return '(' ~ to!string(f0) ~ ',' ~ to!string(f1) ~ ')';
  }
}

/**
 * GenotypeCombinator points to TrueGenoTypes - we have a list 
 * of these for each observed type - as it is defined in the
 * dataset.
 */

class GenotypeCombinator {
  TrueGenotype[] list;
  string name;

  this(string name = "?") { this.name = name; }

  // uniquely add a genotype. Return match.
  TrueGenotype add(TrueGenotype g) { 
    foreach(m; list) { if (m.founders == g.founders) return m; }
    list ~= g; 
    return g;
  }
  GenotypeCombinator opOpAssign(string op)(TrueGenotype g) if (op == "~") {
    add(g);
    return this;
  }
  auto length() { return list.length; }
  bool opEquals(Object other) {
    auto rhs = cast(GenotypeCombinator)other;
    return (list.sort == rhs.list.sort); // probably not the fastest way ;)
  }
  string toString() {
    if (isNA) return "[NA]";
    return to!string(list);
  }
  bool isNA() { return length == 0; }
}

/** 
 * ObservedGenotypes tracks all the observed genotypes in a dataset, for the
 * full set, or at a marker position.  This is a convenience class, mostly.
 */

class ObservedGenotypes {
  GenotypeCombinator[] list;
  /// Pass in a combination of genotypes, add if not already in list and
  /// return the actual match.
  GenotypeCombinator add(GenotypeCombinator c) { 
    foreach(m; list) { if (m == c) return m; }
    list ~= c; 
    return c;
  }
  GenotypeCombinator find(GenotypeCombinator combinator) {
    return null;
  }
  ObservedGenotypes opOpAssign(string op)(GenotypeCombinator c) if (op == "~") {
    add(c);
    return this;
  }
  auto length() { return list.length; }
}

/**
 * New style genotypes - support for observed genotypes
 */

import std.typecons; 

/**
 * Unit tests for TrueGenotype and GenotypeCombinators
 */

unittest {
  writeln("Unit test" ~ __FILE__);
  // at this point founders are simply numbers
  FounderIndex[] founder = [ 1, 2, 3, 4, 5 ];
  // create a few genotypes (at a marker location)
  auto g1  = new TrueGenotype(founder[0],founder[2]);
  auto g2  = new TrueGenotype(founder[1],founder[1]); // could be a RIL
  auto g2a = new TrueGenotype(founder[1],founder[1]); // duplicate
  auto g3  = new TrueGenotype(founder[2],founder[3]);
  auto g4  = new TrueGenotype(founder[4],founder[3]);
  assert(g1.heterozygous());
  assert(g2.homozygous());
  assert(g2.founders[0] == 2);
  // observed genotypes can be any combination of true genotypes
  auto observed_combi1 = new GenotypeCombinator;
  observed_combi1 ~= g1;
  observed_combi1 ~= g2;
  observed_combi1 ~= g2; // tests also for duplicate add
  writeln(observed_combi1.list);
  assert(!observed_combi1.isNA);
  assert(observed_combi1.list.length == 2);
  observed_combi1 ~= g2; // tests also for duplicate add
  assert(observed_combi1.list.length == 2);
  auto testg = observed_combi1.add(g2a); // we want the return type
  assert(testg == g2);
  assert(observed_combi1.list.length == 2);
  auto observed_combi2 = new GenotypeCombinator;
  observed_combi2.list ~= g2;
  assert(observed_combi2.list.length == 1);
  assert(observed_combi2 == observed_combi2);
  assert(observed_combi1 != observed_combi2);
  auto observed_combi1a = new GenotypeCombinator;
  observed_combi1a ~= g2a;
  observed_combi1a ~= g1;
  assert(observed_combi1 == observed_combi1a);
  // observed_combi already acts as an index, so now we only need 
  // to store the observed genotypes with the marker. The marker/genotype
  // matrix stores references to observed_combiX by reference. Say
  // we have a marker column, and 2 individuals (2 observed genotypes):
  auto observed = new GenotypeCombinator[](2);
  observed[0] = new GenotypeCombinator;
  observed[1] = new GenotypeCombinator;
  // currently there is no content:
  assert(to!string(observed[0]) == "[NA]", to!string(observed[0]));
  observed[1] ~= g1;
  observed[1] ~= g2; // add 2nd possible genotype
  observed[1] ~= g2; // duplicate type 
  assert(to!string(observed[1]) == "[(1,3), (2,2)]");
  assert(observed[1].list == [g1,g2]);
  assert(observed[1] != observed[0]);
  // you see, each observed/individual simply contains a list of true (possible)
  // genotypes. To keep track of observed genotypes you may use
  auto tracker = new ObservedGenotypes();
  tracker ~= observed[0];
  tracker ~= observed[1];
  auto test = tracker.add(observed[1]); // add duplicate
  assert(test == observed[1]);
  assert(test == observed_combi1); // No duplication!
  assert(tracker.length == 2);
}

/**
 * RIL set using Genotype combinator
 * RIL { NA, A, B };
 */

unittest {
  auto na = new GenotypeCombinator("NA","-");
  auto a  = new GenotypeCombinator("A");
  a ~= new TrueGenotype(0,0);
  auto b  = new GenotypeCombinator("B");
  b ~= new TrueGenotype(1,1);
  assert(na.name == "NA");
  assert(a.name == "A");
  auto tracker = new ObservedGenotypes();
  tracker ~= na;
  tracker ~= a;
  tracker ~= b;
  auto ril = GenotypeCombinator[];  // one set of observed genotypes
  ril ~= set_genotype(tracker,"-");
  ril ~= set_genotype(tracker,"NA");
  ril ~= set_genotype(tracker,"A");
  ril ~= set_genotype(tracker,"B");
  // ril ~= set_genotype(tracker,"C"); raises exception
  assert(ril[0].value.isNA);
  assert(ril[1].value.name == "A");
  assert(ril[2].value.name == "B");
}

/** 
 * Emulate BC using Genotype combinator
 */

/**
 * Enum (old) style genotypes
 *
 * The follwing are the enum typed genotypes - which are not as
 * flexible as the newer genotyping scheme, explained above.
 */

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
  assert(to!string(ril[0]) == "-");
  assert(to!string(ril[1]) == "A", to!string(ril[1]));
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

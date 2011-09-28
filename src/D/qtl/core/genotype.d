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

   -1        ->   []                              i.e. NA
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
   -1        ->   []       i.e. NA
    0        ->   [ 0 ]    i.e. [1,1]
    1        ->   [ 1,2 ]  i.e. [1,2] or [2,2] 
    2        ->   [ 2 ]    i.e. [2,2]

  Note that types are reused - i.e. there are no duplicates defined.

  Each marker has a list of GenotypeCombinators. The actual GenotypeCombinators
  are also maintained in a separate list - so as to share observed types, a
  tracker, rather than duplicating them for each marker. For this we define
  ObservedGenotypes, which maintains this list. The tracker can be used
  for a full dataset, or for each marker individually.

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
 * (founder) parent. These genotypes can be directional 
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
 * dataset. The list is filled with a combination of true
 * genotypes, while reading in the data(file). This means the 
 * minimal amount of memory is used for the (exhaustive) set of
 * combinations of genotypes.
 *
 * To prevent duplication of stored types, when a new combination
 * is added, the add function only stores new types, and returns
 * the existing one in the list.
 *
 * GenotypeCombinator may be parametrized in the future - i.e.
 * to support other underlying true genotypes, and the accepted
 * encoding.
 */

class GenotypeCombinator {
  alias string Encoding;
  GenotypeCombinator value; 

  TrueGenotype[] list;
  string name; // the display name
  Encoding encoding[]; // a list of encoding types

  // initialize a combinator by name. If code is undefined name 
  // is also used for input encoding
  this(string name, Encoding code = null) { 
    value = this;
    this.name = name; 
    if (code) add_encoding(code);
    else add_encoding(name);
  }
  void add_encoding(Encoding code) { this.encoding ~= code; }
  // see if an input matches
  bool match(in Encoding code) {
    foreach (n; encoding) {
      if (code == n) return true; 
    }
    return false;
  }
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
  GenotypeCombinator opOpAssign(string op)(GenotypeCombinator c) if (op == "~") {
    foreach(g ; c.list ) { add(g); };
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
 * full set, or at a marker position.  This tracker is a convenience class,
 * mostly.
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
  /// Decode an input to an (observed) genotype
  GenotypeCombinator decode(in string s) {
    foreach(m; list) { if (m.match(s)) return m; }
    throw new Exception("Unknown genotype " ~ s);
  }
  /// Define =~ to combine combinators
  ObservedGenotypes opOpAssign(string op)(GenotypeCombinator c) if (op == "~") {
    add(c);
    return this;
  }
  auto length() { return list.length; }
  string toString() {
    return to!string(list);
  }
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
  auto observed1 = new GenotypeCombinator("OBSERVE1");
  observed1 ~= g1;
  observed1 ~= g2;
  observed1 ~= g2; // tests also for duplicate add
  writeln(observed1.list);
  assert(!observed1.isNA);
  assert(observed1.list.length == 2);
  observed1 ~= g2; // tests also for duplicate add
  assert(observed1.list.length == 2);
  auto testg = observed1.add(g2a); // we want the return type
  assert(testg == g2);
  assert(observed1.list.length == 2);
  auto observed2 = new GenotypeCombinator("OBSERVE2");
  observed2.list ~= g2;
  assert(observed2.list.length == 1);
  assert(observed2 == observed2);
  assert(observed1 != observed2);
  auto observed1a = new GenotypeCombinator("OBSERVE1A");
  observed1a ~= g2a;
  observed1a ~= g1;
  assert(observed1 == observed1a);
  // observed already acts as an index, so now we only need 
  // to store the observed genotypes with the marker. The marker/genotype
  // matrix stores references to observedX by reference. Say
  // we have a marker column, and 2 individuals (2 observed genotypes):
  auto observed = new GenotypeCombinator[](2);
  observed[0] = new GenotypeCombinator("0");
  observed[1] = new GenotypeCombinator("1");
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
  assert(test == observed1); // No duplication!
  assert(tracker.length == 2);
}

/**
 * RIL set using Genotype combinator
 * RIL { NA, A, B };
 */

class RIL {
  GenotypeCombinator NA, A, B;
  this() {
    NA = new GenotypeCombinator("NA");
    A  = new GenotypeCombinator("A");
    A ~= new TrueGenotype(0,0);
    B  = new GenotypeCombinator("B");
    B ~= new TrueGenotype(1,1);
    NA.add_encoding("-"); 
    A.add_encoding("AA"); 
    B.add_encoding("BB"); 
  }
}

unittest {
  auto ril = new RIL;
  auto tracker = new ObservedGenotypes();
  tracker ~= ril.NA;
  tracker ~= ril.A;
  tracker ~= ril.B;
  assert(tracker.decode("-") == ril.NA);
  assert(tracker.decode("A") == ril.A);
  assert(tracker.decode("AA") == ril.A);
}

unittest {
  // We can also roll our own types (to create a RIL)
  auto NA = new GenotypeCombinator("NA");
  auto A  = new GenotypeCombinator("A");
  A ~= new TrueGenotype(0,0);
  auto B  = new GenotypeCombinator("B");
  B ~= new TrueGenotype(1,1);
  assert(NA.name == "NA");
  assert(A.name == "A");
  auto tracker = new ObservedGenotypes();
  tracker ~= NA;
  tracker ~= A;
  tracker ~= B;
  GenotypeCombinator ril[];  // create a set of observed genotypes
  // now find them by name
  NA.add_encoding("-"); // also support dash inputs for NA
  ril ~= tracker.decode("-");
  ril ~= tracker.decode("NA");
  ril ~= tracker.decode("A");
  ril ~= tracker.decode("B");
  // ril ~= tracker.decode("C"); // raises exception!
  assert(ril[0].isNA);
  assert(ril[1].isNA);
  assert(ril[2].name == "A");
  assert(to!string(ril[2]) == "[(0,0)]");
  assert(ril[3].name == "B");
  assert(to!string(ril[3]) == "[(1,1)]");
}

class ObservedRIL {
}

/**
 * BC implementation
 * BC  { NA, A, H };
 */

class BC {
  GenotypeCombinator NA, A, H;
  this() {
    NA = new GenotypeCombinator("NA");
    A  = new GenotypeCombinator("A");
    A ~= new TrueGenotype(0,0);
    H  = new GenotypeCombinator("H");
    H ~= new TrueGenotype(1,0);
    NA.add_encoding("-"); 
  }
}

class ObservedBC {
  BC crosstype;
  ObservedGenotypes tracker;
  this() {
    auto bc = new BC;
    tracker = new ObservedGenotypes();
    tracker ~= bc.NA;
    tracker ~= bc.A;
    tracker ~= bc.H;
    crosstype = bc;
  }
  /// Decode an input to an (observed) genotype
  auto decode(in string s) {
    return tracker.decode(s);
  }
  auto length() { return tracker.length; }
  string toString() {
    return to!string(tracker);
  }

}

unittest {
  auto bc = new BC;
  auto tracker = new ObservedGenotypes();
  tracker ~= bc.NA;
  tracker ~= bc.A;
  tracker ~= bc.H;
  assert(tracker.decode("-") == bc.NA);
  assert(tracker.decode("A") == bc.A);
  assert(tracker.decode("H") == bc.H);
}

/**
 * F2 set using Genotype combinator
 * F2  { NA, A, H, B, HorB, HorA }; 
 */

class F2 {
  GenotypeCombinator NA, A, B, H, HorB, HorA;
  alias HorB C;
  alias HorA D;
  this() {
    auto aa = new TrueGenotype(0,0);
    auto bb = new TrueGenotype(1,1);
    auto ab = new TrueGenotype(0,1);
    NA = new GenotypeCombinator("NA");
    A  = new GenotypeCombinator("A");
    A ~= aa;
    B  = new GenotypeCombinator("B");
    B ~= bb;
    H  = new GenotypeCombinator("H");
    H ~= ab;
    HorB  = new GenotypeCombinator("HorB","C");
    HorB ~= ab;
    HorB ~= bb;
    HorA  = new GenotypeCombinator("HorA","D");
    HorA ~= ab;
    HorA ~= aa;
    NA.add_encoding("-"); 
    A.add_encoding("AA"); 
    B.add_encoding("BB"); 
    H.add_encoding("AB"); 
    H.add_encoding("BA"); 
  }
}

class ObservedF2 {
  F2 crosstype;
  ObservedGenotypes tracker;
  this() {
    auto f2 = new F2;
    tracker = new ObservedGenotypes();
    tracker ~= f2.NA;
    tracker ~= f2.A;
    tracker ~= f2.B;
    tracker ~= f2.H;
    tracker ~= f2.HorB;
    tracker ~= f2.HorA;
    crosstype = f2;
  }
  /// Decode an input to an (observed) genotype
  auto decode(in string s) {
    return tracker.decode(s);
  }
  auto length() { return tracker.length; }
  string toString() {
    return to!string(tracker);
  }
}

unittest {
  auto f2 = new F2;
  auto tracker = new ObservedGenotypes();
  tracker ~= f2.NA;
  tracker ~= f2.A;
  tracker ~= f2.B;
  tracker ~= f2.H;
  tracker ~= f2.HorB;
  tracker ~= f2.HorA;
  assert(f2.HorB == f2.C);
  assert(tracker.decode("-") == f2.NA);
  assert(tracker.decode("A") == f2.A);
  assert(tracker.decode("B") == f2.B);
  // assert(tracker.decode("H") == f2.H);
  // assert(tracker.decode("C") == f2.HorB);
  writeln(f2.HorA.encoding);
  writeln(tracker);
  assert(tracker.decode("D") == f2.HorA);
}

unittest {
  // the quick way is to used the predefined ObservedF2
  auto tracker = new ObservedF2;
  auto f2 = tracker.crosstype;
  assert(tracker.decode("-") == f2.NA);
  assert(tracker.decode("A") == f2.A);
  assert(tracker.decode("B") == f2.B);
}


/** 
 * Directional F2 with ambiguous scoring:
 * F2  { NA, A, B, AB, BA, HorB, HorA }; 
 */

unittest {
  auto NA = new GenotypeCombinator("NA","-");
  auto A = new GenotypeCombinator("A");
  A ~= new TrueGenotype(0,0);
  auto B  = new GenotypeCombinator("B");
  B ~= new TrueGenotype(1,1);
  auto AB  = new GenotypeCombinator("AB");
  AB ~= new TrueGenotype(0,1);
  auto BA  = new GenotypeCombinator("BA");
  BA ~= new TrueGenotype(1,0);
  auto HorB  = new GenotypeCombinator("HorB","C");
  HorB ~= B;
  HorB ~= AB;
  HorB ~= BA;
  auto HorA  = new GenotypeCombinator("HorA","D");
  HorA ~= A;
  HorA ~= AB;
  HorA ~= BA;
  auto tracker = new ObservedGenotypes();
  tracker ~= NA;
  tracker ~= A;
  tracker ~= B;
  tracker ~= AB;
  tracker ~= BA;
  tracker ~= HorA;
  tracker ~= HorB;
  writeln(HorA);
  assert(HorA.name == "HorA");
  assert(to!string(HorA) == "[(0,0), (0,1), (1,0)]");
}



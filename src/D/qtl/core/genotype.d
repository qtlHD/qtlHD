/**
 * Genotype module
 */

module qtl.core.genotype;

import std.conv;
import std.stdio;
import std.string;
import std.typecons;
import std.exception;
import std.algorithm;
import qtl.core.primitives;

/**

  Genotypes:

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
  symbol tracker named symbols, rather than duplicating them for each marker.
  For this we define ObservedGenotypes, which maintains this list. The symbols
  can be used for a full dataset, or for each marker individually.

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
alias Tuple!(FounderIndex,FounderIndex) Alleles;

/**
 * A true genotype consists of a Tuple of genotypes/alleles - one from each
 * (founder) parent. These genotypes can be phase known/directional 
 */

class TrueGenotype {
  Alleles founders;  // simple Tuple!(i,i)
  this(FounderIndex founder1,FounderIndex founder2) {
    founders = Alleles(founder1, founder2);
  }
  this(TrueGenotype g) { founders = g.founders; }
  // read true genotypes as a comma separated string, e.g. "0,0" or "1,0".
  // make sure not to pass in "" or "NA" - these are simply illegal TrueGenotypes.
  this(in string str) {
    auto fields = split(strip(str),",");
    if (fields.length != 2)
      throw new Exception("Can not parse field with value \"" ~ str ~ "\"");
    auto field1 = to!uint(fields[0]);
    auto field2 = to!uint(fields[1]);
    founders = Alleles(field1,field2);
  }
  bool homozygous()   { return founders[0] == founders[1]; };
  bool heterozygous() { return !homozygous(); };
  // Comparison for consistent sorting of true genotypes, based on founder 
  // genotype numbering scheme
  // e.g. 2,2 > 2,1 > 2,0 > 1,2 > 1,1 > 1,0 > 0,1 > 0,0
  int opCmp(Object other) {
    TrueGenotype _other = cast(TrueGenotype)other;
    auto founders2 = _other.founders;
    if (founders[0] > founders2[0]) return 1;
    if (founders[0] < founders2[0]) return -1;
    if (founders[1] > founders2[1]) return 1;
    if (founders[1] < founders2[1]) return -1;
    return 0;
  }
  // If both founders are equal, the genotype is the same
  const bool opEquals(Object other) {
    auto founders2 = (cast(TrueGenotype)other).founders;
    return founders[0] == founders2[0] && founders[1] == founders2[1];
  }
  // convert to short string, i.e. 1,0 => "1,0"
  const string toTrueGenotype() {
    uint f0 = founders[0]; uint f1 = founders[1];
    return to!string(f0) ~ ',' ~ to!string(f1);
  }
  // string representation of true genotype, i.e. 1,0 -> "(1,0)"
  const string toString() {
    return '(' ~ toTrueGenotype ~ ')';
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

  TrueGenotype[] list;
  string name; // the display name
  Encoding encoding[]; // a list of encoding types

  // initialize a combinator by name. name
  // is also used for input encoding
  this(string name, Encoding code = null) {
    this.name = name;
    add_encoding(name);
    if (code) add_encoding(code);
  }
  // add another encoding string (aliases)
  void add_encoding(Encoding code) { this.encoding ~= code; }
  // see if a true genotype is compatible
  bool match(in TrueGenotype truegen) const { 
    return canFind(cast(TrueGenotype[])list,truegen); 
  }
  // see if an input matches
  bool match(in Encoding code) { return canFind(encoding,code); }
  // uniquely add a genotype to the observed list. Return match.
  TrueGenotype add(TrueGenotype g) {
    foreach(m; list) { if (m.founders == g.founders) return m; }
    list ~= g;
    return g;
  }
  // The ~ operator can add a true genotype to the list
  GenotypeCombinator opOpAssign(string op)(TrueGenotype g) if (op == "~") {
    add(g);
    return this;
  }
  // The ~ operator can add another combinator to the list
  GenotypeCombinator opOpAssign(string op)(GenotypeCombinator c) if (op == "~") {
    foreach(g ; c.list ) { add(g); };
    return this;
  }
  // Return the number of true genotypes available
  const auto length() { return list.length; }
  // Test for equality of two combinators
  bool opEquals(Object other) {
    auto rhs = cast(GenotypeCombinator)other;
    return (list.sort == rhs.list.sort); // probably not the fastest way ;)
  }
  /**
   * Return a string format of the combinator in the form "[(p1,p2)]"
   */
  const string toString() {
    if (isNA) return "[NA]";
    return to!string(list);
  }
  /**
   * Return a string format using the encoder types, e.g, "A B H"
   */
  const string toEncodings() {
    auto a = map!"to!string(a)"(encoding);
    return join(a," ");
  }
  /**
   * Return the primary encoder as a string
   */
  const string toEncoding() {
    return name;
  }
  /**
   * Return the combinator true genotypes as a string, e.g. "0,0 1,0"
   */
  const string toTrueGenotypes() {
    if (list.length == 0) return "None";
    auto a = map!"a.toTrueGenotype"(list);
    return join(a," ");
  }
  const bool isNA() { return length == 0; }
  // compatibility function - deprecate
  auto value() { return this; }
}

alias GenotypeCombinator Gref;  // short name for referencing a combinator

/**
 * ObservedGenotypes tracks all the observed genotypes in a dataset, for the
 * full set, or at a marker position.  This symbol tracker is a convenience
 * class, mostly.
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
  /**
   * Decode an input to an (observed) genotype. This function walks the
   * list of true genotypes and returns the first match.
   *
   * returns NULL on NA.
   */
  GenotypeCombinator decode(in string s) {
    // try to decode by symbol name
    foreach(m; list) { if (m.match(s)) return m; }
    // try to decode by genotype - an empty should return NA
    if (s == "") return decode("NA");
    auto geno = new TrueGenotype(s);
    foreach(m; list) {
      // writeln(m,geno);
      if (!m.isNA) {
        foreach(g ; m.list) {
          if (g == geno) return m;
        }
      }
    }
    throw new Exception("Fail to decode genotype \"" ~ s ~ "\"");
  }
  /// Define =~ to combine combinators
  ObservedGenotypes opOpAssign(string op)(GenotypeCombinator c) if (op == "~") {
    add(c);
    return this;
  }
  auto length() { return list.length; }
  const string toString() {
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
  // create TrueGenotype from string and compare
  assert(new TrueGenotype("0,0") == new TrueGenotype(0,0));
  assert(new TrueGenotype("1,0") == new TrueGenotype(1,0));
  assert(new TrueGenotype("1,2") == new TrueGenotype(1,2));

  // observed genotypes can be any combination of true genotypes
  auto observed1 = new GenotypeCombinator("OBSERVE1");
  observed1 ~= g1;
  observed1 ~= g2;
  observed1 ~= g2; // tests also for duplicate add
  // writeln(observed1.list);
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
  auto symbols = new ObservedGenotypes();
  symbols ~= observed[0];
  symbols ~= observed[1];
  auto test = symbols.add(observed[1]); // add duplicate
  assert(test == observed[1]);
  assert(test == observed1); // No duplication!
  assert(symbols.length == 2);
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

class ObservedRIL {
  RIL crosstype;
  ObservedGenotypes symbols;
  this() {
    auto ril = new RIL;
    symbols = new ObservedGenotypes();
    symbols ~= ril.NA;
    symbols ~= ril.A;
    symbols ~= ril.B;
    crosstype = ril;
  }
  /// Decode an input to an (observed) genotype
  auto decode(in string s) {
    return symbols.decode(s);
  }
  auto length() { return symbols.length; }
  string toString() {
    return to!string(symbols);
  }
}

unittest {
  auto ril = new RIL;
  auto symbols = new ObservedGenotypes();
  symbols ~= ril.NA;
  symbols ~= ril.A;
  symbols ~= ril.B;
  assert(symbols.decode("-") == ril.NA);
  assert(symbols.decode("A") == ril.A);
  assert(symbols.decode("AA") == ril.A);
  assert(symbols.decode("0,0") == ril.A);
  assert(symbols.decode("1,1") == ril.B);
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
  auto symbols = new ObservedGenotypes();
  symbols ~= NA;
  symbols ~= A;
  symbols ~= B;
  GenotypeCombinator ril[];  // create a set of observed genotypes
  // now find them by name
  NA.add_encoding("-"); // also support dash inputs for NA
  ril ~= symbols.decode("-");
  ril ~= symbols.decode("NA");
  ril ~= symbols.decode("A");
  ril ~= symbols.decode("B");
  // ril ~= symbols.decode("C"); // raises exception!
  assert(ril[0].isNA);
  assert(ril[1].isNA);
  assert(ril[2].name == "A");
  assert(to!string(ril[2]) == "[(0,0)]");
  assert(ril[3].name == "B");
  assert(to!string(ril[3]) == "[(1,1)]");
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
  ObservedGenotypes symbols;
  this() {
    auto bc = new BC;
    symbols = new ObservedGenotypes();
    symbols ~= bc.NA;
    symbols ~= bc.A;
    symbols ~= bc.H;
    crosstype = bc;
  }
  /// Decode an input to an (observed) genotype
  auto decode(in string s) {
    return symbols.decode(s);
  }
  auto length() { return symbols.length; }
  string toString() {
    return to!string(symbols);
  }

}

unittest {
  auto bc = new BC;
  auto symbols = new ObservedGenotypes();
  symbols ~= bc.NA;
  symbols ~= bc.A;
  symbols ~= bc.H;
  assert(symbols.decode("-") == bc.NA);
  assert(symbols.decode("A") == bc.A);
  assert(symbols.decode("H") == bc.H);
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
    auto ba = new TrueGenotype(1,0);
    NA = new GenotypeCombinator("NA");
    A  = new GenotypeCombinator("A");
    A ~= aa;
    B  = new GenotypeCombinator("B");
    B ~= bb;
    H  = new GenotypeCombinator("H");
    H ~= ab;
    H ~= ba;
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
  ObservedGenotypes symbols;
  this() {
    auto f2 = new F2;
    symbols = new ObservedGenotypes();
    symbols ~= f2.NA;
    symbols ~= f2.A;
    symbols ~= f2.B;
    symbols ~= f2.H;
    symbols ~= f2.HorB;
    symbols ~= f2.HorA;
    crosstype = f2;
  }
  /// Decode an input to an (observed) genotype
  auto decode(in string s) {
    return symbols.decode(s);
  }
  auto length() { return symbols.length; }
  string toString() {
    return to!string(symbols);
  }
}

unittest {
  auto f2 = new F2;
  auto symbols = new ObservedGenotypes();
  symbols ~= f2.NA;
  symbols ~= f2.A;
  symbols ~= f2.B;
  symbols ~= f2.H;
  symbols ~= f2.HorB;
  symbols ~= f2.HorA;
  assert(f2.HorB == f2.C);
  assert(symbols.decode("-") == f2.NA);
  assert(symbols.decode("A") == f2.A);
  assert(symbols.decode("B") == f2.B);
  // assert(symbols.decode("H") == f2.H);
  // assert(symbols.decode("C") == f2.HorB);
  // writeln(f2.HorA.encoding);
  // writeln(symbols);
  assert(symbols.decode("D")   == f2.HorA);
  assert(symbols.decode("NA")  == f2.NA);
  assert(symbols.decode("")    == f2.NA);
  assert(symbols.decode("1,1") == f2.B);
  assert(symbols.decode("0,1") == f2.H);
  assert(symbols.decode("1,0") == f2.H);
}

unittest {
  // the quick way is to used the predefined ObservedF2
  auto symbols = new ObservedF2;
  auto f2 = symbols.crosstype;
  assert(symbols.decode("-") == f2.NA);
  assert(symbols.decode("A") == f2.A);
  assert(symbols.decode("B") == f2.B);
}


/**
 * Phase known/directional F2 with ambiguous scoring:
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
  auto symbols = new ObservedGenotypes();
  symbols ~= NA;
  symbols ~= A;
  symbols ~= B;
  symbols ~= AB;
  symbols ~= BA;
  symbols ~= HorA;
  symbols ~= HorB;
  assert(HorA.name == "HorA");
  assert(to!string(HorA) == "[(0,0), (0,1), (1,0)]");

  auto gA = new TrueGenotype(0,0);
  auto gB = new TrueGenotype(1,1);
  auto gAB = new TrueGenotype(0,1);
  auto gBA = new TrueGenotype(1,0);

  assert(HorA.match(gA));
  assert(HorA.match(gBA));
  assert(HorA.match(gAB));
  assert(!HorA.match(gB));

  assert(!HorB.match(gA));
  assert(HorB.match(gBA));
  assert(HorB.match(gAB));
  assert(HorB.match(gB));

  // Here a quick test to see if canFind needs sorted input
  int as[] = [ 6,5,4,3,4,7 ];
  assert(canFind(as,4));
  assert(canFind(as,7));
}

/**
 * Flex implementation
 */

/*
class Flex {
  GenotypeCombinator NA;
  this() {
    NA = new GenotypeCombinator("NA");
    NA.add_encoding("-");
  }
}

class ObservedFlex {
  Flex crosstype;
  ObservedGenotypes symbols;
  this() {
    auto flex = new Flex;
    symbols = new ObservedGenotypes();
    symbols ~= flex.NA;
    crosstype = flex;
  }
  /// Decode an input to an (observed) genotype
  auto decode(in string s) {
    return symbols.decode(s);
  }
  auto length() { return symbols.length; }
  string toString() {
    return to!string(symbols);
  }

}
*/

/**
 * Structure for reading encoded CSV files. In above examples the Genotype is
 * predefined at compile time. Here we define the genotype at run time.  The
 * encoding can be in a separate file, or within file. An encoding, numbers
 * referring to founders, simply reads:
 *
 * GENOTYPE A as 0,0
 * GENOTYPE BB as 1,1
 * GENOTYPE C,CC as 2,2         # C and CC both represent 2,2
 * GENOTYPE AB as 1,0           # 1,0 phase known/directional
 * GENOTYPE BA as 0,1           # 0,1 phase known/directional
 * GENOTYPE AC,CA as 0,2 2,0    # 2,0 phase unknown/not directional
 * GENOTYPE CA,F as 2,0 0,2     # F is an alias for AC,CA
 * GENOTYPE AorB as 0,0 1,1
 * GENOTYPE AorABorAC as 0,0 1,0 0,1 0,2 2,0
 * GENOTYPE NA,- as None
 *
 * and should be in the header, close to the start, of the file. The
 * identifier (A,B, etc) can be any string.
 *
 * Note: duplicates are undefined, so don't do
 *
 * GENOTYPE B as 1,1
 * GENOTYPE BB as 1,1
 *
 * but
 *
 * GENOTYPE B,BB as 1,1
 */

/**
 * Fetch an observed genotype definition from a string.
 * Returns name (aliases) and true genotypes. The input string
 * is of format
 *
 *   GENOTYPE name1,name2[,...] as uint,uint uint,uint [...]
 *
 * On error this function will throw an exception.
 *
 * See genotype.d for more information
 *
 *
 */

Tuple!(string[],TrueGenotype[]) parse_observed_genotype_string(string s)
out(result) {
  auto names = result[0];
  auto tgs = result[0];
  assert(names.length >= 1);
  assert(tgs.length >= 1);
}
body {
  auto tokens = split(s," ");
  if (tokens[0] != "GENOTYPE")
    throw new Exception("Expected GENOTYPE for " ~ s);
  int i = 0;
  string names[];
  foreach(name ; tokens[1..$]) {
    if (name == "as") break;
    foreach(n2 ; split(name,",")) {
      names ~= n2;
    }
    i++;
  }
  if (i > tokens.length-2)
    throw new Exception("Expected 'as' for " ~ s);
  TrueGenotype[] tgs;
  foreach(tt ; tokens[i+2..$]) {
    // writeln("<",tt,">");
    if (tt == "#" || tt == "" || tt == "None") break;
    auto alleles = split(tt,",");
    if (alleles.length != 2)
      throw new Exception("Malformed genotype in " ~ s);
    uint a1 = to!uint(alleles[0]);
    uint a2 = to!uint(alleles[1]);
    tgs ~= new TrueGenotype(a1,a2);
  }
  return tuple(names, tgs);
}

/**
 * Convenience object for parsing observed genotype definition
 * strings. The result is stored in names and genotypes arrays.
 */
class EncodedGenotype {
  string names[];
  TrueGenotype genotypes[];
  this(string s) {
    auto tuple = parse_observed_genotype_string(s);
    names = tuple[0];
    genotypes = tuple[1];
  }
  GenotypeCombinator combinator() {
    return null;
  }
}

unittest {
  auto symbols = new ObservedGenotypes();
  auto tuple = parse_observed_genotype_string("GENOTYPE A as 0,0");
  assert(tuple[0] == ["A"]);
  assert(to!string(tuple[1]) == "[(0,0)]");
  assertThrown(parse_observed_genotype_string("GENOTYPE A as (0,0)"));

  auto eg = new EncodedGenotype("GENOTYPE A as 0,0");
  symbols.add(eg.combinator);
  assert(eg.names == ["A"]);
  assert(to!string(eg.genotypes) == "[(0,0)]");
  eg = new EncodedGenotype("GENOTYPE AB as 1,0           # 1,0 phase known/directional");
  assert(eg.names == ["AB"]);
  assert(to!string(eg.genotypes) == "[(1,0)]", to!string(eg.genotypes));
  eg = new EncodedGenotype("GENOTYPE AC,CA as 0,2 2,0    # phase unknown / not directional");
  assert(eg.names == ["AC","CA"]);
  assert(to!string(eg.genotypes) == "[(0,2), (2,0)]", to!string(eg.genotypes));
  eg = new EncodedGenotype("GENOTYPE NA,- as None  ");
  assert(eg.names == ["NA","-"]);
  assert(to!string(eg.genotypes) == "[]", to!string(eg.genotypes));
  eg = new EncodedGenotype("GENOTYPE A AA as 0,0");
  assert(eg.names == ["A","AA"]);
  eg = new EncodedGenotype("GENOTYPE AorABorAC as 0,0 1,0 0,1 0,2 2,0");
  assert(eg.names == ["AorABorAC"]);
  assert(to!string(eg.genotypes) == "[(0,0), (1,0), (0,1), (0,2), (2,0)]", to!string(eg.genotypes));
  // faulty ones are more interesting
  assertThrown(new EncodedGenotype("GENOTYPE A as (0,0)"));
  assertThrown(new EncodedGenotype("ENOTYPE A as 0,0"));
  assertThrown(new EncodedGenotype(" GENOTYPE A as 0,0"));
  assertThrown(new EncodedGenotype("GENOTYPE A 0,0"));
  assertThrown(new EncodedGenotype("GENOTYPE A 0"));
  assertThrown(new EncodedGenotype("GENOTYPE A 0,A"));

  // now start decoding
 * GENOTYPE A as 0,0
 * GENOTYPE BB as 1,1
 * GENOTYPE C,CC as 2,2         # C and CC both represent 2,2
 * GENOTYPE AB as 1,0           # 1,0 phase known/directional
 * GENOTYPE BA as 0,1           # 0,1 phase known/directional
 * GENOTYPE AC,CA as 0,2 2,0    # 2,0 phase unknown/not directional
 * GENOTYPE CA,F as 2,0 0,2     # F is an alias for AC,CA
 * GENOTYPE AorB as 0,0 1,1
 * GENOTYPE AorABorAC as 0,0 1,0 0,1 0,2 2,0
 * GENOTYPE NA,- as None
  assert(to!string(symbols.decode("NA")) == "(NA)");
  assert(symbols.decode("-") == symbols.decode("NA"));
  assert(to!string(symbols.decode("A")) == "[(0,0)]");
  assert(to!string(symbols.decode("C")) == "[(2,2)]");
  assert(symbols.decode("CC") == symbols.decode("C"));
  assert(symbols.decode("BB") == symbols.decode("B"); 
  assert(symbols.decode("BA") == cross.combinator["BA"]);
  assert(to!string(symbols.decode("AB")) == "[(1,0)]"); // phase known!
  assert(to!string(symbols.decode("BA")) == "[(0,1)]"); // phase known!
  assert(to!string(symbols.decode("AorABorAC")) == "[(0,0),(1,0),(0,1),(0,2),(2,0)]");
}

/**
 * Takes a list of genotype encoding string lines, e.g.
 *
 * ["GENOTYPE NA,- as None","GENOTYPE A as 0,0","GENOTYPE B,BB as 1,1"]
 *
 * and turns it into an ObservedGenotypes object
 */

/*
class EncodedCross {
  Tuple!(string,GenotypeCombinator) combinator[];
  this(string list[]) {
    // parse list
    foreach(line ; list) {
       // writeln(line);
       if (line.strip() == "") continue;
       add(line);
    }
  }

  EncodedGenotype add(string line) {
    auto line_item = new EncodedGenotype(line);
    auto n = line_item.names[0];
    if (n in combinator)
      throw new Exception("Duplicate " ~ line);

    combinator[n] = new GenotypeCombinator(n);
    foreach (tt ; line_item.genotypes) {
       combinator[n] ~= tt;
    }
    foreach (n_alias ; line_item.names[1..$]) {
       combinator[n].add_encoding(n_alias);
    }
    // writeln("--->",combinator[n].encoding,combinator[n]);
    return line_item;
  }

}

unittest {
  auto encoded = "
GENOTYPE NA,- as None
GENOTYPE A as 0,0
GENOTYPE B,BB as 1,1
GENOTYPE C,CC as 2,2         # alias
GENOTYPE AB as 1,0           # phase known/directional
GENOTYPE BA as 0,1           # phase known/directional
GENOTYPE AC,CA as 0,2 2,0    # phase unknown/not directional
GENOTYPE CA,F as 2,0 0,2     # alias
GENOTYPE AorB as 0,0 1,1
GENOTYPE AorABorAC as 0,0 1,0 0,1 0,2 2,0";

  auto cross = new EncodedCross(split(encoded,"\n"));
  auto symbols = new ObservedGenotypes();
  // add genotypes to symbols
  foreach (legaltype; cross.combinator) {
    writeln(legaltype);
    symbols ~= legaltype;
  }
  assert(symbols.decode("NA") == cross.combinator["NA"]);
  assert(symbols.decode("-") == cross.combinator["NA"]);
  assert(symbols.decode("A") == cross.combinator["A"]);
  assert(symbols.decode("B") == cross.combinator["B"]);
  assert(symbols.decode("CC") == cross.combinator["C"]);
  // assert(symbols.decode("BB") == cross.combinator["BB"]); <- error
  assert(symbols.decode("BB") == cross.combinator["B"]); //  <- OK
  assert(symbols.decode("AB") == cross.combinator["AB"]);
  assert(symbols.decode("BA") == cross.combinator["BA"]);
  assert(symbols.decode("AorB") == cross.combinator["AorB"]);
  assert(symbols.decode("0,0") == cross.combinator["A"]);
  assert(symbols.decode("1,1") == cross.combinator["B"]);
  assert(symbols.decode("1,0") == cross.combinator["AB"]);
  writeln("Cross symbols: ",symbols);
  assert(symbols.decode("0,1") == cross.combinator["BA"], to!string(symbols.decode("0,1")) ~ " expected " ~ to!string(cross.combinator["BA"]));
  // assert(symbols.decode("0,0|1,1") == cross.combinator["AorB"]);
  // writeln(cross["AorB"].encoding);
  // writeln(symbols);
  // Test for duplicates
  encoded = "
GENOTYPE A as 0,0
GENOTYPE A as 1,1
";
  assertThrown(new EncodedCross(split(encoded,"\n")));
  // This is legal, though probably undefined
  encoded = "
GENOTYPE A as 0,0
GENOTYPE AA as 0,0
";
  assertNotThrown(new EncodedCross(split(encoded,"\n")));
}
*/

// CrossType : allowable cross types we will handle
enum CrossType { BC, F2, RILself, RILsib };

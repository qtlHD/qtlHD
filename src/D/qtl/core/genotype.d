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

alias GenotypeCombinator[][] GenotypeMatrix; // = new double[][][](n_markers,n_ind);

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

  X chromosome (vs autosome) should be queried at the chromosome level.

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
    if (str == "NA" || str == "" || str == "-")
      throw new Exception("Can not initialize genotype field with value \"" ~ str ~ "\"");
    auto fields = split(strip(str),",");
    if (fields.length != 2)
      throw new Exception("Can not parse true genotype field with value \"" ~ str ~ "\"");
    auto field1 = to!uint(fields[0]);
    auto field2 = to!uint(fields[1]);
    founders = Alleles(field1,field2);
  }
  const bool homozygous()   { return founders[0] == founders[1]; };
  const bool heterozygous() { return !homozygous(); };
  TrueGenotype reversed() { return new TrueGenotype(founders[1],founders[0]); }
  // Comparison for consistent sorting of true genotypes, based on founder 
  // genotype numbering scheme
  // e.g. 2,2 > 2,1 > 2,0 > 1,2 > 1,1 > 1,0 > 0,1 > 0,0
  override int opCmp(Object other) {
    TrueGenotype _other = cast(TrueGenotype)other;
    auto founders2 = _other.founders;
    if (founders[0] > founders2[0]) return 1;
    if (founders[0] < founders2[0]) return -1;
    if (founders[1] > founders2[1]) return 1;
    if (founders[1] < founders2[1]) return -1;
    return 0;
  }
  // If both founders are equal, the genotype is the same
  override const bool opEquals(Object other) {
    auto founders2 = (cast(TrueGenotype)other).founders;
    return founders[0] == founders2[0] && founders[1] == founders2[1];
  }
  // convert to short string, i.e. 1,0 => "1,0"
  const string toTrueGenotype() {
    uint f0 = founders[0]; uint f1 = founders[1];
    return to!string(f0) ~ ',' ~ to!string(f1);
  }
  // string representation of true genotype, i.e. 1,0 -> "(1,0)"
  override const string toString() {
    return '(' ~ toTrueGenotype ~ ')';
  }
  const FounderIndex get_allele(in uint allele_index) {
    if(allele_index) return(this.founders[1]);
    else return(this.founders[0]);
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
  Encoding encoding[]; // a list of encoding types (strings)
  bool phase_known = true;

  // initialize a combinator by name. name
  // is also used for input encoding
  this(string name, Encoding code = null, bool phase_known = true) {
    this.name = name;
    add_encoding(name);
    if (code) add_encoding(code);
    this.phase_known = phase_known;
  }
  this(string name, TrueGenotype g) {
    this.name = name;
    add(g);
  }
  this(string name, TrueGenotype gs[]) {
    this.name = name;
    foreach(g ; gs) {
      add(g);
    }
  }
  // initialize a combinator with lists
  this(string names[], TrueGenotype gs[]) {
    this.name = names[0];
    add_encoding(name);
    foreach(n ; names[1..$]) { 
      add_encoding(n);
    }
    foreach(g ; gs) {
      add(g);
    }
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
  // If phase_unknown is true, the method will also add the reversed
  // founders to the list.
  TrueGenotype add(TrueGenotype g) {
    // writeln("Adding ",g);
    foreach(m; list) { if (m.founders == g.founders) return m; }
    list ~= g;
    if (!phase_known) add(g.reversed);
    return g;
  }
  // The ~ operator can add a true genotype to the list
  GenotypeCombinator opOpAssign(string op)(TrueGenotype g) if (op == "~") {
    add(g);
    return this;
  }
  // The ~ operator can add another combinator to the list
  GenotypeCombinator opOpAssign(string op)(GenotypeCombinator c) if (op == "~") {
    foreach(g ; c.list ) { add(g); }
    return this;
  }
  // Return the number of true genotypes available
  const auto length() { return list.length; }
  // Test for equality of two combinators
  override bool opEquals(Object other) {
    auto rhs = cast(GenotypeCombinator)other;
    return (list.sort == rhs.list.sort); // probably not the fastest way ;)
  }
  /**
   * Return a string format of the combinator in the form "[(p1,p2)]"
   */
  override const string toString() {
    if (isNA) return "[NA]";
    return to!string(list);
  }
  /**
   * Return a string format using the encoder types, e.g, "A B H"
   */
  const string toEncodingString() {
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
   * Decode an input to an (observed) genotype. This function first walks the
   * list of true genotype *names*, and next the actual genotypes. 
   * Returns the first matching combinator.
   *
   * returns NULL on NA.
   */
  GenotypeCombinator decode(in string s) {
    // try to decode by symbol name
    foreach(m; list) { if (m.match(s)) return m; }
    // try to decode by genotype - an empty should return NA
    if (s == "") return decode("NA");
    try {
      auto geno = new TrueGenotype(s);
      foreach(m; list) {
        // writeln(m,geno);
        if (!m.isNA) {
          foreach(g ; m.list) {
            if (g == geno) return m;
          }
        }
      }
    }
    catch (Exception e) {
      throw new Exception("Failed to decode genotype \"" ~ s ~ "\"");
    }
    throw new Exception("Failed to decode genotype \"" ~ s ~ "\"");
  }
  /// Define =~ to combine combinators
  ObservedGenotypes opOpAssign(string op)(GenotypeCombinator c) if (op == "~") {
    add(c);
    return this;
  }
  auto length() { return list.length; }
  override const string toString() { return to!string(list); }
}

/**
 * New style genotypes - support for observed genotypes
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
    writeln(names);
    writeln(genotypes);
    return new GenotypeCombinator(names,genotypes);
  }
}

/**
 * Convert a genotype matrix, in string representation, to a genotype combinator
 * matrix
 */

GenotypeMatrix convert_to_combinator_matrix(string[][] g,ObservedGenotypes observed) {
  GenotypeMatrix gm;
  foreach(row ; g) {
    auto nrow = new GenotypeCombinator[](row.length);
    foreach(i, symbolname; row) {
      nrow[i] = observed.decode(symbolname); 
    }
    gm ~= nrow;
  }

  return gm;
}

// omit individuals from genotype matrix
GenotypeCombinator[][] omit_ind_from_genotypes(GenotypeCombinator[][] geno, bool[] to_omit)
{
  if(geno.length != to_omit.length)
    throw new Exception("no. individuals in geno (" ~ to!string(geno.length) ~
                        ") doesn't match length of to_omit (" ~ to!string(to_omit.length) ~ ")");

  GenotypeCombinator[][] ret;

  foreach(i; 0..to_omit.length) {
    if(!to_omit[i])
      ret ~= geno[i];
  }

  return ret;
}

/**
 * Parse the genotype_id string and turn it into a combinator
 *
 *   F2  : "A H B D C" where D is HorB and C is HorA
 *   BC  : "A H"
 *   RIL : "A B"
 *   NA  : "- NA"
 */

ObservedGenotypes parse_genotype_ids(string cross,string genotype_ids,string na_ids) {
  auto symbols = new ObservedGenotypes();
  if (cross == "F2") {
    auto aa = new TrueGenotype(0,0);
    auto bb = new TrueGenotype(1,1);
    auto ab = new TrueGenotype(0,1);
    auto ba = new TrueGenotype(1,0);
    auto na_ids2 = split(na_ids," ");
    auto NA = new GenotypeCombinator(na_ids2,null);
    auto ids = split(genotype_ids," ");
    auto A  = new GenotypeCombinator(ids[0],aa);
    auto H  = new GenotypeCombinator(ids[1],[ab,ba]);
    auto B  = new GenotypeCombinator(ids[2],bb);
    auto D  = new GenotypeCombinator(ids[3],[ab,aa]);
    auto C  = new GenotypeCombinator(ids[4],[ab,bb]);
    symbols ~= NA;
    symbols ~= A;
    symbols ~= B;
    symbols ~= H;
    symbols ~= D;
    symbols ~= C;
  }
  return symbols;
}


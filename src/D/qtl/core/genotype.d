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

  To introduce sex chromosomes, we can set one type to a special value (say -1
  for the missing counterpart of X). The actual type in human would be [1,-1].
  Sex should be queried at the chromosome level.

  Note that outbreeding can be supported by creating artificial founders
  consisting of singular genotypes (similar to a RIL), with one maternal and
  one paternal genotype.

  With SNP markers, two representations can be supported using this genotyping
  system. In the first, each SNP is treated as a marker - therefore there are
  always two genotypes involved. This is the simplest approach. The second
  approach works at the allele level. SNPs are combined into ORFs, or genes.
  Multiple allele combinations can be supported with this genotyping
  representation.

 */


/**
 * The follwing is the enum typed genotypes - which are not as
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
/**
 * New style genotypes - support for observed genotypes
 */

import std.typecons; 

unittest {
  writeln("Unit test" ~ __FILE__);
  // at this point founders are simply numbers
  FounderIndex[] founder = [ 1, 2, 3, 4, 5 ];
  // create a few genotypes (at a marker location)
  auto g1 = new TrueGenotype(founder[0],founder[2]);
  auto g2 = new TrueGenotype(founder[1],founder[1]); // could be a RIL
  auto g3 = new TrueGenotype(founder[2],founder[3]);
  auto g4 = new TrueGenotype(founder[4],founder[3]);
  assert(g1.heterozygous());
  assert(g2.homozygous());
  // observed genotypes can be any combination of true genotypes
  GenotypeCombinator observed_combi1;
  observed_combi1 ~= g1;
  observed_combi1 ~= g2;
  GenotypeCombinator observed_combi2;
  observed_combi2 ~= g2;
  // observed_combiX already acts as an index, so now we only need 
  // to store the observed genotypes with the marker. The marker/genotype
  // matrix stores references to observed_combiX.
}

/**
 * Enum style genotypes
 */

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

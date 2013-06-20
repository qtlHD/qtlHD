import std.typecons;
import std.stdio;
import std.string;
import std.exception;
import std.conv;
import qtl.core.genotype;

/**
 * Unit tests for TrueGenotype and GenotypeCombinators
 */

unittest {
  writeln("Unit test " ~ __FILE__);
  // at this point founders are simply numbers
  FounderIndex[] founder = [ 1, 2, 3, 4, 5 ];
  // create a few genotypes (at a marker location)
  auto g1  = new TrueGenotype(founder[0],founder[2]);
  auto g2  = new TrueGenotype(founder[1],founder[1]); // could be a RISELF
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
  symbols.add(eg.combinator);
  assert(eg.names == ["AB"]);
  assert(to!string(eg.genotypes) == "[(1,0)]", to!string(eg.genotypes));
  eg = new EncodedGenotype("GENOTYPE BA as 0,1           # 1,0 phase known/directional");
  symbols.add(eg.combinator);
  eg = new EncodedGenotype("GENOTYPE AC,CA as 0,2 2,0    # phase unknown / not directional");
  symbols.add(eg.combinator);
  assert(eg.names == ["AC","CA"]);
  assert(to!string(eg.genotypes) == "[(0,2), (2,0)]", to!string(eg.genotypes));
  eg = new EncodedGenotype("GENOTYPE C,CC as 2,2");
  symbols.add(eg.combinator);
  // eg = new EncodedGenotype("GENOTYPE NA,- as None  ");
  // assert(eg.names == ["NA","-"]);
  // assert(to!string(eg.genotypes) == "[]", to!string(eg.genotypes));
  eg = new EncodedGenotype("GENOTYPE A AA as 0,0");
  symbols.add(eg.combinator);
  assert(eg.names == ["A","AA"]);
  eg = new EncodedGenotype("GENOTYPE AorABorAC as 0,0 1,0 0,1 0,2 2,0");
  symbols.add(eg.combinator);
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
  // assert(to!string(symbols.decode("NA")) == "(NA)");
  // assert(symbols.decode("-") == symbols.decode("NA"));
  foreach (symbol ; symbols.list) {
    writeln(symbol.toEncodingString," --> ",symbol);
  }
  assert(to!string(symbols.decode("A")) == "[(0,0)]");
  assert(to!string(symbols.decode("C")) == "[(2,2)]");
  assert(symbols.decode("CC") == symbols.decode("C"));
  assert(to!string(symbols.decode("AB")) == "[(1,0)]"); // phase known!
  assert(to!string(symbols.decode("BA")) == "[(0,1)]"); // phase known!
  assert(to!string(symbols.decode("AorABorAC")) == "[(0,0), (0,1), (0,2), (1,0), (2,0)]");
}



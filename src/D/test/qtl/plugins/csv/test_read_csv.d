/**
 * Plugin for reading simple CSV files
 *
 * Note: this implementation reads everything into an object. We may switch to a Tuple
 * implementation later.
 */

module test.qtl.plugins.csv.test_read_csv;

import core.memory;
import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.typecons;

import qtl.core.primitives;
import qtl.core.chromosome;
// import qtl.core.util.matrix;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.core.individual;
import example.genotype_examples;
import qtl.plugins.csv.read_csv;

// ==== UNIT TESTS ===================================================

unittest {
  writeln("Unit test " ~ __FILE__);
  alias std.path.buildPath buildPath;
  auto fn = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","..","test","data","input","listeria.csv"));
  writeln("  - reading CSV " ~ fn);
  Marker m2 = new Marker(4.8);
  auto markers = [ m2 ];
  auto data = new ReadSimpleCSV!(F2,ObservedF2)(fn);
  auto F2 = data.crosstype;
  assert(data.markers.length == 133, to!string(data.markers.length));
  assert(data.phenotypenames[0] == "T264");
  assert(data.markers[0].name == "D10M44");
  assert(data.markers[0].id == 0);
  assert(data.markers[0].chromosome.name == "1",data.markers[0].chromosome.name);
  assert(data.markers[1].id == 1);
  // Check chromosomes
  assert(data.chromosomes.length == 20, to!string(data.chromosomes.length));
  assert(data.chromosomes["X"].id == ID_UNKNOWN);
  assert(data.chromosomes["7"].id == 7);
  assert(data.markers[2].position == 24.84773, "Marker position not matching");
  // Check phenotype
  assert(isNA(data.phenotypes[29][0]), to!string(data.phenotypes[29][0].value));
  assert(data.phenotypes[30][0].value == 74.417);
  // Check genotype (hard coded)
  assert(data.genotypecombinator[1][0] == F2.NA);
  assert(data.genotypecombinator[1][1] == F2.B);
  // This should also work
  assert(data.genotypecombinator[1][0] == data.symbols.decode("NA"));
  assert(data.genotypecombinator[1][1] == data.symbols.decode("B"));
  assert(data.individuals.length == 120);

  // foreach(name; data.chromosomes.keys.sort) {
  // }
}

unittest {
  alias std.path.buildPath buildPath;
  auto fn = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","..","test","data","input","hyper.csv"));
  writeln("  - reading CSV " ~ fn);
  auto data = new ReadSimpleCSV!(F2,ObservedF2)(fn);
  assert(data.markers.length == 174, to!string(data.markers.length));
  assert(data.phenotypenames == ["bp", "sex"]);
  assert(data.markers[3].name == "D1Mit178");
  assert(data.markers[3].id == 3);
  assert(data.markers[4].id == 4);
  // Check chromosomes
  assert(data.chromosomes.length == 20, to!string(data.chromosomes.length));
  assert(data.chromosomes["X"].id == ID_UNKNOWN);
  assert(data.chromosomes["7"].id == 7);
  assert(data.markers[2].position == 32.8000000002, "Marker position not matching");
  // Check phenotype
  assert(data.phenotypes[29][0].value == 116.3);
  assert(data.phenotypes[30][0].value == 110.2);
  assert(data.phenotypes[29][1].value == 1);
  assert(data.phenotypes[30][1].value == 1);
  // Check genotype
  auto F2 = data.crosstype;
  assert(data.genotypecombinator[1][0] == F2.H);
  assert(data.genotypecombinator[1][1] == F2.H);
  assert(data.genotypecombinator[2][3] == F2.NA);
  assert(data.genotypecombinator[2][4] == F2.B);
}

unittest {
  alias std.path.buildPath buildPath;
  auto fn = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","..","test","data","input","hyper_noX.csv"));
  writeln("  - reading CSV " ~ fn);
  auto data = new ReadSimpleCSV!(BC,ObservedBC)(fn);
  assert(data.markers.length == 170, to!string(data.markers.length));
  assert(data.phenotypenames == ["bp", "sex"]);
  assert(data.markers[3].name == "D1Mit178");
  assert(data.markers[3].id == 3);
  assert(data.markers[4].id == 4);
  // Check chromosomes
  assert(data.chromosomes.length == 19, to!string(data.chromosomes.length));
  assert(data.chromosomes["7"].id == 7);
  assert(data.markers[2].position == 32.8000000002, "Marker position not matching");
  // Check phenotype
  assert(data.phenotypes[29][0].value == 116.3);
  assert(data.phenotypes[30][0].value == 110.2);
  assert(data.phenotypes[29][1].value == 1);
  assert(data.phenotypes[30][1].value == 1);
  // Check genotype
  auto BC = data.crosstype;
  assert(data.genotypecombinator[1][0] == BC.H);
  assert(data.genotypecombinator[1][1] == BC.H);
  assert(data.genotypecombinator[2][3] == BC.NA);
  assert(data.genotypecombinator[2][4] == BC.A);
}

/**
 * Unit tests for the Flex cross. A flex cross can (potentially) read a data
 * file without assuming what is in the file. Observed genotypes are added to
 * the symbols beforehand, as the parser has to know what an A, B or H symbol
 * means.
 */

unittest {
  writeln("Unit test " ~ __FILE__);
  auto encoded = "
GENOTYPE NA,- as None        
GENOTYPE A as 0,0
GENOTYPE B as 1,1
GENOTYPE H as 0,1
GENOTYPE C as 0,0 0,1
GENOTYPE D as 1,1 0,1";

  auto types = new EncodedCross(split(encoded,"\n"));
  auto observed = new ObservedFlex();
  // add genotypes to symbols
  foreach (legaltype; types.combinator) {
    observed.symbols ~= legaltype;
  }
  alias std.path.buildPath buildPath;
  auto fn = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","..","test","data","input","listeria.csv"));
  writeln("  - reading CSV (for Flex) " ~ fn);
  auto data = new ReadSimpleCSV!(Flex,ObservedFlex)(fn, observed);
  auto cross = data.crosstype;
  assert(data.markers.length == 133, to!string(data.markers.length));
  assert(data.phenotypenames[0] == "T264");
  assert(data.markers[0].name == "D10M44");
  assert(data.markers[0].id == 0);
  assert(data.markers[0].chromosome.name == "1",data.markers[0].chromosome.name);
  assert(data.markers[1].id == 1);
  // Check chromosomes
  assert(data.chromosomes.length == 20, to!string(data.chromosomes.length));
  assert(data.chromosomes["X"].id == ID_UNKNOWN);
  assert(data.chromosomes["7"].id == 7);
  assert(data.markers[2].position == 24.84773, "Marker position not matching");
  // Check phenotype
  assert(isNA(data.phenotypes[29][0]), to!string(data.phenotypes[29][0].value));
  assert(data.phenotypes[30][0].value == 74.417);
  // Check genotype
  assert(data.genotypecombinator[1][0] == data.symbols.decode("NA"));
  assert(data.genotypecombinator[1][1] == data.symbols.decode("B"));
  assert(data.individuals.length == 120);
}



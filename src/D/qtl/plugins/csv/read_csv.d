/**
 * Plugin for reading simple CSV files
 *
 * Note: this implementation reads everything into an object. We may switch to a Tuple
 * implementation later.
 */

module qtl.plugins.csv.read_csv;

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
import qtl.plugins.qtab.read_qtab : InputInfo;


/**
 * Read a simple CSV file containing marker names, chromosome nrs, position,
 * phenotype and genotype - such as the listeria.csv file used in R/qtl.
 *
 * The cross type is injected as XType. An XType works as long as it is a known
 * Genotype class (e.g. BC, F2, RIL, Flex). Observed genotypes are handled as
 * 'fixed' in this implementation. That is, we already know the observed types
 * before reading the genotype data file. The symbols hands back matching
 * (observed) genotypes.
 *
 * The file is parsed once on class instantiation. Elements can be queried.
 */

class ReadSimpleCSV(XType,ObservedXType) {
  private File f;
  XType crosstype;
  ObservedXType symbols;
  string[] phenotypenames;
  Marker[] markers;
  Individuals individuals;
  Chromosome[string] chromosomes;
  Phenotype!double[][] phenotypes;
  GenotypeCombinator[][] genotypes;
  size_t n_phenotypes;

  this(in string fn, ObservedXType observed = null) {
    f = File(fn,"r");
    scope(exit) f.close();
    if (observed is null)
      symbols = new ObservedXType;
    else
      symbols = observed;
    writeln(symbols);
    // symbols = new ObservedXType; // this may become a parameter
    individuals = new Individuals();
    // read markers
    Marker[] ms;
    auto header = split(f.readln(),",");
    immutable n_columns = header.length;
    ms.reserve(n_columns);  // pre-allocate memory (just good practise)

    // read chromosome info
    auto cnames = split(f.readln(),",");
    if (cnames.length != n_columns) throw new Exception("# Chromosomes out of range");

    // find first non-blank (indicating last phenotype columns)

    foreach (i, cname; cnames)
    {
      if(cname != "") {
        n_phenotypes = i;
        break;
      }
    }

    // create Marker objects; have m.id be column within phenotypes or within genotypes
    foreach (i, mname; header)
    {
      Marker m = new Marker(MARKER_POSITION_UNKNOWN, mname, cast(uint)(i-n_phenotypes));
      ms ~= m;
    }

    // split into phenotypes and markers
    phenotypenames = new string[n_phenotypes];
    foreach (i; 0..n_phenotypes) {
      phenotypenames[i] = ms[i].name;
    }
    markers = ms[n_phenotypes..$];
    cnames = cnames[n_phenotypes..$];

    foreach (i, cname; cnames)
    {
      immutable cname2 = strip(cname);
      if (!(cname in chromosomes))
        chromosomes[cname2] = get_chromosome_with_id(cname2);
      markers[i].chromosome = chromosomes[cname2];
    }

    // read position info
    auto positions = split(f.readln(),",");
    if (positions.length != n_columns) throw new Exception("Positions out of range");
    positions = positions[n_phenotypes..$];
    foreach (i, pos; positions)
    {
      immutable pos2 = strip(pos);
      markers[i].position = to!double(pos2);
    }
    // read remaining data
    string buf;
    uint n_individual = 0;
    while (f.readln(buf)) {
      n_individual++;
      auto fields = split(buf,",");
      if (fields.length != n_columns) throw new Exception("Field # out of range in ", buf);
      Phenotype!double[] ps;
      ps.reserve(n_columns);
      foreach(field; fields[0..n_phenotypes]) {
        ps ~= set_phenotype!double(field);
      }
      phenotypes ~= ps;
      // set genotype
      crosstype = new XType;
      GenotypeCombinator[] gs;
      // we use the predefined crosstype symbols
      gs.reserve(n_columns);  // pre-allocate memory (just good practise)
      foreach (field; fields[n_phenotypes..$]) {
        gs ~= symbols.decode(strip(field));
      }
      genotypes ~= gs;
      individuals.list ~= new Individual(n_individual);
    }
  }
}


InputInfo load_csv(string fn) {
  SymbolSettings s;
  Founders f;
  Marker[] ms;
  Inds i;
  PhenotypeMatrix p;
  string[][] g;
  ObservedGenotypes observed;
  foreach (fn ; fns) {
    auto res = load_qtab(fn);
    auto t = res[0].get!QtabFileType;
    auto d = res[1];
    with (QtabFileType) {
      switch(t) {
        case symbols: 
          s = d.get!SymbolSettings; 
          observed = res[2].get!ObservedGenotypes;
          break;
        case founder: f = d.get!Founders; writeln(f); break;
        case location: 
          ms = d.get!(Marker[]); 
          writeln(ms);
          break;
        case genotype: 
          i = d.get!Inds;
          g = res[2].get!(string[][]);
          break;
        case phenotype: 
          // auto pids = d; ignored, for now
          p = res[2].get!PhenotypeMatrix; break;
        default: throw new Exception("Unsupported file type for " ~ fn);
      }
    }
  }
  // Turn the genotype matrix into a genotype combinator matrix
  auto gc = convert_to_combinator_matrix(g,observed);
  return tuple(s,f,ms,i,p,observed,gc);
}

// ==== UNIT TESTS ===================================================

unittest {
  writeln("Unit test " ~ __FILE__);
  alias std.path.buildPath buildPath;
  auto fn = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","test","data","input","listeria.csv"));
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
  assert(data.phenotypes[29][0].value == PHENOTYPE_NA, to!string(data.phenotypes[29][0].value));
  assert(data.phenotypes[30][0].value == 74.417);
  // Check genotype (hard coded)
  assert(data.genotypes[1][0] == F2.NA);
  assert(data.genotypes[1][1] == F2.B);
  // This should also work
  assert(data.genotypes[1][0] == data.symbols.decode("NA"));
  assert(data.genotypes[1][1] == data.symbols.decode("B"));
  assert(data.individuals.length == 120);

  // foreach(name; data.chromosomes.keys.sort) {
  // }
}

unittest {
  alias std.path.buildPath buildPath;
  auto fn = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","test","data","input","hyper.csv"));
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
  assert(data.genotypes[1][0] == F2.H);
  assert(data.genotypes[1][1] == F2.H);
  assert(data.genotypes[2][3] == F2.NA);
  assert(data.genotypes[2][4] == F2.B);
}

unittest {
  alias std.path.buildPath buildPath;
  auto fn = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","test","data","input","hyper_noX.csv"));
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
  assert(data.genotypes[1][0] == BC.H);
  assert(data.genotypes[1][1] == BC.H);
  assert(data.genotypes[2][3] == BC.NA);
  assert(data.genotypes[2][4] == BC.A);
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
  auto fn = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","test","data","input","listeria.csv"));
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
  assert(data.phenotypes[29][0].value == PHENOTYPE_NA, to!string(data.phenotypes[29][0].value));
  assert(data.phenotypes[30][0].value == 74.417);
  // Check genotype
  assert(data.genotypes[1][0] == data.symbols.decode("NA"));
  assert(data.genotypes[1][1] == data.symbols.decode("B"));
  assert(data.individuals.length == 120);
}


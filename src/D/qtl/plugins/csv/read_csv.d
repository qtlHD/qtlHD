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
import std.algorithm;
import std.array;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.core.individual;
import example.genotype_examples;

/**
 * Read a simple CSV file containing marker names, chromosome nrs, position,
 * phenotype and genotype - such as the listeria.csv file used in R/qtl.
 *
 * The cross type is injected as XType. An XType works as long as it is a known
 * Genotype class (e.g. BC, F2, RIL, Flex). Observed genotypecombinator are handled as
 * 'fixed' in this implementation. That is, we already know the observed types
 * before reading the genotype data file. The symbols hands back matching
 * (observed) genotypecombinator.
 *
 * The file is parsed once on class instantiation. Elements can be queried.
 */

class ReadSimpleCSV(XType,CrossType) {
  private File f;
  XType crosstype;
  CrossType symbols;
  string[] phenotypenames;
  Marker[] markers;
  Individuals individuals;
  Chromosome[string] chromosomes;
  Phenotype!double[][] phenotypes;
  GenotypeCombinator[][] genotypecombinator;
  size_t n_phenotypes;

  this(in string fn, CrossType observed = null) {
    f = File(fn,"r");
    scope(exit) f.close();
    if (observed is null)
      symbols = new CrossType;
    else
      symbols = observed;
    writeln(symbols);
    // symbols = new CrossType; // this may become a parameter
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

    // create Marker objects; have m.id be column within phenotypes or within genotypecombinator
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
      genotypecombinator ~= gs;
      individuals.list ~= new Individual(n_individual);
    }
  }
}

mixin RealizePhenotypeMatrix!double;

Tuple!(Marker[],Inds,PhenotypeMatrix,ObservedGenotypes,GenotypeCombinator[][]) 
  load_csv(string fn) {
  PhenotypeMatrix p;
  // FIXME: we force an F2 here
  auto data = new ReadSimpleCSV!(F2,ObservedF2)(fn);
  // Convert individuals to string[]
  auto iter = new ListIter!Individuals(data.individuals);
  string[] inds = map!(ind => to!string(ind.name))(iter).array();
  // Turn the genotype matrix into a genotype combinator matrix
  ObservedGenotypes observed;
  // auto gc = convert_to_combinator_matrix(g,observed);

  return tuple(data.markers,inds,data.phenotypes,observed,data.genotypecombinator);
}



/**
 * Plugin for reading simple CSV files
 *
 *   rm read_csv ; dmd -unittest qtl/plugins/input/read_csv.d qtl/core/*.d ; time ./read_csv
 */

module qtl.plugins.input.read_csv;

import core.memory;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.matrices;
import qtl.core.phenotype;
import qtl.core.deprecate.genotype_enum;

import std.stdio;
import std.conv;
import std.string;
import std.path;

/**
 * Read a simple CSV file containing marker names, chromosome nrs, position,
 * phenotype and genotype - such as the listeria.csv file used in R/qtl.
 *
 * The cross type is injected as XType. An XType works as long as it is
 * a known Genotype class (BC, F2, RIL).
 *
 * The file is parsed once on class instantiation. Elements can be queried.
 */

class ReadSimpleCSV(XType) {
  private File f;
  string[] phenotypenames;
  Marker[] markers;
  Individuals individuals;
  Chromosome[string] chromosomes;
  Phenotype!double[][] phenotypes;
  Genotype!XType[][] genotypes;
  size_t n_phenotypes;

  double** getPhenotypesForMQM(){
    double** pheno = newmatrix!double(n_phenotypes,individuals.length);
    for(auto p=0;p<n_phenotypes;p++){
      for(auto i=0;i<individuals.length;i++){
        //writefln("%d %d %f",p,i,phenotypes[i][p].value);
        pheno[p][i] = phenotypes[i][p].value;
      }
    }
    return pheno;
  }

  char** getGenotypesForMQM(){
    char** geno = newmatrix!char(markers.length,individuals.length);
    for(int m=0;m<markers.length;m++){
      for(int i=0;i<individuals.length;i++){
        //writefln("%d %d %f",m,i,genotypes[i][m].value);
        geno[m][i] = 'A';
      }
    }
    return geno;
  }

  int* getChromosomesForMQM(){
    int* chromo = newvector!int(markers.length);
    for(int i=0;i<markers.length;i++){
      chromo[i] = markers[i].chromosome.id;
    }
    return chromo;
  }

  double* getDistancesForMQM(){
    double* dist = newvector!double(markers.length);
    for(int i=0;i<markers.length;i++){
      dist[i] = markers[i].position;
    }
    return dist;
  }

  this(in string fn) {
    f = File(fn,"r");
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
      Genotype!XType[] gs;
      gs.reserve(n_columns);  // pre-allocate memory (just good practise)
      foreach (field; fields[n_phenotypes..$]) {
        gs ~= set_genotype!XType(strip(field));
      }
      genotypes ~= gs;
      individuals.list ~= new Individual(n_individual);
    }
    f.close();
  }
}

unittest {
  writeln("Unit test " ~ __FILE__);
  alias std.path.join join;
  auto fn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","listeria.csv");
  writeln("  - reading CSV " ~ fn);
  Marker m2 = new Marker(4.8);
  auto markers = [ m2 ];
  auto data = new ReadSimpleCSV!F2(fn);
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
  assert(data.genotypes[1][0].value == F2.NA);
  assert(data.genotypes[1][1].value == F2.B);
  assert(data.individuals.length == 120);

  // foreach(name; data.chromosomes.keys.sort) {
  // }
}

unittest {
  alias std.path.join join;
  auto fn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","hyper.csv");
  writeln("  - reading CSV " ~ fn);
  auto data = new ReadSimpleCSV!F2(fn);
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
  assert(data.genotypes[1][0].value == F2.H);
  assert(data.genotypes[1][1].value == F2.H);
  assert(data.genotypes[2][3].value == F2.NA);
  assert(data.genotypes[2][4].value == F2.B);
}

unittest {
  alias std.path.join join;
  auto fn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","hyper_noX.csv");
  writeln("  - reading CSV " ~ fn);
  auto data = new ReadSimpleCSV!BC(fn);
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
  assert(data.genotypes[1][0].value == BC.H);
  assert(data.genotypes[1][1].value == BC.H);
  assert(data.genotypes[2][3].value == BC.NA);
  assert(data.genotypes[2][4].value == BC.A);
}


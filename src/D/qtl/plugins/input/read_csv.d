/**
 * Plugin for reading simple CSV files
 *
 *   rm read_csv ; dmd -unittest qtl/plugins/input/read_csv.d qtl/core/*.d ; time ./read_csv
 */

module qtl.plugins.input.read_csv;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;

import std.stdio;
import std.conv;
import std.string;
import std.path;

/** 
 * Read a simple CSV file containing marker names, chromosome nrs, position, 
 * phenotype and genotype - such as the listeria.csv file used in R/qtl.
 *
 * The file is parsed once on class instantiation. Elements can be queried.
 */

class ReadSimpleCSV {

  private File f;
  string[] phenotypenames;
  Marker[] markers;
  Chromosome[string] chromosomes;
  Phenotype!double[][] phenotypes;
  Genotype!F2[][] genotypes;      // FIXME: currently fixated for F2

  this(in string fn) {
    f = File(fn,"r");
    // read markers
    Marker[] ms;
    auto header = split(f.readln(),",");
    immutable size = header.length;
    ms.reserve(size);  // pre-allocate memory (just good practise)
    foreach (i, mname; header)
    {
       // writeln(mname);
       Marker m = { name:mname, id:i-1};
       ms ~= m;
    }

    // read chromosome info
    auto cnames = split(f.readln(),",");
    if (cnames.length != size) throw new Exception("# Chromosomes out of range");

    // find first non-blank (indicating last phenotype columns)
    size_t n_phenotypes;
    foreach (i, cname; cnames)
    {
      if(cname != "") {
        n_phenotypes = i;
        break;
      }
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
        chromosomes[cname2] = new Chromosome(cname2);
      markers[i].chromosome = chromosomes[cname2];
    }

    // read position info
    auto positions = split(f.readln(),",");
    if (positions.length != size) throw new Exception("Positions out of range");
    positions = positions[n_phenotypes..$];
    foreach (i, pos; positions)
    {
      immutable pos2 = strip(pos);
      markers[i].position = to!double(pos2);
    }
    // read remaining data
    string buf;
    uint individual = 0;
    while (f.readln(buf)) {
      auto fields = split(buf,",");
      if (fields.length != size) throw new Exception("Field # out of range in ", buf);
      Phenotype!double[] ps;
      ps.reserve(size);
      foreach(field; fields[0..n_phenotypes]) {
	ps ~= set_phenotype!double(field);
      }
      phenotypes ~= ps;
      // set genotype
      Genotype!F2[] gs;
      gs.reserve(size);  // pre-allocate memory (just good practise)
      foreach (field; fields[n_phenotypes..$]) {
        gs ~= set_genotype!F2(strip(field));
      }
      genotypes ~= gs;
      individual++;
    }
    f.close();
  }
}

unittest {
  writeln("Unit test " ~ __FILE__);
  alias std.path.join join;
  auto fn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","listeria.csv");
  writeln("  - reading CSV " ~ fn);
  Marker m2 = { id:2, position:4.8};
  assert(m2.id == 2);
  auto markers = [ m2 ];
  auto data = new ReadSimpleCSV(fn);
  assert(data.markers.length == 133, to!string(data.markers.length));
  assert(data.phenotypenames[0] == "T264");
  assert(data.markers[0].name == "D10M44");
  assert(data.markers[0].id == 0);
  assert(data.markers[1].id == 1);
  // Check chromosomes
  assert(data.chromosomes.length == 20, to!string(data.chromosomes.length));
  assert(data.chromosomes["X"].id == 0);
  assert(data.chromosomes["7"].id == 7);
  assert(data.markers[2].position == 24.84773, "Marker position not matching");
  // Check phenotype
  assert(data.phenotypes[29][0].value == PHENOTYPE_NA, to!string(data.phenotypes[29][0].value));
  assert(data.phenotypes[30][0].value == 74.417);
  // Check genotype
  assert(data.genotypes[1][0].value == F2.NA);
  assert(data.genotypes[1][1].value == F2.B);

  // foreach(name; data.chromosomes.keys.sort) {
  // }
}

unittest {
  writeln("Unit test " ~ __FILE__);
  alias std.path.join join;
  auto fn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","hyper.csv");
  writeln("  - reading CSV " ~ fn);
  auto data = new ReadSimpleCSV(fn);
  assert(data.markers.length == 174, to!string(data.markers.length));
  assert(data.phenotypenames == ["bp", "sex"]);
  assert(data.markers[3].name == "D1Mit178");
  assert(data.markers[3].id == 4);
  assert(data.markers[4].id == 5);
  // Check chromosomes
  assert(data.chromosomes.length == 20, to!string(data.chromosomes.length));
  assert(data.chromosomes["X"].id == 0);
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

void main() { }

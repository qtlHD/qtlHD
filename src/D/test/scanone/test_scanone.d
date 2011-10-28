/**
 * Test scanone routines, using listeria (CSV) set
 */

module test.scanone.test_scanone;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.deprecate.genotype_enum;
import qtl.core.marker;
import qtl.core.map;
import qtl.core.make_map;
import qtl.plugins.deprecate.read_csv;
import qtl.core.scanone_hk;
import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.algorithm;
alias std.algorithm.find find;


import qtl.core.genetic_map_functions;
import qtl.core.hmm_f2;

// The comments are based on the Ruby/Biolib-R/qtl integration

static bool VERBOSE = false;

unittest {

  writeln("Unit test " ~ __FILE__);
  alias std.path.buildPath buildPath;
  auto fn = to!string(dirName(__FILE__) ~ sep ~ buildPath("..","..","..","..","test","data","input","listeria.csv"));
  writeln("  - reading CSV " ~ fn);
  Marker m2 = new Marker(4.8);
  auto markers = [ m2 ];
  auto data = new ReadSimpleCSV!F2(fn);
  assert(data.markers.length == 133, to!string(data.markers.length));
  assert(data.phenotypenames[0] == "T264");
  assert(data.markers[0].name == "D10M44");
  assert(data.markers[0].id == 0);
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

  // Markers per chromosome
  /*
    => [["1", 13], ["10", 5], ["11", 6], ["12", 6], ["13", 12], ["14", 4], ["15", 8], ["16", 4], ["17", 4], ["18", 4], ["19", 4], ["2", 6], ["3", 6], ["4", 4], ["5", 13], ["6", 13], ["7", 6], ["8", 6], ["9", 7], ["X", 2]]
  */
  auto c_mslist = get_markers_by_chromosome(data.markers);
  foreach(c_ms ; c_mslist) {
    auto c = c_ms[0];  // first part of Tuple
    auto ms = c_ms[1]; // second part of Tuple
    // look for X chromosome with 2 markers
    if (c.name == "X") {
      assert(ms.length == 2, to!string(ms.length));
      assert(ms[0].position == 0, to!string(ms[0].position));
      assert(ms[0].toString == "DXM186~0", ms[0].toString);
      // writeln(ms[1]);
      assert(ms[1].name == "DXM64", ms[1].name);
      assert(ms[1].toString == "DXM64~42.3459", ms[1].toString);
    }
    if (VERBOSE) {
      writeln(c.name);
      foreach (m; ms) {
        writeln(m.name,'-',m.position);
      }
    }
  }
  // Find marker by name
  Marker m = find!("a.name == \"D10M44\"")(data.markers)[0];
  assert(m.name == "D10M44");
  assert(m.id == 0);
  // Find marker by id
  Marker m1 = find!("a.id == 1")(data.markers)[0];
  assert(m1.id == 1);
  assert(m1.position == 0.99675);
  Marker m14 = find!("a.id == 14")(data.markers)[0];
  assert(m14.chromosome.name == "2");

  // Calculate map size X Chromosome
  auto cx = find!("a[0].name == \"X\"")(c_mslist)[0];
  auto msx = cx[1];
  assert(msx.length == 2);
  assert(to!string(map_size(msx))=="42.3459",to!string(map_size(msx)));
  // Calculate map size all Chromosomes
  double totalsize = 0.0;
  foreach(c_ms ; c_mslist) {
    totalsize += map_size(c_ms[1]);
  }
  assert(to!string(totalsize) == "1104.29");

  // Chromosome 1
  auto c1 = find!("a[0].name == \"1\"")(c_mslist)[0];
  auto ms1 = c1[1].sort;
  writeln(ms1);
  writeln(" -- Chromosome 1");
  auto rfs = recombination_fractions(ms1);
  assert(rfs.length==12);
  assert(to!string(rfs[0])[0..9]=="0.0098688");
  assert(to!string(rfs[1])=="0.189685");
  // Chromosome X
  writeln(" -- Chromosome X");
  writeln("    msx.length: ", msx.length);
  writeln("    ", recombination_fractions(msx));
  assert(to!string(recombination_fractions(msx)[0])=="0.285633");
  // expand map for each chromosome
  double totalexpsize = 0.0;
  foreach(c_ms ; c_mslist) {
    auto msin = new Markers!Marker(c_ms[1]);
    auto ms = add_stepped_markers_autosome(msin,2.5,0);
    totalexpsize += map_size(ms.list);
  }
  // after expansion the size is still:
  assert(to!string(totalexpsize) == "1104.29",to!string(totalexpsize));
 
  // test rfs after expansion
  writeln(" -- Chr 1 with pseudomarkers");
  auto msin1 = new Markers!Marker(ms1);
  auto expanded_ms1 = add_stepped_markers_autosome(msin1,2.5,0);
  auto expanded_recombfs1 = recombination_fractions(expanded_ms1.list);
  assert(expanded_recombfs1.length==49,to!string(expanded_recombfs1.length));
  assert(to!string(expanded_recombfs1[0])[0..9]=="0.0098688");
  assert(to!string(expanded_recombfs1[2])=="0.0243853");
  // for chromosome 1 (c1, ms1, rfs)
  // calc genotype probabilities (using data.genotypes)
  // here using map.d's Haldane - which should merge with hmm_f2 etc.
  // auto rec_frac = mapFunction(dist_cM);
  // auto rec_frac = expanded_recombfs1;
  Marker[] markers_on_chr_4;
  writeln("    Grab markers");
  foreach(marker; data.markers) {
    if(marker.chromosome.name=="4") {
      markers_on_chr_4 ~= marker;
    }
  }

  writeln("      - Subset genotype data");
  Genotype!F2[][] chr_4_genotypes;
  chr_4_genotypes.reserve(data.genotypes.length);
  foreach(i; 0..data.genotypes.length) {
    Genotype!F2[] an_individuals_genotype;
    an_individuals_genotype.reserve(markers_on_chr_4.length);
    foreach(j; 0..markers_on_chr_4.length) {
      an_individuals_genotype ~= data.genotypes[i][markers_on_chr_4[j].id];
    }
    chr_4_genotypes ~= an_individuals_genotype;
  }

  writeln("      - Get recombination fractions");
  double[] dist_cM;
  foreach(i; 1..markers_on_chr_4.length) {
    dist_cM ~= markers_on_chr_4[i].position - markers_on_chr_4[i-1].position;
  }
  auto rec_frac = dist_to_recfrac(dist_cM, GeneticMapFunc.Haldane);

  writeln("      - Run calcGenoprob for F2");
  // auto genoprobs = calc_geno_prob(chr_4_genotypes, rec_frac, 0.002);

  // writeln("rf: ",rfs);
  // auto genoprobs = calcGenoprob(data.genotypes, rfs, 0.002);
  // GenoProbs gprobs;
  // Getting ready for scanone
  // auto result = scanone_hk(data.individuals,data.phenotypes,markers_on_chr_4,genoprobs); 
  /*
We are going to scan for QTL's. The first R equivalent here is:

    > mr = scanone(listeria,method='mr')
      Warning message:
      Dropping 4 individuals with missing phenotypes.
        in: checkcovar(cross, pheno.col, addcovar, intcovar, perm.strata,  

    > mr
              chr      pos         lod
      D10M44    1  0.00000 0.457256443
      D1M3      1  0.99675 0.687377692
      D1M75     1 24.84773 0.244500303
      D1M215    1 40.41361 0.070120158
      (...)
      D19M117  19 16.36398 0.436675089
      D19M65   19 32.82935 0.007363598
      D19M10   19 44.49432 0.000000000
      DXM186    X  0.00000 0.678435461
      DXM64     X 42.34593 0.002756074

Execute single QTL mapping using R/qtl's marker regression passing in
the @qtl dataset. We need the RQTL convenience class to map against biolib:

    >> require 'qtl/rqtl'
    >> rqtl = RQTL.new(qtl)

Now execute single QTL mapping with scanone using multiple regression analysis

    >> mr = rqtl.scanone_mr()

Which returns and array of markers containing the marker name, chromosome,
position and lod score:

    >> mr[0].lod
    => 0.457256442704442
    >> mr[0].name
    => 'D10M44'
    >> mr[0].chromosome
    => '1'
    >> mr[0].position
    => 0
    >> mr[0].to_a
    => ["D10M44", "1", 0, 0.457256442704442]

    >> mr[3].lod
    => 0.0701201576892285
    >> mr[3].name
    => 'D1M215'
    >> mr[3].chromosome
    => '1'
    >> mr[3].position
    => 40.41361 

*/

}

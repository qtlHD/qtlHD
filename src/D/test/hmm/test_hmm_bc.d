/**
 * test reading data and then apply HMM
 * (backcross example: hyper with no X chr)
 */

module test.hmm.test_hmm_bc;

import std.math, std.stdio, std.path;
import std.exception, std.conv;
import std.random;
import std.algorithm;

import qtl.core.primitives;
import qtl.core.phenotype, qtl.core.chromosome, qtl.core.genotype;
import qtl.plugins.qtab.read_qtab;
import qtl.core.map.make_map;
import qtl.core.map.genetic_map_functions;
import qtl.core.hmm.hmm_bc;

unittest {
  writeln("Unit test " ~ __FILE__);
}

unittest {
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ sep ~
                       buildPath("..","..","..","..","test","data", "input", "hyper_noX_qtab"));

  // Read founder info
  auto founder_fn = to!string(buildPath(dir, "hyper_noX_founder.qtab"));
  writeln("reading ", founder_fn);
  auto info = read_founder_settings_qtab(founder_fn);
  assert(info["Cross"] == "BC");


  // First read symbol information (the GenotypeCombinators)
  auto symbol_fn = to!string(buildPath(dir,"hyper_noX_symbol.qtab"));
  writeln("reading ",symbol_fn);
  auto fs = File(symbol_fn,"r");
  auto symbols = read_genotype_symbol_qtab(fs);

  // Test working of symbols
  assert(to!string(symbols.decode("NA")) == "[NA]");
  assert(to!string(symbols.decode("A")) == "[(0,0)]");
  assert(to!string(symbols.decode("H")) == "[(1,0)]");

  // Read genotype matrix
  auto genotype_fn = to!string(buildPath(dir,"hyper_noX_genotype.qtab"));
  writeln("reading ",genotype_fn);
  auto fg = File(genotype_fn,"r");
  auto ret = read_genotype_qtab(fg, symbols);
  auto individuals = ret[0];
  auto genotype_matrix = ret[1];

  // Show the first individual and genotypes
  assert(individuals.list.length == 250);
  assert(individuals.list[244].name == "245");
  assert(genotype_matrix[244].length == 170);

  // by symbol
  assert(genotype_matrix[239][0] == symbols.decode("NA"));
  assert(genotype_matrix[239][2] == symbols.decode("A"));
  assert(genotype_matrix[239][37] == symbols.decode("H"));

  assert(genotype_matrix[244][0] == symbols.decode("NA"));
  assert(genotype_matrix[244][2] == symbols.decode("H"));
  assert(genotype_matrix[244][169] == symbols.decode("NA"));

  writeln("genotype_matrix[244][0]: ", genotype_matrix[244][0]);
  writeln("genotype_matrix[244][0].length: ", genotype_matrix[244][0].length);
  writeln("genotype_matrix[244][0].list.length: ", genotype_matrix[244][0].list.length);
  foreach(g; genotype_matrix[244][0].list)
    writeln("genotype_matrix[244][0].list[i]: ", g, " ", typeid(g));

  assert(genotype_matrix[19][0] == symbols.decode("H"));
  assert(genotype_matrix[19][2] == symbols.decode("A"));
  assert(genotype_matrix[19][3] == symbols.decode("NA"));
  assert(genotype_matrix[19][169] == symbols.decode("A"));

  writeln("genotype_matrix[19][2]: ", genotype_matrix[19][2]);
  writeln("genotype_matrix[19][2].length: ", genotype_matrix[19][2].length);
  writeln("genotype_matrix[19][2].list.length: ", genotype_matrix[19][2].list.length);
  foreach(g; genotype_matrix[19][2].list)
    writeln("genotype_matrix[19][2].list[i]: ", g, " ", typeid(g));

  writeln("genotype_matrix[19][0]: ", genotype_matrix[19][0]);
  writeln("genotype_matrix[19][0].length: ", genotype_matrix[19][0].length);
  writeln("genotype_matrix[19][0].list.length: ", genotype_matrix[19][0].list.length);
  foreach(g; genotype_matrix[19][0].list)
    writeln("genotype_matrix[19][0].list[i]: ", g, " ", typeid(g));

  // by founders
  assert(genotype_matrix[239][2].list[0].homozygous == true);
  assert(genotype_matrix[239][2].list[0].heterozygous == false);
  assert(genotype_matrix[239][2].list[0].founders[0] == 0);
  assert(genotype_matrix[239][2].list[0].founders[1] == 0);

  assert(genotype_matrix[244][2].list[0].homozygous == false);
  assert(genotype_matrix[244][2].list[0].heterozygous == true);
  assert(genotype_matrix[244][2].list[0].founders[0] == 1);
  assert(genotype_matrix[244][2].list[0].founders[1] == 0);

  // reading phenotypes
  auto pheno_fn = to!string(buildPath(dir,"hyper_noX_phenotype.qtab"));
  writeln("reading ",pheno_fn);
  auto p_res = read_phenotype_qtab!(Phenotype!double)(pheno_fn);
  Phenotype!double[][] pheno = p_res[0];

  assert(pheno.length == 250);
  foreach(p; pheno) assert(p.length == 2);
  // 1st ind, 1st phenotype
  assert(pheno[0][0].value == 109.6);
  // 3rd ind, 1st phenotype
  assert(pheno[2][0].value == 110.1);
  // 133rd ind
  assert(pheno[132][0].value == 97);
  // all have 2nd phenotype == 1
  foreach(ph; pheno)
    assert(ph[1].value == 1);

  // Marker map reader
  auto marker_map_fn = to!string(buildPath(dir,"hyper_noX_marker_map.qtab"));
  writeln("reading ",marker_map_fn);
  auto markers = read_marker_map_qtab!(Marker)(marker_map_fn);
  assert(markers.length == 170);
  assert(markers[0].name == "D1Mit296");
  assert(markers[0].chromosome.name == "1");
  assert(markers[0].get_position == 3.3);
  assert(markers[3].name == "D1Mit178");
  assert(markers[3].chromosome.name == "1");
  assert(markers[3].get_position == 35);
  assert(markers[164].name == "D18Mit50");
  assert(markers[164].chromosome.name == "18");
  assert(markers[164].get_position == 26.2);
  assert(markers[169].name == "D19Mit137");
  assert(markers[169].chromosome.name == "19");
  assert(markers[169].get_position == 55.7);

  // marker id == numeric index
  foreach(i, m; markers) assert(m.id == i);

  // test splitting up of markers into chromosomes
  auto markers_by_chr = get_markers_by_chromosome(markers);
  writeln("markers_by_chr.length: ", markers_by_chr.length);
  writeln("markers_by_chr[0].length: ", markers_by_chr[0].length);
  writeln("markers_by_chr[0][0].name: ", markers_by_chr[0][0].name);
  writeln("markers_by_chr[0][1][0].name: ", markers_by_chr[0][1][0].name);
  foreach(chr; markers_by_chr) {
    writefln("%2s (%2d): %-9s (%3d)", chr[0].name, chr[1].length, chr[1][0].name, chr[1][0].id);
    // check that markers within chromosome are in order:
    //    contiguous ids; non-decreasing position
    assert(chr[1][0].chromosome.name == chr[0].name);
    for(auto i=1; i<chr[1].length; i++) {
      assert(chr[1][i].id == chr[1][i-1].id+1);
      assert(chr[1][i].get_position >= chr[1][i-1].get_position);
      assert(chr[1][i].chromosome.name == chr[0].name);
    }
  }

  // sorted chromosomes
  writeln("sorted chromosomes: ");
  auto markers_by_chr_sorted = sort_chromosomes_by_marker_id(markers_by_chr);
  foreach(chr; markers_by_chr_sorted) {
    writef("%2s (%2d): %-9s (%3d)", chr[0].name, chr[1].length, chr[1][0].name, chr[1][0].id);
    // check that markers within chromosome are in order:
    //    contiguous ids; non-decreasing position
    assert(chr[1][0].chromosome.name == chr[0].name);
    for(auto i=1; i<chr[1].length; i++) {
      assert(chr[1][i].id == chr[1][i-1].id+1);
      assert(chr[1][i].get_position >= chr[1][i-1].get_position);
      assert(chr[1][i].chromosome.name == chr[0].name);
    }

    auto pmap_stepped = add_stepped_markers_autosome(chr[1], 1, 0);
    auto pmap_minimal = add_minimal_markers_autosome(chr[1], 1, 0);
    writefln("\tmarkers: %3d\tpmar (stepped): %3d\tpmar (minimal): %3d", chr[1].length,
             pmap_stepped.length-chr[1].length, pmap_minimal.length-chr[1].length);
  }

  // test calc_geno_prob with hyper data
  writeln("Test calc_geno_prob with hyper data, chr 4");
  auto chr4_map = markers_by_chr_sorted[4][1];
  sort(chr4_map); // sort in place
  auto pmap_stepped_chr4 = add_stepped_markers_autosome(chr4_map, 1.0, 0.0);

  double[] dist_cM;
  foreach(i; 1..pmap_stepped_chr4.length)
    dist_cM ~= pmap_stepped_chr4[i].get_position() - pmap_stepped_chr4[i-1].get_position();
  auto rec_frac = dist_to_recfrac(dist_cM, GeneticMapFunc.Carter_Falconer);
  writeln("dist_cM.length: ", dist_cM.length);
  writeln("rec_frac.length: ", rec_frac.length);
  writeln("pmap_stepped_chr4.length: ", pmap_stepped_chr4.length);

  auto chr4probs = calc_geno_prob_BC(genotype_matrix, pmap_stepped_chr4, rec_frac, 0.001);
}

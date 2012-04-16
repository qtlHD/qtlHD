/**
 * Test scanone routines
 */

module test.scanone.test_scanone;

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.algorithm;
alias std.algorithm.find find;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.marker;
import qtl.core.genotype;
import qtl.core.phenotype;
import qtl.plugins.qtab.read_qtab;
import qtl.core.map.map;
import qtl.core.map.make_map;
import qtl.core.scanone.scanone_hk;

import qtl.core.map.genetic_map_functions;
import qtl.core.hmm.bc;
import qtl.core.scanone.scanone_hk;

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

  // Read genotype matrix
  auto genotype_fn = to!string(buildPath(dir,"hyper_noX_genotype.qtab"));
  writeln("reading ",genotype_fn);
  auto fg = File(genotype_fn,"r");
  auto ret = read_genotype_qtab(fg, symbols);
  auto individuals = ret[0];
  auto genotype_matrix = ret[1];

  // reading phenotypes
  auto pheno_fn = to!string(buildPath(dir,"hyper_noX_phenotype.qtab"));
  writeln("reading ",pheno_fn);
  auto p_res = read_phenotype_qtab!(Phenotype!double)(pheno_fn);
  Phenotype!double[][] pheno = p_res[0];

  // Marker map reader
  auto marker_map_fn = to!string(buildPath(dir,"hyper_noX_marker_map.qtab"));
  writeln("reading ",marker_map_fn);
  auto markers = read_marker_map_qtab!(Marker)(marker_map_fn);

  // split up markers into chromosomes
  auto markers_by_chr = get_markers_by_chromosome(markers);

  // sorted chromosomes
  auto markers_by_chr_sorted = sort_chromosomes_by_marker_id(markers_by_chr);

  // chr 4 with pseudomarkers
  auto chr4_map = markers_by_chr_sorted[3][1];
  sort(chr4_map); // sort in place
  auto pmap_stepped_chr4 = add_stepped_markers_autosome(chr4_map, 1.0, 0.0);

  // calc_geno_prob
  auto rec_frac = recombination_fractions(pmap_stepped_chr4, GeneticMapFunc.Carter_Falconer);
  auto chr4probs = calc_geno_prob_BC(genotype_matrix, pmap_stepped_chr4, rec_frac, 0.001);

  // empty covariate matrices
  auto addcovar = new double[][](genotype_matrix.length, 0);
  auto intcovar = new double[][](genotype_matrix.length, 0);
  auto weights = new double[](0);

  // run scanone
  auto rss = scanone_hk(chr4probs, pheno, addcovar, intcovar, weights);

  foreach(i, pos; pmap_stepped_chr4) {
    writef("%-2s %5.1f  ", pos.chromosome.name, pos.position);
    foreach(j; 0..rss[i].length)
      writef("%8.4f ", rss[i][j]);
    writeln;
  }
}

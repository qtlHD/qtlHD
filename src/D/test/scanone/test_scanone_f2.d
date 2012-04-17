/**
 * Test scanone routines
 */

module test.scanone.test_scanone_f2;

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.algorithm;
import std.math;
import std.algorithm;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.marker;
import qtl.core.genotype;
import qtl.core.phenotype;
import qtl.plugins.qtab.read_qtab;
import qtl.core.map.map;
import qtl.core.map.make_map;

import qtl.core.map.genetic_map_functions;
import qtl.core.hmm.f2;
import qtl.core.scanone.scanone_hk;
import qtl.core.util.data_manip;


unittest {
  writeln("Unit test " ~ __FILE__);
}

unittest {
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ sep ~
                       buildPath("..","..","..","..","test","data", "input", "listeria_qtab"));

  // load founder info
  auto founder_fn = to!string(buildPath(dir, "listeria_founder.qtab"));
  writeln("reading ", founder_fn);
  auto info = read_founder_settings_qtab(founder_fn);

  // load symbols
  auto symbol_fn = to!string(buildPath(dir,"listeria_symbol.qtab"));
  // First read symbol information (the GenotypeCombinators)
  writeln("reading ",symbol_fn);
  auto fs = File(symbol_fn,"r");
  auto symbols = read_genotype_symbol_qtab(fs);

  // Read genotype matrix
  auto genotype_fn = to!string(buildPath(dir,"listeria_genotype.qtab"));
  writeln("reading ",genotype_fn);
  auto fg = File(genotype_fn,"r");
  auto ret = read_genotype_qtab(fg, symbols);
  auto individuals = ret[0];
  auto genotype_matrix = ret[1];

  // reading phenotypes
  auto pheno_fn = to!string(buildPath(dir,"listeria_phenotype.qtab"));
  writeln("reading ",pheno_fn);
  auto p_res = read_phenotype_qtab!(Phenotype!double)(pheno_fn);
  Phenotype!double[][] pheno = p_res[0];

  auto ind_to_omit = is_any_phenotype_missing(pheno);
  auto n_to_omit = count(ind_to_omit, true);
  writeln("Omitting ", n_to_omit, " individuals with missing phenotype");

  // Marker map reader
  auto marker_map_fn = to!string(buildPath(dir,"listeria_marker_map.qtab"));
  writeln("reading ",marker_map_fn);
  auto markers = read_marker_map_qtab!(Marker)(marker_map_fn);

  // split up markers into chromosomes
  auto markers_by_chr = get_markers_by_chromosome(markers);

  // sorted chromosomes
  auto markers_by_chr_sorted = sort_chromosomes_by_marker_id(markers_by_chr);

}

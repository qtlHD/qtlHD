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

  genotype_matrix = omit_ind_from_genotypes(genotype_matrix, ind_to_omit);
  pheno = omit_ind_from_phenotypes(pheno, ind_to_omit);

  // Marker map reader
  auto marker_map_fn = to!string(buildPath(dir,"listeria_marker_map.qtab"));
  writeln("reading ",marker_map_fn);
  auto markers = read_marker_map_qtab!(Marker)(marker_map_fn);

  // split up markers into chromosomes
  auto markers_by_chr = get_markers_by_chromosome(markers);

  // sorted chromosomes
  auto markers_by_chr_sorted = sort_chromosomes_by_marker_id(markers_by_chr);

  // chr 5
  auto chr5_map = markers_by_chr_sorted[4][1];
  sort(chr5_map); // sort in place

  // add pseudomarkers
  auto chr5_map_wpmark = add_minimal_markers_autosome(chr5_map, 2);

  // calc_genoprob
  auto rec_frac = recombination_fractions(chr5_map_wpmark, GeneticMapFunc.Kosambi);
  auto chr5probs = calc_geno_prob_F2(genotype_matrix, chr5_map_wpmark, rec_frac, 0.01);

  // empty covariate matrices
  auto addcovar = new double[][](0, 0);
  auto intcovar = new double[][](0, 0);
  auto weights = new double[](0);

  // run scanone and calculate LOD scores
  auto rss = scanone_hk(chr5probs, pheno, addcovar, intcovar, weights);
  auto rss0 = scanone_hk_null(pheno, addcovar, weights);
  auto lod = rss_to_lod(rss, rss0, pheno.length);

  /*******************************
   * in R:

   data(listeria)
   listeria <- listeria["-X",]
   listeria <- calc.genoprob(listeria, step=2, stepwidth="max", err=0.01, map="kosambi")
   out <- scanone(listeria, method="hk", chr=5)

   *
   ******************************/

    auto Rlod = [1.75082264702388, 1.90066849718966, 2.01704712283993, 2.08712140364622, 2.10492717358017, 2.70630273019253,
                 3.37761230506607, 4.05942123556486, 4.67519573823500, 5.16235351860865, 5.49575699915243, 5.68655510311692,
                 6.03045154338031, 6.13583408757336, 5.99580084184214, 5.64987645100132, 6.37031242544765, 6.67654669747520,
                 6.55134165341997, 6.06133769955562, 6.06133826161380, 6.03243273069080, 5.84745187330492, 5.78282966210315,
                 5.41716088137809, 4.81528999604944, 4.52366494344466, 4.00127126667383, 3.35266038068522, 3.40387254237714,
                 3.38560912296015, 3.29543897580413, 3.14207936272970, 3.07376784424855, 2.99184546370128, 2.90164906063137,
                 2.80647516232494, 2.70834876275660, 2.60898335637313];

  assert(lod.length == Rlod.length);
  foreach(i; 0..lod.length)
    assert(abs(lod[i][0] - Rlod[i]) < 1e-5);
}

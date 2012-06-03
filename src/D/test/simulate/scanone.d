/* *
 * test simulation with scanone
 */

module test.simulate.scanone;

import qtl.core.simulate.crosses;
import math.distributions.random;
import qtl.core.simulate.meiosis;
import qtl.core.simulate.map;
import qtl.core.genotype;
import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.genotype;
import qtl.core.map.make_map;
import qtl.core.map.genetic_map_functions;
import qtl.core.map.map;
import qtl.core.hmm.bc;
import qtl.core.scanone.scanone_hk;

import std.algorithm;
import std.stdio;
import std.string;
import std.conv;
import std.random;
import std.math;

unittest {
  writeln("Unit test " ~ __FILE__);
  writeln("test backcross");

  Random gen;
  gen.seed(unpredictableSeed);

  auto chromosome = get_chromosome_with_id("1");
  auto marker_map = generate_map_eqspacing_onechr(100.0, 11, 1, chromosome);

  size_t n_ind = 100;
  auto BCgen = simulate_backcross_autosome(marker_map, n_ind, [0, 1], 10, 0.0, gen);

  auto phenotype = new Phenotype!double[][](n_ind,1);

  auto het = new TrueGenotype(0,0);

  foreach(i; 0..n_ind) {
    phenotype[i][0].value = rnorm(0.0, 1.0, gen);
    if(BCgen[i][5].match(het))
      phenotype[i][0].value += 0.8;
  }

  // pseudomarkers and HMM
  auto pmap = add_stepped_markers_autosome(marker_map, 1.0, 0.0);
  auto rec_frac = recombination_fractions(pmap, GeneticMapFunc.Carter_Falconer);
  auto genoprobs = calc_geno_prob_BC(BCgen, pmap, rec_frac, 0.001);

  // empty covariate matrices
  auto addcovar = new double[][](0, 0);
  auto intcovar = new double[][](0, 0);
  auto weights = new double[](0);

  // scanone
  writeln("phenotype.length: ", phenotype.length);
  writeln("phenotype[0].length: ", phenotype[0].length);
  writeln("genoprobs.length: ", genoprobs.length);
  writeln("genoprobs[0].length: ", genoprobs[0].length);
  auto rss = scanone_hk(genoprobs, phenotype, addcovar, intcovar, weights);
  auto rss0 = scanone_hk_null(phenotype, addcovar, weights);
  auto lod = rss_to_lod(rss, rss0, phenotype.length);

  foreach(i, m; pmap)
    writefln("%20s %5.1f %8.3f", m.name, m.get_position, lod[i][0]);
}
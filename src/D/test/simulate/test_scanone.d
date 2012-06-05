/* *
 * test simulation with scanone
 */

module test.simulate.test_scanone;

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
import qtl.core.hmm.f2;
import qtl.core.scanone.scanone_hk;
import qtl.core.scanone.peaks;

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

  auto phenotype = new Phenotype!double[](n_ind);

  auto het = new TrueGenotype(1,0);

  foreach(i; 0..n_ind) {
    phenotype[i].value = rnorm(0.0, 1.0, gen);
    if(BCgen[i][5].match(het))
      phenotype[i].value += 0.8;
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
  auto rss = scanone_hk(genoprobs, phenotype, addcovar, intcovar, weights);
  auto rss0 = scanone_hk_null(phenotype, addcovar, weights);
  auto lod = rss_to_lod(rss, rss0, phenotype.length);

  auto peak = get_peak_scanone(lod, pmap);
  writefln("Max lod = %6.2f at %7.2f", peak[0], peak[1].get_position);
}

unittest {
  writeln("test cross");

  Random gen;
  gen.seed(unpredictableSeed);

  auto chromosome = get_chromosome_with_id("1");
  auto marker_map = generate_map_eqspacing_onechr(100.0, 11, 1, chromosome);

  size_t n_ind = 100;
  auto F2gen = simulate_intercross_autosome(marker_map, n_ind, [0, 1], 10, 0.0, gen);

  auto phenotype = new Phenotype!double[](n_ind);

  auto het = new TrueGenotype(1,0);
  auto homB = new TrueGenotype(1,1);

  foreach(i; 0..n_ind) {
    phenotype[i].value = rnorm(0.0, 1.0, gen);
    if(F2gen[i][5].match(het))
      phenotype[i].value += 0.5;
    if(F2gen[i][5].match(homB))
      phenotype[i].value += 1.0;
  }

  // pseudomarkers and HMM
  auto pmap = add_stepped_markers_autosome(marker_map, 1.0, 0.0);
  auto rec_frac = recombination_fractions(pmap, GeneticMapFunc.Carter_Falconer);
  auto genoprobs = calc_geno_prob_F2(F2gen, pmap, rec_frac, 0.001);

  // empty covariate matrices
  auto addcovar = new double[][](0, 0);
  auto intcovar = new double[][](0, 0);
  auto weights = new double[](0);

  // scanone
  auto rss = scanone_hk(genoprobs, phenotype, addcovar, intcovar, weights);
  auto rss0 = scanone_hk_null(phenotype, addcovar, weights);
  auto lod = rss_to_lod(rss, rss0, phenotype.length);

  // maximum LOD score
  auto peak = get_peak_scanone(lod, pmap);
  writefln("Max lod = %6.2f at %7.2f", peak[0], peak[1].get_position);
}

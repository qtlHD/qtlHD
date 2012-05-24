/**
 * Simulate simple crosses
 */

module qtl.core.simulate.crosses;

import qtl.core.simulate.meiosis;
import qtl.core.simulate.map;
import qtl.core.genotype;
import qtl.core.primitives;
import qtl.core.chromosome;
import math.distributions.random;
import std.algorithm;
import std.stdio;
import std.string;
import std.conv;
import std.random;
import std.math;

// simulate a backcross
GenotypeCombinator[][] simulate_backcross(in Marker[] marker_map, in size_t n_individuals, in FounderIndex[] founders,
                                          uint m, double p, ref Random gen)
in {
  assert(founders.length == 2, "founders must have length 2");
  assert(marker_map.length > 0, "marker_map must have length > 0");
  assert(n_individuals > 0, "n_individuals must be > 0");
  assert(p >= 0.0 && p <= 1.0, "p must be in [0,1]");
}
body {
  GenotypeCombinator[][] genotypes;
  genotypes.reserve(marker_map.length*n_individuals);

  // parents
  auto parent0chr = generate_founder_chromosome(marker_map, founders[0]);
  auto f1 = generate_f1_genotypes_onechr(marker_map, founders);

  foreach(i; 0..n_individuals) {
      auto meiotic_product = simulate_meiotic_product(marker_map, f1, m, p, gen);
      auto backcross_truegen = combine_meiotic_products(meiotic_product, parent0chr);
      genotypes ~= convert_truegenotype_to_genotypecombinator(backcross_truegen, false);
  }

  return(genotypes);
}

unittest {
  writeln("Unit test " ~ __FILE__);
  writeln("test backcross");

  Random gen;
  gen.seed(unpredictableSeed);

  auto chromosome = get_chromosome_with_id("1");
  auto marker_map = generate_map_eqspacing_onechr(100.0, 11, 1, chromosome);

  auto BCgen = simulate_backcross(marker_map, 100, [0, 1], 10, 0.0, gen);

  // print the 50th individual
  foreach(g; BCgen[49]) {
    write(g, " ");
  }
  writeln;
}

// simulate an intercross
GenotypeCombinator[][] simulate_intercross(in Marker[] marker_map, in size_t n_individuals, in FounderIndex[] founders,
                                          uint m, double p, ref Random gen)
in {
  assert(founders.length == 2, "founders must have length 2");
  assert(marker_map.length > 0, "marker_map must have length > 0");
  assert(n_individuals > 0, "n_individuals must be > 0");
  assert(p >= 0.0 && p <= 1.0, "p must be in [0,1]");
}
body {
  GenotypeCombinator[][] genotypes;
  genotypes.reserve(marker_map.length*n_individuals);

  // parents
  auto f1 = generate_f1_genotypes_onechr(marker_map, founders);

  foreach(i; 0..n_individuals) {
      auto egg = simulate_meiotic_product(marker_map, f1, m, p, gen);
      auto sperm = simulate_meiotic_product(marker_map, f1, m, p, gen);
      auto intercross_truegen = combine_meiotic_products(egg, sperm);
      genotypes ~= convert_truegenotype_to_genotypecombinator(intercross_truegen, true);
  }

  return(genotypes);
}

unittest {
  writeln("test intercross");

  Random gen;
  gen.seed(unpredictableSeed);

  auto chromosome = get_chromosome_with_id("1");
  auto marker_map = generate_map_eqspacing_onechr(100.0, 11, 1, chromosome);

  auto F2gen = simulate_intercross(marker_map, 100, [0, 1], 10, 0.0, gen);

  // print the 75th individual
  foreach(g; F2gen[74]) {
    write(g, " ");
  }
  writeln;
}

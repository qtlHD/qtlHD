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
GenotypeCombinator[][] simulate_backcross_autosome(in Marker[] marker_map, in size_t n_individuals, in FounderIndex[] founders,
                                                   in uint m, in double p, ref Random gen)
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

  auto BCgen = simulate_backcross_autosome(marker_map, 100, [0, 1], 10, 0.0, gen);

  // print the 50th individual
  foreach(g; BCgen[49]) {
    write(g, " ");
  }
  writeln;
}

// simulate an intercross
GenotypeCombinator[][] simulate_intercross_autosome(in Marker[] marker_map, in size_t n_individuals, in FounderIndex[] founders,
                                                    in uint m, in double p, ref Random gen)
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

  auto F2gen = simulate_intercross_autosome(marker_map, 100, [0, 1], 10, 0.0, gen);

  // print the 75th individual
  foreach(g; F2gen[74]) {
    write(g, " ");
  }
  writeln;
}

// simulate RIL by selfing, to complete inbreeding
GenotypeCombinator[][] simulate_riself(in Marker[] marker_map, in size_t n_individuals, in FounderIndex[] founders,
                                       in uint m, in double p, ref Random gen)
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
    auto fk = f1;

    while(count_het_genotypes(fk)>0) { // check if inbred
      auto egg = simulate_meiotic_product(marker_map, fk, m, p, gen);
      auto sperm = simulate_meiotic_product(marker_map, fk, m, p, gen);
      fk = combine_meiotic_products(egg, sperm);
    }

    genotypes ~= convert_truegenotype_to_genotypecombinator(fk, false);
  }

  return(genotypes);
}

unittest {
  writeln("test RIL selfing");

  Random gen;
  gen.seed(unpredictableSeed);

  auto chromosome = get_chromosome_with_id("1");
  auto marker_map = generate_map_eqspacing_onechr(100.0, 11, 1, chromosome);

  auto RIself_geno = simulate_riself(marker_map, 100, [0, 1], 10, 0.0, gen);

  // print the 33rd individual
  foreach(g; RIself_geno[32]) {
    write(g, " ");
  }
  writeln;
}

uint count_het_genotypes(in TrueGenotype[] genotypes)
{
  uint number_het_genotypes=0;

  foreach(g; genotypes)
    if(g.heterozygous())
      number_het_genotypes++;

  return(number_het_genotypes);
}

// simulate RIL by sib-mating, to complete inbreeding
GenotypeCombinator[][] simulate_risib_autosome(in Marker[] marker_map, in size_t n_individuals, in FounderIndex[] founders,
                                               in uint m, in double p, ref Random gen)
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
    auto male = f1;
    auto female = f1;

    while(count_het_genotypes(male)>0 || count_het_genotypes(female) > 0 || // check if inbred
          count_genotype_differences(male, female)>0) {
      auto egg1 = simulate_meiotic_product(marker_map, female, m, p, gen);
      auto egg2 = simulate_meiotic_product(marker_map, female, m, p, gen);
      auto sperm1 = simulate_meiotic_product(marker_map, male, m, p, gen);
      auto sperm2 = simulate_meiotic_product(marker_map, male, m, p, gen);

      female = combine_meiotic_products(egg1, sperm1);
      male = combine_meiotic_products(egg2, sperm2);
    }

    genotypes ~= convert_truegenotype_to_genotypecombinator(female, false);
  }

  return(genotypes);
}

unittest {
  writeln("test RIL sibling mating");

  Random gen;
  gen.seed(unpredictableSeed);

  auto chromosome = get_chromosome_with_id("1");
  auto marker_map = generate_map_eqspacing_onechr(100.0, 11, 1, chromosome);

  auto RIsib_geno = simulate_risib_autosome(marker_map, 100, [0, 1], 10, 0.0, gen);

  // print the 33rd individual
  foreach(g; RIsib_geno[32]) {
    write(g, " ");
  }
  writeln;
}


uint count_genotype_differences(in TrueGenotype[] genotypes1, in TrueGenotype[] genotypes2)
in {
  assert(genotypes1.length == genotypes2.length, "genotypes1 and genotypes2 must have the same length");
}
body {
  uint number_different=0;

  foreach(i, g; genotypes1)
    if(g.get_allele(0) != genotypes2[i].get_allele(0) ||
       g.get_allele(1) != genotypes2[i].get_allele(1))
      number_different++;

  return(number_different);
}

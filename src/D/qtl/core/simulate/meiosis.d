/**
 * Code to simulate meiosis and basic crosses
 */

module qtl.core.simulate.meiosis;

import qtl.core.genotype;
import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.simulate.map;
import math.distributions.random;
import std.algorithm;
import std.stdio;
import std.string;
import std.conv;
import std.random;
import std.math;

// simulate crossover locations on a chromosome under no interference
// chrlen_cM = chromosome length in cM
double[] meiosisNI(in double chrlen_cM, ref Random gen)
in {
  assert(chrlen_cM > 0, "chrlen_cM must be > 0");
}
body {
  uint n_xo = rpois(chrlen_cM/100.0, gen);

  auto xo_locations = new double[](n_xo);
  foreach(i; 0..n_xo)
    xo_locations[i] = uniform(0.0, chrlen_cM, gen);

  sort(xo_locations);
  return(xo_locations);
}



/*
 Simulate crossover locations on a chromosome using Stahl model
   chrlen_cM = chromosome length in cM
   m = interference parameter
   p = proportion of crossovers coming from no interference model

 p=0 gives the chi-square model
 m=0 or p=1 gives the no-interference model
*/
double[] meiosis(in double chrlen_cM, in uint m, in double p, ref Random gen)
in {
  assert(chrlen_cM > 0, "chrlen_cM must be > 0");
  assert(p >= 0 && p <= 1, "p must be in [0,1]");
}
body {
  // no interference case is simple
  if(m==0 || p==1) return(meiosisNI(chrlen_cM, gen));

  // number of crossovers and intermediates (on 4-strand bundle)
  uint n_points_interf = rpois(chrlen_cM*(m+1)/50.0*(1.0-p), gen );
  // number of crossovers from no interference mechanism (on final product)
  uint n_points_ni;
  if(p==0) n_points_ni = 0;
  else n_points_ni = rpois(chrlen_cM/100.0*p, gen);

  double[] xo_locations;
  xo_locations.reserve(n_points_interf/(m+1) + n_points_ni);

  // locations of crossovers and intermediates
  auto pts_interf = new double[](n_points_interf);
  foreach(i; 0..n_points_interf)
    pts_interf[i] = uniform(0.0, chrlen_cM, gen);
  sort(pts_interf);

  // every (m+1)st point is chiasma; thin to give crossovers
  uint first = cast(uint)uniform(0.0, m+1, gen);
  for(auto i=first; i<n_points_interf; i += m+1)
    if(dice(gen, 1, 1)) xo_locations ~= pts_interf[i];

  // crossovers from no interference mechanism
  foreach(i; 0..n_points_ni)
    xo_locations ~= uniform(0.0, chrlen_cM, gen);

  sort(xo_locations);
  return(xo_locations);
}

// meiosis with start and end; basically meiosis(end-start)+start
double[] meiosis(in double start_cM, in double end_cM, in uint m, in double p, ref Random gen)
in {
  assert(end_cM > start_cM, "end_cM must be > start_cM");
  assert(p >= 0 && p <= 1, "p must be in [0,1]");
}
body {
  auto xo_locations = meiosis(end_cM - start_cM, m, p, gen);

  foreach(ref xoloc; xo_locations) {
    xoloc += start_cM;
  }
  return(xo_locations);
}

unittest {
  writeln("Unit test " ~ __FILE__);

  Random gen;
  gen.seed(unpredictableSeed);

  double[] product;

  writeln("Crossover locations in (0, 200):");
  foreach(m; 0..10) {
    product = meiosis(200.0, m, 0.05, gen);
    writeln("No. crossovers = ", product.length);
    foreach(p; product)
      writef("%6.3f ", p);
    writeln();
  }

  writeln("Crossover locations in (-200, 100):");
  foreach(m; 0..10) {
    product = meiosis(-200, 100, m, 0.0, gen);
    writeln("No. crossovers = ", product.length);
    foreach(p; product)
      writef("%6.3f ", p);
    writeln();
  }
}


// generate genotypes for an inbred line
TrueGenotype[] generate_founder_genotypes_onechr(in Marker[] marker_map, in FounderIndex founder)
{
  TrueGenotype[] genotypes;
  genotypes.reserve(marker_map.length);

  foreach(m; marker_map) {
    auto g = new TrueGenotype(founder, founder);
    genotypes ~= g;
  }

  return(genotypes);
}

// generate single chromosome of founder alleles
FounderIndex[] generate_founder_chromosome(in Marker[] marker_map, in FounderIndex founder)
{
  FounderIndex[] founders;
  founders.reserve(marker_map.length);

  foreach(m; marker_map) {
    founders ~= founder;
  }
  return(founders);
}

// generate genotypes for the F1 hybrid of two inbred lines
TrueGenotype[] generate_f1_genotypes_onechr(in Marker[] marker_map, in FounderIndex[] founders)
in {
  assert(founders.length == 2, "founders must have length 2");
}
body {
  TrueGenotype[] genotypes;
  genotypes.reserve(marker_map.length);

  foreach(m; marker_map) {
    auto g = new TrueGenotype(founders[0], founders[1]);
    genotypes ~= g;
  }

  return(genotypes);
}


// convert true genotypes to GenotypeSymbolMapper version
GenotypeSymbolMapper[] convert_truegenotype_to_genotypecombinator(TrueGenotype[] genotypes, in bool make_phase_unknown=true)
{
  GenotypeSymbolMapper[] obs_genotypes;
  obs_genotypes.reserve(genotypes.length);

  foreach(g; genotypes) {
    auto observed = new GenotypeSymbolMapper("observed");
    observed ~= g;
    if(make_phase_unknown && g.heterozygous()) {
      observed ~= g.reversed();
    }
    obs_genotypes ~= observed;
  }

  return(obs_genotypes);
}



unittest {
  auto chromosome = get_chromosome_with_id("1");
  auto marker_map = generate_map_eqspacing_onechr(100.0, 11, 1, chromosome);

  FounderIndex[] founders = [0, 1];

  auto inbred1 = generate_founder_genotypes_onechr(marker_map, founders[0]);
  auto inbred2 = generate_founder_genotypes_onechr(marker_map, founders[1]);
  auto inbredchr1 = generate_founder_chromosome(marker_map, founders[0]);
  auto inbredchr2 = generate_founder_chromosome(marker_map, founders[1]);
  auto f1 = generate_f1_genotypes_onechr(marker_map, founders);

  writeln("Inbred 1:");
  foreach(g; inbred1) {
    write(to!string(g) ~ " ");
  }
  writeln;

  writeln("Inbred 1 chromosome:");
  foreach(g; inbredchr1) {
    write(to!string(g) ~ " ");
  }
  writeln;

  writeln("Inbred 2:");
  foreach(g; inbred2) {
    write(to!string(g) ~ " ");
  }
  writeln;

  writeln("Inbred 2 chromosome:");
  foreach(g; inbredchr2) {
    write(to!string(g) ~ " ");
  }
  writeln;

  writeln("F1:");
  foreach(g; f1) {
    write(to!string(g) ~ " ");
  }
  writeln;

  // phase-unknown versions
  auto inbred1_puk = convert_truegenotype_to_genotypecombinator(inbred1);
  auto inbred2_puk = convert_truegenotype_to_genotypecombinator(inbred2);
  auto f1_puk = convert_truegenotype_to_genotypecombinator(f1);

  writeln("Inbred 1 (phase 'unknown'): ");
  foreach(g; inbred1_puk) {
    write(to!string(g) ~ " ");
  }
  writeln;

  writeln("Inbred 2 (phase 'unknown'): ");
  foreach(g; inbred2_puk) {
    write(to!string(g) ~ " ");
  }
  writeln;

  writeln("F1 (phase unknown): ");
  foreach(g; f1_puk) {
    write(to!string(g) ~ " ");
  }
  writeln;

}


// simulate a sperm/egg cell from parent
// marker_map must be sorted
FounderIndex[] simulate_meiotic_product(in Marker[] marker_map, in TrueGenotype[] parent,
                                        in uint m, in double p, ref Random gen)
in {
  assert(marker_map.length == parent.length, "marker_map and parent must be the same length");
  assert(p >= 0 && p <= 1, "p must be in [0,1]");
}
body {
  auto genotypes_on_product = new FounderIndex[](parent.length);

  double start_cM = marker_map[0].get_position;
  double end_cM = marker_map[$-1].get_position;

  // locations of crossovers
  double[] xo_locations = meiosis(start_cM, end_cM, m, p, gen);

  // draw genotype at first position
  uint current_strand = cast(uint)dice(gen, 1, 1);
  genotypes_on_product[0] = parent[0].get_allele(current_strand);

  if(xo_locations.length == 0) { // no XOs: fill in with that strand
    foreach(i; 1..parent.length)
      genotypes_on_product[i] = parent[i].get_allele(current_strand);
    return(genotypes_on_product);
  }

  auto current_marker = 1;
  foreach(xoloc; xo_locations) {
    // loop markers up to that crossover; paste with current strand
    // then switch strand
    while(marker_map[current_marker].get_position < xoloc) {
      genotypes_on_product[current_marker] = parent[current_marker].get_allele(current_strand);
      current_marker++;

      assert(current_marker < parent.length,
             "We shouldn't get here; all xo_locations should be < end_cM");
    }
    current_strand = 1-current_strand;
  }

  // fill in the rest
  foreach(i; current_marker..parent.length) {
    genotypes_on_product[i] = parent[i].get_allele(current_strand);
  }
  return(genotypes_on_product);
}

// combine two meiotic products to create a new individual
//     (I was tempted to say "fertilize")
TrueGenotype[] combine_meiotic_products(in FounderIndex[] egg, in FounderIndex[] sperm)
in {
  assert(egg.length == sperm.length, "egg and sperm must have same length");
}
body {
  TrueGenotype[] genotypes;

  foreach(i, eggv; egg) {
    auto g = new TrueGenotype(eggv, sperm[i]);
    genotypes ~= g;
  }

  return(genotypes);
}


unittest {
  auto chromosome = get_chromosome_with_id("1");
  auto marker_map = generate_map_eqspacing_onechr(100.0, 26, 1, chromosome);

  FounderIndex[] founders = [0, 1];

  auto inbred1 = generate_founder_genotypes_onechr(marker_map, founders[0]);
  auto inbred2 = generate_founder_genotypes_onechr(marker_map, founders[1]);
  auto f1 = generate_f1_genotypes_onechr(marker_map, founders);

  Random gen;
  gen.seed(unpredictableSeed);

  writeln("Simulated meiotic products:");
  auto meiotic_product = simulate_meiotic_product(marker_map, f1, 0, 0.0, gen);
  foreach(g; meiotic_product) write(g);
  writeln;
  meiotic_product = simulate_meiotic_product(marker_map, f1, 0, 0.0, gen);
  foreach(g; meiotic_product) write(g);
  writeln;

  writeln("Simulate backcross:");
  auto parent0_chr = generate_founder_chromosome(marker_map, founders[0]);
  meiotic_product = simulate_meiotic_product(marker_map, f1, 10, 0.0, gen);
  auto backcross = combine_meiotic_products(meiotic_product, parent0_chr);
  foreach(g; backcross) write(g, " ");
  writeln;
  meiotic_product = simulate_meiotic_product(marker_map, f1, 10, 0.0, gen);
  backcross = combine_meiotic_products(meiotic_product, parent0_chr);
  foreach(g; backcross) write(g, " ");
  writeln;

  writeln("Simulate intercross:");
  auto meiotic_product1 = simulate_meiotic_product(marker_map, f1, 10, 0.0, gen);
  auto meiotic_product2 = simulate_meiotic_product(marker_map, f1, 10, 0.0, gen);
  auto intercross = convert_truegenotype_to_genotypecombinator(combine_meiotic_products(meiotic_product1, meiotic_product2));
  foreach(g; intercross) write(g, " ");
  writeln;
  meiotic_product1 = simulate_meiotic_product(marker_map, f1, 10, 0.0, gen);
  meiotic_product2 = simulate_meiotic_product(marker_map, f1, 10, 0.0, gen);
  intercross = convert_truegenotype_to_genotypecombinator(combine_meiotic_products(meiotic_product1, meiotic_product2));
  foreach(g; intercross) write(g, " ");
  writeln;

  marker_map = generate_map_eqspacing_onechr(10000.0, 26, 1, chromosome);

  inbred1 = generate_founder_genotypes_onechr(marker_map, founders[0]);
  inbred2 = generate_founder_genotypes_onechr(marker_map, founders[1]);
  f1 = generate_f1_genotypes_onechr(marker_map, founders);

  writeln("Simulated meiotic products:");
  meiotic_product = simulate_meiotic_product(marker_map, f1, 0, 0.0, gen);
  foreach(g; meiotic_product) write(g);
  writeln;
  meiotic_product = simulate_meiotic_product(marker_map, f1, 0, 0.0, gen);
  foreach(g; meiotic_product) write(g);
  writeln;

  writeln("Simulate backcross:");
  parent0_chr = generate_founder_chromosome(marker_map, founders[0]);
  meiotic_product = simulate_meiotic_product(marker_map, f1, 10, 0.0, gen);
  backcross = combine_meiotic_products(meiotic_product, parent0_chr);
  foreach(g; backcross) write(g, " ");
  writeln;
  meiotic_product = simulate_meiotic_product(marker_map, f1, 10, 0.0, gen);
  backcross = combine_meiotic_products(meiotic_product, parent0_chr);
  foreach(g; backcross) write(g, " ");
  writeln;

  writeln("Simulate intercross:");
  meiotic_product1 = simulate_meiotic_product(marker_map, f1, 10, 0.0, gen);
  meiotic_product2 = simulate_meiotic_product(marker_map, f1, 10, 0.0, gen);
  intercross = convert_truegenotype_to_genotypecombinator(combine_meiotic_products(meiotic_product1, meiotic_product2));
  foreach(g; intercross) write(g, " ");
  writeln;
  meiotic_product1 = simulate_meiotic_product(marker_map, f1, 10, 0.0, gen);
  meiotic_product2 = simulate_meiotic_product(marker_map, f1, 10, 0.0, gen);
  intercross = convert_truegenotype_to_genotypecombinator(combine_meiotic_products(meiotic_product1, meiotic_product2));
  foreach(g; intercross) write(g, " ");
  writeln;
}

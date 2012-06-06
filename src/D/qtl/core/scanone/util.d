/**
 * Utility functions for scanone
 */

module qtl.core.scanone.util;

import std.stdio;
import std.algorithm;
import std.container;
import std.variant;
import std.conv;
import std.typecons;
import std.random;

import qtl.core.chromosome;
import qtl.core.primitives;
import qtl.core.genotype;

// fill up X matrix for scanone, one position
double[] create_scanone_Xmatrix(in Probability[][] genoprobs, in double[][] addcovar,
                                in double [][] intcovar, in double[] weights)
{
  if(genoprobs.length != weights.length)
    throw new Exception("Mismatch in no. individuals in genoprobs and weights");
  if(addcovar.length > 0 && addcovar[0].length > 0 && genoprobs.length != addcovar.length)
    throw new Exception("Mismatch in no. individuals in genoprobs and addcovar");
  if(intcovar.length > 0 && intcovar[0].length > 0 && genoprobs.length != intcovar.length)
    throw new Exception("Mismatch in no. individuals in genoprobs and intcovar");
  // genoprobs is [individuals][genotypes]
  // addcovar, intcovar are [individuals][covariates]

  auto n_ind = genoprobs.length;
  auto n_gen = genoprobs[0].length;

  size_t n_addcovar=0, n_intcovar=0;

  if(addcovar.length > 0) n_addcovar = addcovar[0].length;
  if(intcovar.length > 0) n_intcovar = intcovar[0].length;

  // X matrix
  auto ncolx = n_gen + n_addcovar + (n_gen-1)*n_intcovar;
  auto Xmatrix = new double[](n_ind * ncolx);

  size_t i, j, k1, k2, s;

  for(i=0; i<n_ind; i++) {
    // genotype probability columns
    for(j=0; j<n_gen; j++)
      Xmatrix[i+j*n_ind] = genoprobs[i][j]*weights[i];

    // addcovar columns
    for(j=0; j<n_addcovar; j++)
      Xmatrix[i+(j+n_gen)*n_ind] = addcovar[i][j]*weights[i];

    // intcovar columns
    for(k1=0,s=0; k1<n_gen-1; k1++)
      for(k2=0; k2<n_intcovar; k2++,s++)
        Xmatrix[i+(n_gen+n_addcovar+s)*n_ind] = genoprobs[i][k1]*intcovar[i][k2]*weights[i];
  }

  return Xmatrix;
}


FounderIndex[] get_sorted_founder_alleles(in TrueGenotype[] all_true_genotypes)
{
  // create vector of founder alleles
  size_t[FounderIndex] founders_hash;

  foreach(g; all_true_genotypes) {
    founders_hash[g.founders[0]] = 1;
    founders_hash[g.founders[1]] = 1;
  }

  FounderIndex[] founders_unique = founders_hash.keys;
  sort(founders_unique);

  return founders_unique;
}

unittest {
  import qtl.core.hmm.bc;
  import qtl.core.hmm.f2;

  auto f2_geno = allTrueGeno_F2();
  auto f2_founders = get_sorted_founder_alleles(f2_geno);
  write("F2 founder alleles:\t");
  foreach(f; f2_founders)
    write(f, " ");
  writeln;

  auto f2pk_geno = allTrueGeno_F2PK();
  auto f2pk_founders = get_sorted_founder_alleles(f2pk_geno);
  write("F2PK founder alleles:\t");
  foreach(f; f2pk_founders)
    write(f, " ");
  writeln;

  auto bc_geno = allTrueGeno_BC();
  auto bc_founders = get_sorted_founder_alleles(bc_geno);
  write("BC founder alleles:\t");
  foreach(f; bc_founders)
    write(f, " ");
  writeln;
}

// create table of possible genotyeps vs frequency of founder alleles
double[][] create_allele_freq_table(in TrueGenotype[] all_true_genotypes)
{
  // sorted founder alleles, and hash with indices
  auto founders = get_sorted_founder_alleles(all_true_genotypes);
  size_t[FounderIndex] founders_index;
  foreach(i, f; founders)
    founders_index[f] = i;

  // begin with matrix of 0's
  auto allele_freq_table = new double[][](all_true_genotypes.length, founders.length);
  foreach(i; 0..allele_freq_table.length)
    foreach(j; 0..allele_freq_table[i].length)
      allele_freq_table[i][j] = 0.0;

  foreach(ig, g; all_true_genotypes) {
    allele_freq_table[ig][founders_index[g.founders[0]]] += 0.5;
    allele_freq_table[ig][founders_index[g.founders[1]]] += 0.5;
  }

  return allele_freq_table;
}

unittest {
  import qtl.core.hmm.bc;
  import qtl.core.hmm.f2;

  auto f2_geno = allTrueGeno_F2();
  auto f2_allelefreq = create_allele_freq_table(f2_geno);
  writeln("F2 founder allele freq:");
  foreach(i; 0..f2_allelefreq.length) {
    foreach(j; 0..f2_allelefreq[i].length)
      writef("%3.1f ", f2_allelefreq[i][j]);
    writeln;
  }

  auto f2pk_geno = allTrueGeno_F2PK();
  auto f2pk_allelefreq = create_allele_freq_table(f2pk_geno);
  writeln("F2PK founder allele freq:");
  foreach(i; 0..f2pk_allelefreq.length) {
    foreach(j; 0..f2pk_allelefreq[i].length)
      writef("%3.1f ", f2pk_allelefreq[i][j]);
    writeln;
  }

  auto bc_geno = allTrueGeno_BC();
  auto bc_allelefreq = create_allele_freq_table(bc_geno);
  writeln("BC founder allele freq:");
  foreach(i; 0..bc_allelefreq.length) {
    foreach(j; 0..bc_allelefreq[i].length)
      writef("%3.1f ", bc_allelefreq[i][j]);
    writeln;
  }

}

// Collapse genotype probabilities to allele probabilities
double[][][] collapse_geno_prob_to_allele_prob(in double[][][] genoprobs, in double[][] allele_freq_table)
in {
  assert(genoprobs[0][0].length == allele_freq_table.length);
}
body {
  auto n_pos = genoprobs.length;
  auto n_ind = genoprobs[0].length;
  auto n_geno = allele_freq_table.length;
  auto n_allele = allele_freq_table[0].length;

  auto alleleprobs = new double[][][](n_pos, n_ind, n_allele);
  foreach(i; 0..n_pos) {
    foreach(j; 0..n_ind) {
      foreach(k; 0..n_allele) {
        alleleprobs[i][j][k] = 0.0;
        foreach(g; 0..n_geno)
          alleleprobs[i][j][k] += genoprobs[i][j][g]*allele_freq_table[g][k];
      }
    }
  }

  return alleleprobs;
}

// Create allele table then collapse genotype probabilities to allele probabilities
double[][][] collapse_geno_prob_to_allele_prob(in double[][][] genoprobs, in TrueGenotype[] all_true_genotypes)
in {
  assert(genoprobs[0][0].length == all_true_genotypes.length);
}
body {
  // founder allele frequency table
  auto allele_freq_table = create_allele_freq_table(all_true_genotypes);

  return collapse_geno_prob_to_allele_prob(genoprobs, allele_freq_table);
}


// find maximum lod score and the position at which it occurred
// if multiple positions share maximum, return a random one
Tuple!(double, Marker) get_peak_scanone(double[] lod, Marker[] map, ref Random gen)
in {
  assert(lod.length > 0, "lod should have length > 0");
  assert(lod.length == map.length, "lod and map should have the same length");
}
body {
  if(lod.length==1)
    return Tuple!(double, Marker)(lod[0], map[0]);

  double maxlod = lod[0];
  size_t[] index_maxlod = [0];

  foreach(i; 1..lod.length) {
    if(lod[i] == maxlod)
      index_maxlod ~= i;
    if(lod[i] > maxlod) {
      index_maxlod = [i];
      maxlod = lod[i];
    }
  }

  if(index_maxlod.length > 1)
    index_maxlod[0] = index_maxlod[uniform(0, index_maxlod.length, gen)];

  return Tuple!(double, Marker)(lod[index_maxlod[0]], map[index_maxlod[0]]);
}

// same with self-initiated random seed
Tuple!(double, Marker) get_peak_scanone(double[] lod, Marker[] map)
{
  Random gen;
  gen.seed(unpredictableSeed);
  return get_peak_scanone(lod, map, gen);
}


// same for doubly-index LOD scores
Tuple!(double, Marker)[] get_peak_scanone(double[][] lod, Marker[] map, ref Random gen)
in {
  assert(lod.length == map.length, "lod and map should have the same length");
}
body {
  Tuple!(double, Marker)[] peaks;

  foreach(i; 0..lod[0].length) {
    auto thislod = new double[](lod.length);
    foreach(j; 0..lod.length)
      thislod[j] = lod[j][i];

    peaks ~= get_peak_scanone(thislod, map, gen);
  }

  return peaks;
}

// same with self-initiated random seed
Tuple!(double, Marker)[] get_peak_scanone(double[][] lod, Marker[] map)
{
  Random gen;
  gen.seed(unpredictableSeed);
  return get_peak_scanone(lod, map, gen);
}

unittest {
  writeln("Unit test " ~ __FILE__);

  import qtl.core.simulate.map;
  import qtl.core.chromosome;

  // singly indexed lod
  double[] lod = [1.0, 3.0, 2.3, 5.6, 1.8];

  auto chromosome = get_chromosome_with_id("1");
  auto map = generate_map_eqspacing_onechr(100.0, cast(uint)lod.length, 1, chromosome);

  auto peak = get_peak_scanone(lod, map);
  assert(peak[0] == 5.6 && peak[1].name == "m4");

  // doubly indexed lod
  double[][] lod_2d = [[4.0, 5.0], [5.0, 2.6], [4.9, 0.3], [4.9, 5.3], [2.8, 2.8]];
  auto peak_2d = get_peak_scanone(lod_2d, map);
  assert(peak_2d[0][0] == 5.0 && peak_2d[0][1].name == "m2");
  assert(peak_2d[1][0] == 5.3 && peak_2d[1][1].name == "m4");

  // multiple positions sharing maximum LOD
  lod = [1.0, 5.6, 2.3, 5.6, 1.8];

  Random gen;
  gen.seed(51118652);
  peak = get_peak_scanone(lod, map, gen);
  assert(peak[0] == 5.6 && peak[1].name == "m4");

  gen.seed(20276563);
  peak = get_peak_scanone(lod, map, gen);
  assert(peak[0] == 5.6 && peak[1].name == "m2");
}

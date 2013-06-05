/**
 * calc_genoprob
 */

module qtl.core.hmm.calcgenoprob;

import std.string;
import std.conv;
import std.stdio;
import std.math;
import qtl.core.primitives, qtl.core.genotype;
import qtl.core.hmm.util;
import qtl.core.hmm.forwardbackward;
import qtl.core.hmm.cross;

// calculate QTL genotype probabilities
double[][][] calc_geno_prob(T)(Cross cross,
                               GenotypeCombinator[][] genotypes,
                               bool is_X_chr,
                               bool[] is_female,
                               T[] cross_direction,
                               Marker[] marker_map,
                               double[] rec_frac,
                               double error_prob)
{
  if(marker_map.length != rec_frac.length+1) {
    throw new Exception("no. positions in marker map doesn't match rec_frac length");
  }
  if(error_prob < 0.0 || error_prob > 1.0)
    throw new Exception("error_prob out of range");
  foreach(rf; rec_frac) {
    if(rf < 0 || rf > 0.5)
      throw new Exception("rec_frac must be >= 0 and <= 0.5");
  }
  if(is_female.length != genotypes.length)
    throw new Exception("is_female should be same length as genotypes");
  if(cross_direction.length != genotypes.length)
    throw new Exception("cross_direction should be same length as genotypes");

  size_t n_individuals = genotypes.length;
  size_t n_positions = marker_map.length;

  auto all_true_geno = cross.all_true.geno(is_X_chr);

  auto alpha = new double[][](all_true_geno.length, n_positions);
  auto beta = new double[][](all_true_geno.length, n_positions);
  auto genoprobs = new double[][][](n_positions, n_individuals, all_true_geno.length);

  foreach(ind; 0..n_individuals) {
    alpha = forwardEquations(cross, genotypes[ind], is_X_chr, is_female[ind], cross_direction[ind], marker_map, rec_frac, error_prob);
    beta = backwardEquations(cross, genotypes[ind], is_X_chr, is_female[ind], cross_direction[ind], marker_map, rec_frac, error_prob);

    // calculate genotype probabilities
    double sum_at_pos;
    auto possible_true_geno_index = cross.possible_true_geno_index(is_X_chr, is_female[ind], cross_direction[ind]);
    auto first = possible_true_geno_index[0];
    foreach(pos; 0..n_positions) {
      sum_at_pos = genoprobs[pos][ind][0] = alpha[first][pos] + beta[first][pos];
      foreach(i; possible_true_geno_index[1 .. possible_true_geno_index.length]) {
        auto true_geno = all_true_geno[i];
        genoprobs[pos][ind][i] = alpha[i][pos] + beta[i][pos];
        sum_at_pos = addlog(sum_at_pos, genoprobs[pos][ind][i]);
      }
      foreach(i; possible_true_geno_index) {
        genoprobs[pos][ind][i] = exp(genoprobs[pos][ind][i] - sum_at_pos);
      }
    }
  }
  return genoprobs;
}

// version for autosome with no cross direction information
double[][][] calc_geno_prob(Cross cross,
                            GenotypeCombinator[][] genotypes,
                            Marker[] marker_map,
                            double[] rec_frac,
                            double error_prob)
{
  auto is_female  = new bool[](genotypes.length);
  auto cross_direction = new bool[](genotypes.length);

  return(calc_geno_prob(cross, genotypes, 
                        false, is_female, cross_direction,
                        marker_map, rec_frac, error_prob));
}

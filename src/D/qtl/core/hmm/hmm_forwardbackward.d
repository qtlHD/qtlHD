/**
 * hmm_forwardbackward
 */

module qtl.core.hmm.hmm_forwardbackward;

import std.stdio;
import std.math;

import qtl.core.map.genetic_map_functions;
import qtl.core.primitives, qtl.core.genotype;
import qtl.core.hmm.hmm_util;

// forward Equations
double[][] forwardEquations(alias init, alias emit, alias step)(GenotypeCombinator[] genotypes,
                                                                TrueGenotype[] all_true_geno,
                                                                Marker[] marker_map,
                                                                double[] rec_frac,
                                                                double error_prob)
{
  size_t n_positions = marker_map.length;

  auto alpha = new double[][](all_true_geno.length, n_positions);

  // initialize alphas
  foreach(i, true_geno; all_true_geno) {
    if(isPseudoMarker(marker_map[0]))
      alpha[i][0] = init(true_geno);
    else
      alpha[i][0] = init(true_geno) + emit(genotypes[marker_map[0].id], true_geno, error_prob);
  }

  foreach(pos; 1 .. n_positions) {
    foreach(j, true_geno_right; all_true_geno) {
      alpha[j][pos] = alpha[0][pos-1] +
        step(all_true_geno[0], true_geno_right, rec_frac[pos-1]);

     foreach(i, true_geno_left; all_true_geno[1..$]) {
       alpha[j][pos] = addlog(alpha[j][pos],
                              alpha[i+1][pos-1] +
                              step(true_geno_left, true_geno_right, rec_frac[pos-1]));
     }
     if(!isPseudoMarker(marker_map[pos]))
       alpha[j][pos] += emit(genotypes[marker_map[pos].id], true_geno_right, error_prob);
    }
  }

  return alpha;
}



// backward Equations
double[][] backwardEquations(alias init, alias emit, alias step)(GenotypeCombinator[] genotypes,
                                                                 TrueGenotype[] all_true_geno,
                                                                 Marker[] marker_map,
                                                                 double[] rec_frac,
                                                                 double error_prob)
{
  size_t n_positions = marker_map.length;

  auto beta = new double[][](all_true_geno.length,n_positions);

  // initialize beta
  foreach(i, true_geno; all_true_geno) {
    beta[i][n_positions-1] = 0.0;
  }

  // backward equations
  for(int pos = cast(int)n_positions-2; pos >= 0; pos--) {
    foreach(i, true_geno_left; all_true_geno) {
      if(isPseudoMarker(marker_map[pos+1]))
        beta[i][pos] = beta[0][pos+1] +
          step(true_geno_left, all_true_geno[0], rec_frac[pos]);
      else
        beta[i][pos] = beta[0][pos+1] +
          step(true_geno_left, all_true_geno[0], rec_frac[pos]) +
          emit(genotypes[marker_map[pos+1].id], all_true_geno[0], error_prob);

      foreach(j, true_geno_right; all_true_geno[1..$]) {
        if(isPseudoMarker(marker_map[pos+1]))
          beta[i][pos] = addlog(beta[i][pos],
                                beta[j+1][pos+1] +
                                step(true_geno_left, true_geno_right, rec_frac[pos]));
        else
          beta[i][pos] = addlog(beta[i][pos],
                                beta[j+1][pos+1] +
                                step(true_geno_left, true_geno_right, rec_frac[pos])+
                                emit(genotypes[marker_map[pos+1].id], true_geno_right, error_prob));
     }
    }
  }

  return beta;
}

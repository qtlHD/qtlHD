/**
 * forward/backward equations
 */

module qtl.core.hmm.forwardbackward;

import std.stdio;
import std.math;

import qtl.core.map.genetic_map_functions;
import qtl.core.primitives, qtl.core.genotype;
import qtl.core.hmm.util;

// forward Equations
double[][] forwardEquations(alias init, alias emit, alias step)(in GenotypeCombinator[] genotypes,
                                                                in TrueGenotype[] all_true_geno,
                                                                in Marker[] marker_map,
                                                                in Probability[] rec_frac,
                                                                in Probability error_prob)
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
    foreach(ir, true_geno_right; all_true_geno) {
      alpha[ir][pos] = alpha[0][pos-1] +
        step(all_true_geno[0], true_geno_right, rec_frac[pos-1]);

      foreach(il; 1 .. all_true_geno.length) {
        auto true_geno_left = all_true_geno[il];
        alpha[ir][pos] = addlog(alpha[ir][pos],
                               alpha[il][pos-1] +
                               step(true_geno_left, true_geno_right, rec_frac[pos-1]));
     }
     if(!isPseudoMarker(marker_map[pos]))
       alpha[ir][pos] += emit(genotypes[marker_map[pos].id], true_geno_right, error_prob);
    }
  }

  return alpha;
}



// backward Equations
double[][] backwardEquations(alias init, alias emit, alias step)(in GenotypeCombinator[] genotypes,
                                                                 in TrueGenotype[] all_true_geno,
                                                                 in Marker[] marker_map,
                                                                 in Probability[] rec_frac,
                                                                 in Probability error_prob)
{
  size_t n_positions = marker_map.length;

  auto beta = new double[][](all_true_geno.length,n_positions);

  // initialize beta
  foreach(i, true_geno; all_true_geno) {
    beta[i][n_positions-1] = 0.0;
  }

  // backward equations
  for(size_t pos = n_positions-2; pos >= 0; pos--) {
    foreach(il, true_geno_left; all_true_geno) {
      if(isPseudoMarker(marker_map[pos+1]))
        beta[il][pos] = beta[0][pos+1] +
          step(true_geno_left, all_true_geno[0], rec_frac[pos]);
      else
        beta[il][pos] = beta[0][pos+1] +
          step(true_geno_left, all_true_geno[0], rec_frac[pos]) +
          emit(genotypes[marker_map[pos+1].id], all_true_geno[0], error_prob);

      foreach(ir; 1 .. all_true_geno.length) {
        auto true_geno_right = all_true_geno[ir];
        if(isPseudoMarker(marker_map[pos+1]))
          beta[il][pos] = addlog(beta[il][pos],
                                beta[ir][pos+1] +
                                step(true_geno_left, true_geno_right, rec_frac[pos]));
        else
          beta[il][pos] = addlog(beta[il][pos],
                                beta[ir][pos+1] +
                                step(true_geno_left, true_geno_right, rec_frac[pos])+
                                emit(genotypes[marker_map[pos+1].id], true_geno_right, error_prob));
     }
    }
  }

  return beta;
}

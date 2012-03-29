/**
 * hmm_forwardbackward
 */

module qtl.core.hmm.hmm_forwardbackward;

import std.stdio;
import std.math;

import qtl.core.map.genetic_map_functions;
import qtl.core.primitives, qtl.core.genotype;

// forward Equations
double[][] forwardEquations(alias init, alias emit, alias step)(in GenotypeCombinator[] genotypes, 
                                                                in TrueGenotype[] all_true_geno,
                                                                in Marker[] marker_map, 
                                                                in double[] rec_frac, 
                                                                in double error_prob)
{
  size_t n_positions = marker_map.length;

  auto alpha = new double[][](all_true_geno.length, n_positions);

  // initialize alphas
  foreach(true_geno; all_true_geno) {
    if(isPseudoMarker(marker_map[0])) 
      alpha[true_geno][0] = init(true_geno);
    else
      alpha[true_geno][0] = init(true_geno) + emit(genotypes[marker_map[0].id], true_geno, error_prob);
  }

  foreach(pos; 1 .. n_positions) {
    foreach(true_geno_right; all_true_geno) {

     alpha[true_geno_right][pos] = alpha[all_true_geno[0]][pos-1] + 
       step(all_true_geno[0], true_geno_right, rec_frac[pos-1]);

     foreach(true_geno_left; all_true_geno[1..$]) {
       alpha[true_geno_right][pos] = addlog(alpha[true_geno_right][pos], 
                                        alpha[true_geno_left][pos-1] + 
                                        step(true_geno_left, true_geno_right, rec_frac[pos-1]));
     }
     if(!isPseudoMarker(marker_map[pos]))
       alpha[true_geno_right][pos] += emit(genotypes[marker_map[pos].id], true_geno_right, error_prob);
    }
  }
  return alpha;
}



// backward Equations 
double[][] backwardEquations(alias init, alias emit, alias step)(in GenotypeCombinator[] genotypes, 
                                                                 in TrueGenotype[] all_true_geno,
                                                                 in Marker[] marker_map, 
                                                                 in double[] rec_frac, 
                                                                 in double error_prob)
{
  size_t n_positions = marker_map.length;

  auto beta = new double[][](all_true_geno.length,n_positions);

  // initialize beta
  foreach(true_geno; all_true_geno) {
    beta[true_geno][n_markers-1] = 0.0;
  }

  // backward equations
  for(int pos = cast(int)n_positions-2; pos >= 0; pos--) {
    foreach(true_geno_left; all_true_geno) {
      if(isPseudoMarker(marker_map[pos+1])) 
        beta[true_geno_left][pos] = beta[all_true_geno[0]][pos+1] + 
          step(true_geno_left, all_true_geno[0], rec_frac[pos]);
      else
        beta[true_geno_left][pos] = beta[all_true_geno[0]][pos+1] + 
          step(true_geno_left, all_true_geno[0], rec_frac[pos]) + 
          emit(genotypes[marker_map[pos+1].id], all_true_geno[0], error_prob);

     foreach(true_geno_right; all_true_geno[1..$]) {
       if(isPseudoMarker(marker_map[pos+1]))
         beta[true_geno_left][pos] = addlog(beta[true_geno_left][pos], 
                                            beta[true_geno_right][pos+1] + 
                                            step(true_geno_left, true_geno_right, rec_frac[pos]));
       else
         beta[true_geno_left][pos] = addlog(beta[true_geno_left][pos], 
                                            beta[true_geno_right][pos+1] + 
                                            step(true_geno_left, true_geno_right, rec_frac[pos])+
                                            emit(genotypes[marker_map[pos+1].id], true_geno_right, error_prob));
     }
    }
  }

  return beta;
}

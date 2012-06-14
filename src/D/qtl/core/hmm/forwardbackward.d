/**
 * forward/backward equations
 */

module qtl.core.hmm.forwardbackward;

import std.stdio;
import std.math;
import std.string;
import std.conv;

import qtl.core.map.genetic_map_functions;
import qtl.core.primitives, qtl.core.genotype;
import qtl.core.hmm.util;
import qtl.core.hmm.cross;

// forward Equations
double[][] forwardEquations(Cross cross,
                            GenotypeCombinator[] genotypes,
                            Marker[] marker_map,
                            Probability[] rec_frac,
                            Probability error_prob)
in {
  foreach(i, m; marker_map)
    assert(isPseudoMarker(m) || (m.id >= 0 && m.id < genotypes.length),
           "marker " ~ m.name ~ " [" ~ to!string(i) ~ "] has id out of range: " ~ to!string(m.id));
}
body {
  size_t n_positions = marker_map.length;

  auto alpha = new double[][](cross.all_true_geno.length, n_positions);

  // initialize alphas
  foreach(i, true_geno; cross.all_true_geno) {
    if(isPseudoMarker(marker_map[0]))
      alpha[i][0] = cross.init(true_geno);
    else
      alpha[i][0] = cross.init(true_geno) + cross.emit(genotypes[marker_map[0].id], true_geno, error_prob);
  }

  foreach(pos; 1 .. n_positions) {
    foreach(ir, true_geno_right; cross.all_true_geno) {
      alpha[ir][pos] = alpha[0][pos-1] +
        cross.step(cross.all_true_geno[0], true_geno_right, rec_frac[pos-1]);

      foreach(il; 1 .. cross.all_true_geno.length) {
        auto true_geno_left = cross.all_true_geno[il];
        alpha[ir][pos] = addlog(alpha[ir][pos],
                               alpha[il][pos-1] +
                               cross.step(true_geno_left, true_geno_right, rec_frac[pos-1]));
     }
     if(!isPseudoMarker(marker_map[pos]))
       alpha[ir][pos] += cross.emit(genotypes[marker_map[pos].id], true_geno_right, error_prob);
    }
  }

  return alpha;
}



// backward Equations
double[][] backwardEquations(Cross cross,
                             GenotypeCombinator[] genotypes,
                             Marker[] marker_map,
                             Probability[] rec_frac,
                             Probability error_prob)
in {
  foreach(i, m; marker_map)
    assert(isPseudoMarker(m) || (m.id >= 0 && m.id < genotypes.length),
           "marker " ~ m.name ~ " [" ~ to!string(i) ~ "] has id out of range: " ~ to!string(m.id));
}
body {
  size_t n_positions = marker_map.length;

  auto beta = new double[][](cross.all_true_geno.length,n_positions);

  // initialize beta
  foreach(i, true_geno; cross.all_true_geno) {
    beta[i][n_positions-1] = 0.0;
  }

  // backward equations
  for(int pos = cast(int)n_positions-2; pos >= 0; pos--) {
    foreach(il, true_geno_left; cross.all_true_geno) {
      if(isPseudoMarker(marker_map[pos+1]))
        beta[il][pos] = beta[0][pos+1] +
          cross.step(true_geno_left, cross.all_true_geno[0], rec_frac[pos]);
      else
        beta[il][pos] = beta[0][pos+1] +
          cross.step(true_geno_left, cross.all_true_geno[0], rec_frac[pos]) +
          cross.emit(genotypes[marker_map[pos+1].id], cross.all_true_geno[0], error_prob);

      foreach(ir; 1 .. cross.all_true_geno.length) {
        auto true_geno_right = cross.all_true_geno[ir];
        if(isPseudoMarker(marker_map[pos+1]))
          beta[il][pos] = addlog(beta[il][pos],
                                beta[ir][pos+1] +
                                cross.step(true_geno_left, true_geno_right, rec_frac[pos]));
        else
          beta[il][pos] = addlog(beta[il][pos],
                                beta[ir][pos+1] +
                                cross.step(true_geno_left, true_geno_right, rec_frac[pos])+
                                cross.emit(genotypes[marker_map[pos+1].id], true_geno_right, error_prob));
     }
    }
  }

  return beta;
}

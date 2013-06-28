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
                            GenotypeSymbolMapper[] genotypes,
                            bool is_X_chr,
                            bool is_female,
                            int[] cross_direction,
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

  // possible genotypes for this chromosome and then for this individual
  auto all_true_geno = cross.all_true_geno(is_X_chr);
  auto possible_true_geno_index = cross.possible_true_geno_index(is_X_chr, is_female, cross_direction);
  auto first = possible_true_geno_index[0];

  // to contain ln Pr(G_i = g | marker data)
  auto alpha = new double[][](all_true_geno.length, n_positions);

  // if not all genotypes possible, fill array with zeros
  if(possible_true_geno_index.length < all_true_geno.length)
    foreach(i, g; all_true_geno)
      foreach(p; 0..n_positions)
        alpha[i][p] = 0.0;

  // initialize alphas
  foreach(i; possible_true_geno_index) {
    auto true_geno = all_true_geno[i];
    if(isPseudoMarker(marker_map[0]))
      alpha[i][0] = cross.init(true_geno, is_X_chr, is_female, cross_direction);
    else
      alpha[i][0] = cross.init(true_geno, is_X_chr, is_female, cross_direction) + 
        cross.emit(genotypes[marker_map[0].id], true_geno, error_prob, is_X_chr, is_female, cross_direction);
  }

  foreach(pos; 1 .. n_positions) {
    foreach(ir; possible_true_geno_index) {
      auto true_geno_right = all_true_geno[ir];
      alpha[ir][pos] = alpha[first][pos-1] +
        cross.step(all_true_geno[first], true_geno_right, rec_frac[pos-1], is_X_chr, is_female, cross_direction);

      foreach(il; possible_true_geno_index[1..$]) {
        auto true_geno_left = all_true_geno[il];
        alpha[ir][pos] = addlog(alpha[ir][pos],
                                alpha[il][pos-1] +
                                cross.step(true_geno_left, true_geno_right, rec_frac[pos-1], is_X_chr, is_female, cross_direction));
     }
     if(!isPseudoMarker(marker_map[pos]))
       alpha[ir][pos] += cross.emit(genotypes[marker_map[pos].id], true_geno_right, error_prob, is_X_chr, is_female, cross_direction);
    }
  }

  return alpha;
}



// backward Equations
double[][] backwardEquations(Cross cross,
                             GenotypeSymbolMapper[] genotypes,
                             bool is_X_chr,
                             bool is_female,
                             int[] cross_direction,
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

  // possible genotypes for this chromosome and then for this individual
  auto all_true_geno = cross.all_true_geno(is_X_chr);
  auto possible_true_geno_index = cross.possible_true_geno_index(is_X_chr, is_female, cross_direction);
  auto first = possible_true_geno_index[0];

  auto beta = new double[][](all_true_geno.length,n_positions);

  // if not all genotypes possible, fill array with zeros
  if(possible_true_geno_index.length < all_true_geno.length)
    foreach(i, g; all_true_geno)
      foreach(p; 0..n_positions)
        beta[i][p] = 0.0;

  // initialize beta
  foreach(i; possible_true_geno_index) {
    beta[i][n_positions-1] = 0.0;
  }

  // backward equations
  for(int pos = cast(int)n_positions-2; pos >= 0; pos--) {
    foreach(il; possible_true_geno_index) {
      auto true_geno_left = all_true_geno[il];
      if(isPseudoMarker(marker_map[pos+1]))
        beta[il][pos] = beta[first][pos+1] +
          cross.step(true_geno_left, all_true_geno[first], rec_frac[pos], is_X_chr, is_female, cross_direction);
      else
        beta[il][pos] = beta[first][pos+1] +
          cross.step(true_geno_left, all_true_geno[first], rec_frac[pos], is_X_chr, is_female, cross_direction) +
          cross.emit(genotypes[marker_map[pos+1].id], all_true_geno[first], error_prob, is_X_chr, is_female, cross_direction);

      foreach(ir; possible_true_geno_index[1..$]) {
        auto true_geno_right = all_true_geno[ir];
        if(isPseudoMarker(marker_map[pos+1]))
          beta[il][pos] = addlog(beta[il][pos],
                                 beta[ir][pos+1] +
                                 cross.step(true_geno_left, true_geno_right, rec_frac[pos], is_X_chr, is_female, cross_direction));
        else
          beta[il][pos] = addlog(beta[il][pos],
                                 beta[ir][pos+1] +
                                 cross.step(true_geno_left, true_geno_right, rec_frac[pos], is_X_chr, is_female, cross_direction) +
                                 cross.emit(genotypes[marker_map[pos+1].id], true_geno_right, error_prob, is_X_chr, is_female, cross_direction));
     }
    }
  }

  return beta;
}

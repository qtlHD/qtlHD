/**
 * estmap
 */

module qtl.core.hmm.estmap;

import std.string;
import std.conv;
import std.stdio;
import std.math;
import qtl.core.primitives, qtl.core.genotype;
import qtl.core.hmm.util;
import qtl.core.hmm.forwardbackward;
import qtl.core.hmm.cross;

// re-estimate inter-marker recombination fractions
double[] estmap(T)(Cross cross,
                   GenotypeCombinator[][] genotypes,
                   bool is_X_chr,
                   bool[] is_female,
                   T[] cross_direction,
                   Marker[] marker_map,
                   double[] rec_frac,
                   double error_prob,
                   uint max_iterations,
                   double tol,
                   bool verbose)
{
    if(marker_map.length != rec_frac.length+1)
      throw new Exception("no. markers in marker map doesn't match rec_frac length");
    if(error_prob < 0.0 || error_prob > 1.0)
      throw new Exception("error_prob out of range");
    foreach(rf; rec_frac) {
      if(rf < 0 || rf > 0.5)
	throw new Exception("rec_frac must be >= 0 and <= 0.5");
    }
    if(max_iterations < 0)
      throw new Exception("max_iterations should be >= 0");
    if(tol < 0)
      throw new Exception("tol >= 0");
    if(is_female.length != genotypes.length)
      throw new Exception("is_female should be same length as genotypes")
        if(cross_direction.length != genotypes.length)
          throw new Exception("cross_direction should be same length as genotypes")

    auto crossPK = form_cross_phaseknown(cross);

    size_t n_individuals = genotypes.length;
    size_t n_markers = marker_map.length;

    auto cur_rec_frac = rec_frac.dup;
    auto prev_rec_frac = rec_frac.dup;
    auto all_true_geno = crossPK.all_true_geno(is_X_chr);
    double[][] alpha = new double[][](all_true_geno.length, n_markers);
    double[][] beta = new double[][](all_true_geno.length, n_markers);
    double[][] gamma = new double[][](all_true_geno.length, n_markers);
    double sum_gamma;
    foreach(it; 0..max_iterations) {
      foreach(ref rf; cur_rec_frac) {
	rf = 0.0;
      }

      foreach(ind; 0..n_individuals) {

	// forward and backward equations
	alpha = forwardEquations(crossPK, genotypes[ind], is_X_chr, is_female[ind], cross_direction[ind], marker_map, prev_rec_frac, error_prob);
	beta = backwardEquations(crossPK, genotypes[ind], is_X_chr, is_female[ind], cross_direction[ind], marker_map, prev_rec_frac, error_prob);

        // possible genotypes for this individual, indices to all_true_geno
        auto possible_true_geno_index = crossPK.possible_geno_index(is_X_chr, is_female[ind], cross_direction[ind]);

	foreach(j; 0..prev_rec_frac.length) {
	  // calculate gamma = log Pr(v1, v2, O)
	  auto sum_gamma_undef = true;
	  foreach(il; possible_true_geno_index) {
            left_gen = all_true_geno[il];
	    foreach(ir; possible_true_geno_index) {
              right_gen = all_true_geno[ir];
              if(isPseudoMarker(marker_map[j+1]))
                gamma[il][ir] = alpha[il][j] + beta[ir][j+1] +
                  crossPK.step(left_gen, right_gen, prev_rec_frac[j], is_X_chr, is_female[ind], cross_direction[ind]);
              else
                gamma[il][ir] = alpha[il][j] + beta[ir][j+1] +
                  crossPK.emit(genotypes[ind][marker_map[j+1].id], right_gen, error_prob, is_X_chr, is_female[ind], cross_direction[ind]) +
                  crossPK.step(left_gen, right_gen, prev_rec_frac[j], is_X_chr, is_female[ind], cross_direction[ind]);

	      if(sum_gamma_undef) {
		sum_gamma_undef = false;
		sum_gamma = gamma[il][ir];
	      }
	      else {
		sum_gamma = addlog(sum_gamma, gamma[il][ir]);
	      }
	    }
	  }

	  // update cur_rf
	  foreach(il; possible_true_geno_index) {
	    foreach(ir; possible_true_geno_index) {
	      cur_rec_frac[j] += crossPK.nrec(all_true_geno[il], all_true_geno[ir]) * exp(gamma[il][ir] - sum_gamma);
	    }
	  }
	} // loop over marker intervals

      } // loop over individuals

      // rescale
      foreach(ref rf; cur_rec_frac) {
	rf /= n_individuals;
	if(rf < tol/1000.0) rf = tol/1000.0;
	else if(rf > 0.5-tol/1000.0) rf = 0.5-tol/1000.0;
      }

      if(verbose) {
	auto maxdif=0.0;
	double tempdif;
	foreach(j; 0..prev_rec_frac.length) {
	  tempdif = abs(prev_rec_frac[j] - cur_rec_frac[j]);
	  if(tempdif > maxdif) {
	    maxdif = tempdif;
	  }
	}
	writefln("%4d %.12f", it, tempdif);
      }

      // check convergence
      auto converged = true;
      foreach(j; 0..prev_rec_frac.length) {
	if(abs(prev_rec_frac[j] - cur_rec_frac[j]) > tol*(cur_rec_frac[j]+tol*100.0)) {
	  converged = false;
	  break;
	}
      }

      if(converged) break;

      prev_rec_frac = cur_rec_frac.dup;
    }

    // calculate log likelihood
    auto loglik = 0.0;
    double curloglik;
    foreach(ind; 0..n_individuals) {

      alpha = forwardEquations(crossPK, genotypes[ind], is_X_chr, is_female[ind], cross_direction[ind], marker_map, prev_rec_frac, error_prob);
      auto possible_true_geno_index = crossPK.possible_geno_index(is_X_chr, is_female[ind], cross_direction[ind]);

      auto curloglik_undef = true;
      foreach(g_index; possible_true_geno_index) {
	if(curloglik_undef) {
	  curloglik_undef = false;
	  curloglik = alpha[g_index][prev_rec_frac.length-1];
	}
	else {
	  curloglik = addlog(curloglik, alpha[g_index][prev_rec_frac.length-1]);
	}
      }
      loglik += curloglik;
    }

    if(verbose) {
      writefln("loglik = %.3f", loglik);
    }

    return(cur_rec_frac);
}

// version for autosome with no cross direction information
double[] estmap(Cross cross,
                GenotypeCombinator[][] genotypes,
                Marker[] marker_map,
                double[] rec_frac,
                double error_prob,
                uint max_iterations,
                double tol,
                bool verbose)
{
  auto is_female  = new bool[](genotypes.length);
  auto cross_direction = new bool[](genotypes.length);

  return(estmap(cross, genotypes,
                false, is_female, cross_direction,
                marker_map, rec_frac, error_prob,
                max_iterations, tol, verbose));
}

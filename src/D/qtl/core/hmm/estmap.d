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
double[] estmap(in Cross cross,
                in GenotypeCombinator[][] genotypes,
                in Marker[] marker_map,
                in double[] rec_frac,
                in double error_prob,
                in uint max_iterations,
                in double tol,
                in bool verbose)
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

    auto crossPK = form_cross_phaseknown(cross);

    size_t n_individuals = genotypes.length;
    size_t n_markers = marker_map.length;

    auto cur_rec_frac = rec_frac.dup;
    auto prev_rec_frac = rec_frac.dup;
    double[][] alpha = new double[][](crossPK.all_true_geno.length, n_markers);
    double[][] beta = new double[][](crossPK.all_true_geno.length, n_markers);
    double[][] gamma = new double[][](crossPK.all_true_geno.length, n_markers);
    double sum_gamma;
    foreach(it; 0..max_iterations) {
      foreach(ref rf; cur_rec_frac) {
	rf = 0.0;
      }

      foreach(ind; 0..n_individuals) {

	// forward and backward equations
	alpha = forwardEquations(crossPK, genotypes[ind], marker_map, prev_rec_frac, error_prob);
	beta = backwardEquations(crossPK, genotypes[ind], marker_map, prev_rec_frac, error_prob);

	foreach(j; 0..prev_rec_frac.length) {
	  // calculate gamma = log Pr(v1, v2, O)
	  auto sum_gamma_undef = true;
	  foreach(il, left_gen; crossPK.all_true_geno) {
	    foreach(ir, right_gen; crossPK.all_true_geno) {
              if(isPseudoMarker(marker_map[j+1]))
                gamma[il][ir] = alpha[il][j] + beta[ir][j+1] +
                  crossPK.step(left_gen, right_gen, prev_rec_frac[j]);
              else
                gamma[il][ir] = alpha[il][j] + beta[ir][j+1] +
                  crossPK.emit(genotypes[ind][marker_map[j+1].id], right_gen, error_prob) +
                  crossPK.step(left_gen, right_gen, prev_rec_frac[j]);

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
	  foreach(il, left_gen; crossPK.all_true_geno) {
	    foreach(ir, right_gen; crossPK.all_true_geno) {
	      cur_rec_frac[j] += crossPK.nrec(left_gen, right_gen) * exp(gamma[il][ir] - sum_gamma);
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

      alpha = forwardEquations(crossPK, genotypes[ind], marker_map, prev_rec_frac, error_prob);

      auto curloglik_undef = true;
      foreach(g_index, gen; crossPK.all_true_geno) {
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






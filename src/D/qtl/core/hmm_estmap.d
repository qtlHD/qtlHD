/**
 * hmm_estmap
 */

module qtl.core.hmm_estmap;

// things I think I really need
import qtl.core.primitives;
import qtl.core.genotype;
import qtl.core.map_functions;
import std.stdio;
import std.math;
import qtl.core.hmm_f2;
import qtl.core.hmm_bc;
import qtl.core.hmm_util;

// things for the unit tests 
import qtl.plugins.input.read_csv;
import std.path;


// re-estimate inter-marker recombination fractions
double[] estmapF2(Genotype!F2[][] genotypes, double[] rec_frac, double error_prob,
		  int max_iterations, double tol, bool verbose)
{
  if(genotypes[0].length != rec_frac.length+1)
    throw new Exception("no. markers in genotypes doesn't match rec_frac length");
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

  int n_individuals = genotypes.length;
  int n_markers = genotypes[0].length;
  F2pk[] all_true_geno = [F2pk.AA, F2pk.AB, F2pk.BA, F2pk.BB];

  auto cur_rec_frac = rec_frac.dup; 
  double[int][F2pk] alpha, beta;
  double[F2pk][F2pk] gamma;
  double sum_gamma;
  foreach(it; 0..max_iterations) {
    foreach(ref rf; cur_rec_frac) {
      rf = 0.0;
    }

    foreach(ind; 0..n_individuals) {

      // forward and backward equations
      alpha = forwardEquationsF2pk(genotypes[ind], all_true_geno, rec_frac, error_prob);
      beta = backwardEquationsF2pk(genotypes[ind], all_true_geno, rec_frac, error_prob);


      foreach(j; 0..rec_frac.length) {
	// calculate gamma = log Pr(v1, v2, O)
	auto sum_gamma_undef = true;
	foreach(left_gen; all_true_geno) {
	  foreach(right_gen; all_true_geno) {
	    gamma[left_gen][right_gen] = alpha[left_gen][j] + beta[right_gen][j+1] + 
	      emitF2pk(genotypes[ind][j+1], right_gen, error_prob) +
	      stepF2pk(left_gen, right_gen, rec_frac[j]);

	    if(sum_gamma_undef) {
	      sum_gamma_undef = false;
	      sum_gamma = gamma[left_gen][right_gen];
	    }
	    else {
	      sum_gamma = addlog(sum_gamma, gamma[left_gen][right_gen]);
	    }
	  }
	}

	// update cur_rf
	foreach(left_gen; all_true_geno) {
	  foreach(right_gen; all_true_geno) {
	    cur_rec_frac[j] += nrecF2pk(left_gen, right_gen) * exp(gamma[left_gen][right_gen] - sum_gamma);
	  }
	}
      } /* loop over marker intervals */

    } /* loop over individuals */


    /* rescale */
    foreach(ref rf; cur_rec_frac) {
      rf /= n_individuals;
      if(rf < tol/1000.0) rf = tol/1000.0;
      else if(rf > 0.5-tol/1000.0) rf = 0.5-tol/1000.0;
    }

    if(verbose) {
      auto maxdif=0.0;
      double tempdif;
      foreach(j; 0..rec_frac.length) {
	tempdif = abs(rec_frac[j] - cur_rec_frac[j]);
	if(tempdif > maxdif) {
	  maxdif = tempdif;
	}
      }
      writefln("%4d %.12f", it, tempdif);
    }

    /* check convergence */
    auto converged = true;
    foreach(j; 0..rec_frac.length) {
      if(abs(rec_frac[j] - cur_rec_frac[j]) > tol*(cur_rec_frac[j]+tol*100.0)) {
	converged = false; 
	break;
      }
    }

    if(converged) break; 

    rec_frac = cur_rec_frac.dup;
  }
    
  /* calculate log likelihood */
  auto loglik = 0.0;
  double curloglik;
  foreach(ind; 0..n_individuals) {

    alpha = forwardEquationsF2pk(genotypes[ind], all_true_geno, rec_frac, error_prob);

    auto curloglik_undef = true;
    foreach(gen; all_true_geno) {
      if(curloglik_undef) {
	curloglik_undef = false;
	curloglik = alpha[gen][rec_frac.length-1];
      }
      else {
	curloglik = addlog(curloglik, alpha[gen][rec_frac.length-1]);
      }
    }
    loglik += curloglik;
  }

  if(verbose) {
    writefln("loglik = %.3f", loglik);
  }

  return(cur_rec_frac);
}


unittest {
  writeln("Unit test " ~ __FILE__);
  writeln("    unit test estmapF2:");
  alias std.path.join join;
  auto fn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","test","data","input","listeria.csv");
  writeln("      - read CSV " ~ fn);
  auto data = new ReadSimpleCSV!F2(fn);
  
  Marker[] markers_on_chr_4;
  writeln("    Grab markers");
  foreach(marker; data.markers) {
    if(marker.chromosome.name=="4") {
      markers_on_chr_4 ~= marker;
    }
  }

  writeln("      - Subset genotype data");
  Genotype!F2[][] chr_4_genotypes;
  chr_4_genotypes.reserve(data.genotypes.length);
  foreach(i; 0..data.genotypes.length) {
    Genotype!F2[] an_individuals_genotype;
    an_individuals_genotype.reserve(markers_on_chr_4.length);
    foreach(j; 0..markers_on_chr_4.length) {
      an_individuals_genotype ~= data.genotypes[i][markers_on_chr_4[j].id];
    }
    chr_4_genotypes ~= an_individuals_genotype;
  }

  writeln("      - Get recombination fractions");
  double[] dist_cM;
  foreach(i; 1..markers_on_chr_4.length) {
    dist_cM ~= markers_on_chr_4[i].position - markers_on_chr_4[i-1].position;
  }
  auto rec_frac = mapFunction(dist_cM, "kosambi");

  auto rec_frac_rqtl = [0.18274786564786985044,
			0.15620001633845906341,
			0.28772927391795305452];

  foreach(i; 0..rec_frac.length) {
    assert(abs(rec_frac[i] - rec_frac_rqtl[i]) < 1e-14);
  }

  auto rec_frac_rev_rqtl = [0.18726831033621754719,
			    0.15929622543721855266,
			    0.29234793156548449788];

  auto rec_frac_rev = estmapF2(chr_4_genotypes, rec_frac, 0.002, 100, 1e-6, true);

  foreach(i; 0..rec_frac.length) {
    assert(abs(rec_frac_rev[i] - rec_frac_rev_rqtl[i]) < 1e-5);
  }
}

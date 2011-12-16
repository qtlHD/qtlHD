/**
 * hmm_util
 */

module qtl.core.hmm.hmm_util;

import std.stdio;
import std.math;

import qtl.core.map.genetic_map_functions;
import qtl.core.primitives;
import qtl.core.hmm.hmm_f2;
import qtl.core.hmm.hmm_bc;

// Calculate addlog(a,b) = log[exp(a) + exp(b)]
double addlog(double a, double b)
{
  enum TOL = 200.0;

  if(b > a + TOL) return(b);
  else if(a > b + TOL) return(a);
  else return(a + log1p(exp(b-a)));
}

import std.conv;

unittest {
  writeln("Unit test " ~ __FILE__);
  writeln("    unit test addlog");
  double a=50, b=60, d=2;
  auto value = abs(addlog(a,b) - log(exp(a)+exp(b)));
  assert(abs(addlog(a,b) - log(exp(a)+exp(b))) < 1e-15, to!string(value));
  assert(abs(addlog(b,a) - log(exp(a)+exp(b))) < 1e-15);
  assert(abs(addlog(a,d) - log(exp(a)+exp(d))) < 1e-15);
  assert(abs(addlog(d,a) - log(exp(a)+exp(d))) < 1e-15);
  assert(abs(addlog(b,d) - log(exp(b)+exp(d))) < 1e-15);
  assert(abs(addlog(d,b) - log(exp(b)+exp(d))) < 1e-15);
  assert(addlog(a,a+300) == a+300);
  assert(addlog(a,a-300) == a);
  assert(addlog(a+300,a) == a+300);
  assert(addlog(a-300,a) == a);
}

// forward Equations
  double[][] forwardEquations(GT,PKGT)(in Genotype!GT[] genotypes, in PKGT[] all_true_geno, 
                                 in double[] rec_frac, in double error_prob) 
  {
    size_t n_markers = genotypes.length;

    auto alpha = new double[][](all_true_geno.length,n_markers);

    // initialize alphas
    foreach(true_geno; all_true_geno) {
      alpha[true_geno][0] = init(true_geno) + emit(genotypes[0], true_geno, error_prob);
    }

    foreach(pos; 1 .. n_markers) {
      foreach(true_geno_right; all_true_geno) {
 
       alpha[true_geno_right][pos] = alpha[all_true_geno[0]][pos-1] + 
         step(all_true_geno[0], true_geno_right, rec_frac[pos-1]);

       foreach(true_geno_left; all_true_geno[1..$]) {
         alpha[true_geno_right][pos] = addlog(alpha[true_geno_right][pos], 
                                          alpha[true_geno_left][pos-1] + 
                                          step(true_geno_left, true_geno_right, rec_frac[pos-1]));
       }
       alpha[true_geno_right][pos] += emit(genotypes[pos], true_geno_right, error_prob);
      }
    }
    return alpha;
  }



// backward Equations 
  double[][] backwardEquations(GT,PKGT)(in Genotype!GT[] genotypes, in PKGT[] all_true_geno, 
                                  in double[] rec_frac, in double error_prob) 
  {
    size_t n_markers = genotypes.length;

    auto beta = new double[][](all_true_geno.length,n_markers);

    // initialize beta
    foreach(true_geno; all_true_geno) {
      beta[true_geno][n_markers-1] = 0.0;
    }

    // backward equations
    for(int pos = cast(int)n_markers-2; pos >= 0; pos--) {
      foreach(true_geno_left; all_true_geno) {
       beta[true_geno_left][pos] = beta[all_true_geno[0]][pos+1] + 
         step(true_geno_left, all_true_geno[0], rec_frac[pos]) + 
         emit(genotypes[pos+1], all_true_geno[0], error_prob);

       foreach(true_geno_right; all_true_geno[1..$]) {
         beta[true_geno_left][pos] = addlog(beta[true_geno_left][pos], 
                                        beta[true_geno_right][pos+1] + 
                                        step(true_geno_left, true_geno_right, rec_frac[pos])+
                                        emit(genotypes[pos+1], true_geno_right, error_prob));
       }
      }
    }

    return beta;
  }

/**
 * HMM module
 */

module qtl.core.hmm;

// things I think I really need
import qtl.core.primitives;
import qtl.core.genotype;
import std.stdio;
import std.math;
import qtl.core.hmm_f2;

// things for the unit tests 
import qtl.plugins.input.read_csv;
import std.path;


// calculate QTL genotype probabilities
double[F2][int][int] calcGenoprobF2(Genotype!F2[][] genotypes, double[] rec_frac, double error_prob)
in {
  assert(genotypes.length == rec_frac.length-1);
  assert(error_prob >= 0 && error_prob <= 1);
 }
body {
  int n_individuals = genotypes.length;
  int n_markers = genotypes[0].length;
  F2[] all_true_geno = [F2.A, F2.H, F2.B];

  double[int][F2] alpha, beta;
  double[F2][int][int] genoprobs;

  foreach(ind; 0..n_individuals) {
    // initialize alpha and beta
    foreach(true_geno; all_true_geno) {
      alpha[true_geno][0] = initF2(true_geno) + emitF2(genotypes[ind][0], true_geno, error_prob);
      beta[true_geno][n_markers-1] = 0.0;
    }


    // forward equations 
    foreach(pos; 1 .. n_markers) {
      foreach(true_geno_right; all_true_geno) {

	alpha[true_geno_right][pos] = alpha[all_true_geno[0]][pos-1] + 
	  stepF2(all_true_geno[0], true_geno_right, rec_frac[pos-1]);

	foreach(true_geno_left; all_true_geno[1..$]) {
	  alpha[true_geno_right][pos] = addlog(alpha[true_geno_right][pos], 
					       alpha[true_geno_left][pos-1] + 
					       stepF2(true_geno_left, true_geno_right, rec_frac[pos-1]));
	}
	alpha[true_geno_right][pos] += emitF2(genotypes[ind][pos], true_geno_right, error_prob);
      }
    }

	
    // backward equations
    foreach(pos; n_markers-2 .. 0) {
      foreach(true_geno_left; all_true_geno) {
	beta[true_geno_left][pos] = beta[all_true_geno[0]][pos+1] + 
	  stepF2(true_geno_left, all_true_geno[0], rec_frac[pos]) + 
	  emitF2(genotypes[ind][pos+1], all_true_geno[0], error_prob);

	foreach(true_geno_right; all_true_geno[1..$]) {
	  beta[true_geno_left][pos] = addlog(beta[true_geno_left][pos], 
					     beta[true_geno_right][pos+1] + 
					     stepF2(true_geno_left, true_geno_right, rec_frac[pos])+
					     emitF2(genotypes[ind][pos+1], true_geno_right, error_prob));
	}

      }
    }

    /* calculate genotype probabilities */
    double sum_at_pos;
    foreach(pos; 0..n_markers) {
      sum_at_pos = genoprobs[ind][pos][all_true_geno[0]] = alpha[all_true_geno[0]][pos] + beta[all_true_geno[0]][pos];
      foreach(true_geno; all_true_geno[1..$]) {
	genoprobs[ind][pos][true_geno] = alpha[true_geno][pos] + beta[true_geno][pos];
	sum_at_pos = addlog(sum_at_pos, genoprobs[ind][pos][true_geno]);
      }
      foreach(true_geno; all_true_geno) {
	genoprobs[ind][pos][true_geno] = exp(genoprobs[ind][pos][true_geno] - sum_at_pos);
      }
    }
  }
  return genoprobs;
}

unittest {
  writeln("Unit test " ~ __FILE__);
  alias std.path.join join;
  auto fn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","test","data","input","listeria.csv");
  writeln("  - read CSV " ~ fn);
  auto data = new ReadSimpleCSV!F2(fn);
  //  auto cross = new F2Cross(data.genotypes);

  //  auto genoprob = calcGenoprob!F2(cross, 0, 0.01);
}




// Calculate addlog(a,b) = log[exp(a) + exp(b)]
double addlog(double a, double b)
{
  enum TOL = 200.0;

  if(b > a + TOL) return(b);
  else if(a > b + TOL) return(a);
  else return(a + log1p(exp(b-a)));
}

unittest {
  writeln("    unit test addlog");
  double a=50, b=60, d=2;
  assert(abs(addlog(a,b) - log(exp(a)+exp(b))) < 1e-16);
  assert(abs(addlog(b,a) - log(exp(a)+exp(b))) < 1e-16);
  assert(abs(addlog(a,d) - log(exp(a)+exp(d))) < 1e-16);
  assert(abs(addlog(d,a) - log(exp(a)+exp(d))) < 1e-16);
  assert(abs(addlog(b,d) - log(exp(b)+exp(d))) < 1e-16);
  assert(abs(addlog(d,b) - log(exp(b)+exp(d))) < 1e-16);
  assert(addlog(a,a+300) == a+300);
  assert(addlog(a,a-300) == a);
  assert(addlog(a+300,a) == a+300);
  assert(addlog(a-300,a) == a);
}
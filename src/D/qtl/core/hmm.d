/**
 * HMM module
 */

module qtl.core.hmm;

// things I think I really need
import qtl.core.primitives;
import qtl.core.genotype;
import qtl.core.mapFunctions;
import std.stdio;
import std.math;
import qtl.core.hmm_f2;

// things for the unit tests 
import qtl.plugins.input.read_csv;
import std.path;


// calculate QTL genotype probabilities
double[F2][int][int] calcGenoprobF2(Genotype!F2[][] genotypes, double[] rec_frac, double error_prob)
in {
  assert(genotypes[0].length == rec_frac.length+1, "no. markers in genotypes doesn't match rec_frac length");
  assert(error_prob >= 0 && error_prob <= 1, "error_prob out of range");
  foreach(rf; rec_frac) {
    assert(rf >= 0 && rf <= 0.5, "rec_frac must be >= 0 and <= 0.5");
  }
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
    for(auto pos = n_markers-2; pos >= 0; pos--) {
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
  writeln("    unit test calcGenoprobF2:");
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
  auto rec_frac = mapFunction(dist_cM, "haldane");

  writeln("      - Run calcGenoprobF2");
  auto genoprobs = calcGenoprobF2(chr_4_genotypes, rec_frac, 0.002);

  writeln("      - Compare results to R/qtl");
  double[F2][int] genoprobs_from_rqtl;
  /* probs from R/qtl for individual 1 */
  genoprobs_from_rqtl[0] = [F2.A:0.99365258749878116, F2.H:0.005366774350785078, F2.B:0.00098063815043388934];
  genoprobs_from_rqtl[1] = [F2.A:0.01597337476839351, F2.H:0.984010456989931837, F2.B:0.00001616824167473149];
  genoprobs_from_rqtl[2] = [F2.A:0.98819056080900758, F2.H:0.010834274013193156, F2.B:0.00097516517779884938];
  genoprobs_from_rqtl[3] = [F2.A:0.00156449602454221, F2.H:0.998274388090304443, F2.B:0.00016111588515345984];
  
  foreach(i; 0..genoprobs[0].length) {
    foreach(j; [F2.A, F2.H, F2.B]) {
      assert(abs(genoprobs[0][i][j] - genoprobs_from_rqtl[i][j]) < 1e-14);
    }
  }
	
  /* probs from R/qtl for individual 88 */
  genoprobs_from_rqtl[0] = [F2.A:0.0001172399151822272, F2.H:0.0006829310552180061, F2.B:0.999199829029599917];
  genoprobs_from_rqtl[1] = [F2.A:0.0009524442550889619, F2.H:0.0580825863558047800, F2.B:0.940964969389106454];
  genoprobs_from_rqtl[2] = [F2.A:0.0001167398194463191, F2.H:0.0011825266953902206, F2.B:0.998700733485163750];
  genoprobs_from_rqtl[3] = [F2.A:0.0001586426302537041, F2.H:0.9982631735224261060, F2.B:0.001578183847321077];

  foreach(i; 0..genoprobs[87].length) {
    foreach(j; [F2.A, F2.H, F2.B]) {
      assert(abs(genoprobs[87][i][j] - genoprobs_from_rqtl[i][j]) < 1e-14);
    }
  }
	
  
  /* probs from R/qtl for individual 103 */
  genoprobs_from_rqtl[0] = [F2.A:0.3938894726346546, F2.H:0.467329331452057628, F2.B:0.1387811959132879691];
  genoprobs_from_rqtl[1] = [F2.A:0.4722757014939976, F2.H:0.429690530637569956, F2.B:0.0980337678684325559];
  genoprobs_from_rqtl[2] = [F2.A:0.5756148456641348, F2.H:0.365802567955339664, F2.B:0.0585825863805257627];
  genoprobs_from_rqtl[3] = [F2.A:0.9970029970029969, F2.H:0.001998001998001998, F2.B:0.0009990009990009992];

  foreach(i; 0..genoprobs[102].length) {
    foreach(j; [F2.A, F2.H, F2.B]) {
      assert(abs(genoprobs[102][i][j] - genoprobs_from_rqtl[i][j]) < 1e-14);
    }
  }

  /* probs from R/qtl for individual 106 */
  genoprobs_from_rqtl[0] = [F2.A:0.1387811959132879691, F2.H:0.467329331452057628, F2.B:0.3938894726346546];
  genoprobs_from_rqtl[1] = [F2.A:0.0980337678684325975, F2.H:0.429690530637569956, F2.B:0.4722757014939976];
  genoprobs_from_rqtl[2] = [F2.A:0.0585825863805258112, F2.H:0.365802567955339497, F2.B:0.5756148456641349];
  genoprobs_from_rqtl[3] = [F2.A:0.0009990009990009992, F2.H:0.001998001998001998, F2.B:0.9970029970029971];
  
  foreach(i; 0..genoprobs[105].length) {
    foreach(j; [F2.A, F2.H, F2.B]) {
      assert(abs(genoprobs[105][i][j] - genoprobs_from_rqtl[i][j]) < 1e-14);
    }
  }

  /* probs from R/qtl for individual 107 */
  genoprobs_from_rqtl[0] = [F2.A:0.2336319623541089074, F2.H:0.5327360752917820, F2.B:0.2336319623541089074];
  genoprobs_from_rqtl[1] = [F2.A:0.2147748854695731846, F2.H:0.5704502290608535, F2.B:0.2147748854695732956];
  genoprobs_from_rqtl[2] = [F2.A:0.1827669522138612168, F2.H:0.6344660955722773, F2.B:0.1827669522138613833];
  genoprobs_from_rqtl[3] = [F2.A:0.0005005005005005007, F2.H:0.9989989989989990, F2.B:0.0005005005005005007];
  
  foreach(i; 0..genoprobs[106].length) {
    foreach(j; [F2.A, F2.H, F2.B]) {
      assert(abs(genoprobs[106][i][j] - genoprobs_from_rqtl[i][j]) < 1e-14);
    }
  }
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
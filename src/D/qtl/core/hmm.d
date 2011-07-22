/**
 * HMM module
 */

module qtl.core.hmm;

// things I think I really need
import qtl.core.primitives;
import qtl.core.genotype;
import qtl.core.map_functions;
import std.stdio;
import std.math;
import qtl.core.hmm_f2;
import qtl.core.hmm_bc;

// things for the unit tests 
import qtl.plugins.input.read_csv;
import std.path;


// calculate QTL genotype probabilities
double[F2][int][int] calcGenoprobF2(Genotype!F2[][] genotypes, double[] rec_frac, double error_prob)
{
  if(genotypes[0].length != rec_frac.length+1)
    throw new Exception("no. markers in genotypes doesn't match rec_frac length");
  if(error_prob < 0.0 || error_prob > 1.0)
    throw new Exception("error_prob out of range");
  foreach(rf; rec_frac) {
    if(rf < 0 || rf > 0.5)
      throw new Exception("rec_frac must be >= 0 and <= 0.5");
  }

  int n_individuals = genotypes.length;
  int n_markers = genotypes[0].length;
  F2[] all_true_geno = [F2.A, F2.H, F2.B];

  double[int][F2] alpha, beta;
  double[F2][int][int] genoprobs;

  foreach(ind; 0..n_individuals) {
    alpha = forwardEquationsF2(genotypes[ind], all_true_geno, rec_frac, error_prob);
	
    beta = backwardEquationsF2(genotypes[ind], all_true_geno, rec_frac, error_prob);

    // calculate genotype probabilities
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


// forward Equations for F2
double[int][F2] forwardEquationsF2(Genotype!F2[] genotypes, F2[] all_true_geno, 
				   double[] rec_frac, double error_prob) 
{
  int n_markers = genotypes.length;

  double[int][F2] alpha;

  // initialize alphas
  foreach(true_geno; all_true_geno) {
    alpha[true_geno][0] = initF2(true_geno) + emitF2(genotypes[0], true_geno, error_prob);
  }

  foreach(pos; 1 .. n_markers) {
    foreach(true_geno_right; all_true_geno) {

      alpha[true_geno_right][pos] = alpha[all_true_geno[0]][pos-1] + 
	stepF2(all_true_geno[0], true_geno_right, rec_frac[pos-1]);

      foreach(true_geno_left; all_true_geno[1..$]) {
	alpha[true_geno_right][pos] = addlog(alpha[true_geno_right][pos], 
					     alpha[true_geno_left][pos-1] + 
					     stepF2(true_geno_left, true_geno_right, rec_frac[pos-1]));
      }
      alpha[true_geno_right][pos] += emitF2(genotypes[pos], true_geno_right, error_prob);
    }
  }
  return alpha;
}


// backward Equations for F2
double[int][F2] backwardEquationsF2(Genotype!F2[] genotypes, F2[] all_true_geno, 
				    double[] rec_frac, double error_prob) 
{
  int n_markers = genotypes.length;

  double[int][F2] beta;

  // initialize beta
  foreach(true_geno; all_true_geno) {
    beta[true_geno][n_markers-1] = 0.0;
  }

  // backward equations
  for(auto pos = n_markers-2; pos >= 0; pos--) {
    foreach(true_geno_left; all_true_geno) {
      beta[true_geno_left][pos] = beta[all_true_geno[0]][pos+1] + 
	stepF2(true_geno_left, all_true_geno[0], rec_frac[pos]) + 
	emitF2(genotypes[pos+1], all_true_geno[0], error_prob);

      foreach(true_geno_right; all_true_geno[1..$]) {
	beta[true_geno_left][pos] = addlog(beta[true_geno_left][pos], 
					   beta[true_geno_right][pos+1] + 
					   stepF2(true_geno_left, true_geno_right, rec_frac[pos])+
					   emitF2(genotypes[pos+1], true_geno_right, error_prob));
      }
    }
  }

  return beta;
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


// calculate QTL genotype probabilities
double[BC][int][int] calcGenoprobBC(Genotype!BC[][] genotypes, double[] rec_frac, double error_prob)
{
  if(genotypes[0].length != rec_frac.length+1)
    throw new Exception("no. markers in genotypes doesn't match rec_frac length");
  if(error_prob < 0.0 || error_prob > 1.0)
    throw new Exception("error_prob out of range");
  foreach(rf; rec_frac) {
    if(rf < 0 || rf > 0.5)
      throw new Exception("rec_frac must be >= 0 and <= 0.5");
  }

  int n_individuals = genotypes.length;
  int n_markers = genotypes[0].length;
  BC[] all_true_geno = [BC.A, BC.H];

  double[int][BC] alpha, beta;
  double[BC][int][int] genoprobs;

  foreach(ind; 0..n_individuals) {
    alpha = forwardEquationsBC(genotypes[ind], all_true_geno, rec_frac, error_prob);
	
    beta = backwardEquationsBC(genotypes[ind], all_true_geno, rec_frac, error_prob);

    // calculate genotype probabilities
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


// forward Equations for BC
double[int][BC] forwardEquationsBC(Genotype!BC[] genotypes, BC[] all_true_geno, 
				   double[] rec_frac, double error_prob) 
{
  int n_markers = genotypes.length;

  double[int][BC] alpha;

  // initialize alphas
  foreach(true_geno; all_true_geno) {
    alpha[true_geno][0] = initBC(true_geno) + emitBC(genotypes[0], true_geno, error_prob);
  }

  foreach(pos; 1 .. n_markers) {
    foreach(true_geno_right; all_true_geno) {

      alpha[true_geno_right][pos] = alpha[all_true_geno[0]][pos-1] + 
	stepBC(all_true_geno[0], true_geno_right, rec_frac[pos-1]);

      foreach(true_geno_left; all_true_geno[1..$]) {
	alpha[true_geno_right][pos] = addlog(alpha[true_geno_right][pos], 
					     alpha[true_geno_left][pos-1] + 
					     stepBC(true_geno_left, true_geno_right, rec_frac[pos-1]));
      }
      alpha[true_geno_right][pos] += emitBC(genotypes[pos], true_geno_right, error_prob);
    }
  }
  return alpha;
}


// backward Equations for BC
double[int][BC] backwardEquationsBC(Genotype!BC[] genotypes, BC[] all_true_geno, 
				    double[] rec_frac, double error_prob) 
{
  int n_markers = genotypes.length;

  double[int][BC] beta;

  // initialize beta
  foreach(true_geno; all_true_geno) {
    beta[true_geno][n_markers-1] = 0.0;
  }

  // backward equations
  for(auto pos = n_markers-2; pos >= 0; pos--) {
    foreach(true_geno_left; all_true_geno) {
      beta[true_geno_left][pos] = beta[all_true_geno[0]][pos+1] + 
	stepBC(true_geno_left, all_true_geno[0], rec_frac[pos]) + 
	emitBC(genotypes[pos+1], all_true_geno[0], error_prob);

      foreach(true_geno_right; all_true_geno[1..$]) {
	beta[true_geno_left][pos] = addlog(beta[true_geno_left][pos], 
					   beta[true_geno_right][pos+1] + 
					   stepBC(true_geno_left, true_geno_right, rec_frac[pos])+
					   emitBC(genotypes[pos+1], true_geno_right, error_prob));
      }
    }
  }

  return beta;
}



unittest {
  writeln("Unit test " ~ __FILE__);
  writeln("    unit test calcGenoprobBC:");
  alias std.path.join join;
  auto fn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","test","data","input","hyper_noX.csv");
  writeln("      - read CSV " ~ fn);
  auto data = new ReadSimpleCSV!BC(fn);
  
  Marker[] markers_on_chr_5;
  writeln("    Grab markers");
  foreach(marker; data.markers) {
    if(marker.chromosome.name=="5") {
      markers_on_chr_5 ~= marker;
    }
  }

  writeln("      - Subset genotype data");
  Genotype!BC[][] chr_5_genotypes;
  chr_5_genotypes.reserve(data.genotypes.length);
  foreach(i; 0..data.genotypes.length) {
    Genotype!BC[] an_individuals_genotype;
    an_individuals_genotype.reserve(markers_on_chr_5.length);
    foreach(j; 0..markers_on_chr_5.length) {
      an_individuals_genotype ~= data.genotypes[i][markers_on_chr_5[j].id];
    }
    chr_5_genotypes ~= an_individuals_genotype;
  }

  writeln("      - Get recombination fractions");
  double[] dist_cM;
  foreach(i; 1..markers_on_chr_5.length) {
    dist_cM ~= markers_on_chr_5[i].position - markers_on_chr_5[i-1].position;
  }
  auto rec_frac = mapFunction(dist_cM, "carter-falconer");

  double[] rec_frac_from_rqtl = [0.05499838959578960, 0.05399853076176147, 0.03299987476879847,
				 0.01099999948567994, 0.03299987476879847, 0.14181577930458458, 
				 0.18529361912102324, 0.08698405702618413, 0.06599599307740488,
				 0.06499628754224140, 0.01099999948567990, 0.03299987476879847, 
				 0.04399947228314834]; 

  foreach(i; 0..rec_frac.length) {
    assert(abs(rec_frac[i] - rec_frac_from_rqtl[i]) < 1e-12);
  }



  writeln("      - Run calcGenoprobBC");
  auto genoprobs = calcGenoprobBC(chr_5_genotypes, rec_frac, 0.002);

  writeln("      - Compare results to R/qtl");
  double[BC][int] genoprobs_from_rqtl;
  /* probs from R/qtl for individual 1 */
  genoprobs_from_rqtl[0] =  [BC.A:0.999769660491391820578, BC.H:0.0002303395086083074390];
  genoprobs_from_rqtl[1] =  [BC.A:0.996571637787474373660, BC.H:0.0034283622125256579989];
  genoprobs_from_rqtl[2] =  [BC.A:0.999992001840160571469, BC.H:0.0000079981598392023363];
  genoprobs_from_rqtl[3] =  [BC.A:0.999999229385373444856, BC.H:0.0000007706146264155516];
  genoprobs_from_rqtl[4] =  [BC.A:0.999999226845584399115, BC.H:0.0000007731544155547456];
  genoprobs_from_rqtl[5] =  [BC.A:0.999988106195807358034, BC.H:0.0000118938041927743214];
  genoprobs_from_rqtl[6] =  [BC.A:0.998546612868447480693, BC.H:0.0014533871315524147898];
  genoprobs_from_rqtl[7] =  [BC.A:0.000838837664664903023, BC.H:0.9991611623353346960386];
  genoprobs_from_rqtl[8] =  [BC.A:0.000014752882287688144, BC.H:0.9999852471177120838419];
  genoprobs_from_rqtl[9] =  [BC.A:0.000009884987083730629, BC.H:0.9999901150129172355818];
  genoprobs_from_rqtl[10] = [BC.A:0.000006327226189669525, BC.H:0.9999936727738113484421];
  genoprobs_from_rqtl[11] = [BC.A:0.000385258672705765363, BC.H:0.9996147413272942205964];
  genoprobs_from_rqtl[12] = [BC.A:0.000004366470355155772, BC.H:0.9999956335296438236782];
  genoprobs_from_rqtl[13] = [BC.A:0.000092406805046532986, BC.H:0.9999075931949532591858];
  
  foreach(i; 0..genoprobs[0].length) {
    foreach(j; [BC.A, BC.H]) {
      assert(abs(genoprobs[0][i][j] - genoprobs_from_rqtl[i][j]) < 1e-6);
    }
  }
	
  /* probs from R/qtl for individual 88 */
  genoprobs_from_rqtl[0] =  [BC.A:0.9667756240797861222447, BC.H:0.033224375920213690405];
  genoprobs_from_rqtl[1] =  [BC.A:0.0018989493588511250667, BC.H:0.998101050641148845877];
  genoprobs_from_rqtl[2] =  [BC.A:0.0000061727288390725428, BC.H:0.999993827271161594972];
  genoprobs_from_rqtl[3] =  [BC.A:0.0000007694148986088805, BC.H:0.999999230585101317459];
  genoprobs_from_rqtl[4] =  [BC.A:0.0000007727830675603257, BC.H:0.999999227216932462525];
  genoprobs_from_rqtl[5] =  [BC.A:0.0000113389799089367530, BC.H:0.999988661020090652265];
  genoprobs_from_rqtl[6] =  [BC.A:0.0000754209410127622390, BC.H:0.999924579058987128555];
  genoprobs_from_rqtl[7] =  [BC.A:0.0000562473769980956022, BC.H:0.999943752623000836621];
  genoprobs_from_rqtl[8] =  [BC.A:0.0027022235245518625747, BC.H:0.997297776475448261024];
  genoprobs_from_rqtl[9] =  [BC.A:0.9980350818191259243406, BC.H:0.001964918180874547071];
  genoprobs_from_rqtl[10] = [BC.A:0.9999911337652792608921, BC.H:0.000008866234720407318];
  genoprobs_from_rqtl[11] = [BC.A:0.9996128260665828602072, BC.H:0.000387173933416828464];
  genoprobs_from_rqtl[12] = [BC.A:0.9999956281497925925095, BC.H:0.000004371850207192675];
  genoprobs_from_rqtl[13] = [BC.A:0.9999075929709730914396, BC.H:0.000092407029026649721];

  foreach(i; 0..genoprobs[87].length) {
    foreach(j; [BC.A, BC.H]) {
      assert(abs(genoprobs[87][i][j] - genoprobs_from_rqtl[i][j]) < 1e-6);
    }
  }
	
  
  /* probs from R/qtl for individual 103 */
  genoprobs_from_rqtl[0] =  [BC.A:0.0002366636922487103, BC.H:0.99976333630775027217];
  genoprobs_from_rqtl[1] =  [BC.A:0.0036190176590203296, BC.H:0.99638098234097982608];
  genoprobs_from_rqtl[2] =  [BC.A:0.0003819364312548427, BC.H:0.99961806356874582935];
  genoprobs_from_rqtl[3] =  [BC.A:0.0549223423890561752, BC.H:0.94507765761094364443];
  genoprobs_from_rqtl[4] =  [BC.A:0.0722315882263140224, BC.H:0.92776841177368585267];
  genoprobs_from_rqtl[5] =  [BC.A:0.1240732453345061992, BC.H:0.87592675466549385632];
  genoprobs_from_rqtl[6] =  [BC.A:0.3563901310091487917, BC.H:0.64360986899085126378];
  genoprobs_from_rqtl[7] =  [BC.A:0.6573119740963657698, BC.H:0.34268802590363423022];
  genoprobs_from_rqtl[8] =  [BC.A:0.7877017986121217508, BC.H:0.21229820138787824924];
  genoprobs_from_rqtl[9] =  [BC.A:0.8907905734005259202, BC.H:0.10920942659947401043];
  genoprobs_from_rqtl[10] = [BC.A:0.9998456411901881502, BC.H:0.00015435880981190807];
  genoprobs_from_rqtl[11] = [BC.A:0.9994803764517402600, BC.H:0.00051962354825971387];
  genoprobs_from_rqtl[12] = [BC.A:0.9999028296560041884, BC.H:0.00009717034399581742];
  genoprobs_from_rqtl[13] = [BC.A:0.9559119082605704865, BC.H:0.04408809173942938170];

  foreach(i; 0..genoprobs[102].length) {
    foreach(j; [BC.A, BC.H]) {
      assert(abs(genoprobs[102][i][j] - genoprobs_from_rqtl[i][j]) < 1e-6);
    }
  }

  /* probs from R/qtl for individual 106 */
  genoprobs_from_rqtl[0] =  [BC.A:0.00023263214083326644, BC.H:0.9997673678591675];
  genoprobs_from_rqtl[1] =  [BC.A:0.00349747829854584087, BC.H:0.9965025217014539];
  genoprobs_from_rqtl[2] =  [BC.A:0.00014355762173019074, BC.H:0.9998564423782705];
  genoprobs_from_rqtl[3] =  [BC.A:0.01991084477676950681, BC.H:0.9800891552232307];
  genoprobs_from_rqtl[4] =  [BC.A:0.02586339947739364070, BC.H:0.9741366005226065];
  genoprobs_from_rqtl[5] =  [BC.A:0.04268032029532820015, BC.H:0.9573196797046718];
  genoprobs_from_rqtl[6] =  [BC.A:0.09532197837358877268, BC.H:0.9046780216264112];
  genoprobs_from_rqtl[7] =  [BC.A:0.09341117498753083448, BC.H:0.9065888250124690];
  genoprobs_from_rqtl[8] =  [BC.A:0.06741531961276907292, BC.H:0.9325846803872311];
  genoprobs_from_rqtl[9] =  [BC.A:0.03801966338345039859, BC.H:0.9619803366165495];
  genoprobs_from_rqtl[10] = [BC.A:0.00005802675751773932, BC.H:0.9999419732424817];
  genoprobs_from_rqtl[11] = [BC.A:0.00044600617494278827, BC.H:0.9995539938250571];
  genoprobs_from_rqtl[12] = [BC.A:0.00009309207758853352, BC.H:0.9999069079224117];
  genoprobs_from_rqtl[13] = [BC.A:0.04408437235616159688, BC.H:0.9559156276438385];
  
  foreach(i; 0..genoprobs[105].length) {
    foreach(j; [BC.A, BC.H]) {
      assert(abs(genoprobs[105][i][j] - genoprobs_from_rqtl[i][j]) < 1e-6);
    }
  }

  /* probs from R/qtl for individual 107 */
  genoprobs_from_rqtl[0] =  [BC.A:0.00023263214083326644, BC.H:0.9997673678591675];
  genoprobs_from_rqtl[1] =  [BC.A:0.00349747829854584087, BC.H:0.9965025217014539];
  genoprobs_from_rqtl[2] =  [BC.A:0.00014355762173019074, BC.H:0.9998564423782705];
  genoprobs_from_rqtl[3] =  [BC.A:0.01991084477676950681, BC.H:0.9800891552232307];
  genoprobs_from_rqtl[4] =  [BC.A:0.02586339947739364070, BC.H:0.9741366005226065];
  genoprobs_from_rqtl[5] =  [BC.A:0.04268032029532820015, BC.H:0.9573196797046718];
  genoprobs_from_rqtl[6] =  [BC.A:0.09532197837358877268, BC.H:0.9046780216264112];
  genoprobs_from_rqtl[7] =  [BC.A:0.09341117498753083448, BC.H:0.9065888250124690];
  genoprobs_from_rqtl[8] =  [BC.A:0.06741531961276907292, BC.H:0.9325846803872311];
  genoprobs_from_rqtl[9] =  [BC.A:0.03801966338345039859, BC.H:0.9619803366165495];
  genoprobs_from_rqtl[10] = [BC.A:0.00005802675751773932, BC.H:0.9999419732424817];
  genoprobs_from_rqtl[11] = [BC.A:0.00044600617494278827, BC.H:0.9995539938250571];
  genoprobs_from_rqtl[12] = [BC.A:0.00009309207758853352, BC.H:0.9999069079224117];
  genoprobs_from_rqtl[13] = [BC.A:0.04408437235616159688, BC.H:0.9559156276438385];
  
  foreach(i; 0..genoprobs[106].length) {
    foreach(j; [BC.A, BC.H]) {
      assert(abs(genoprobs[106][i][j] - genoprobs_from_rqtl[i][j]) < 1e-6);
    }
  }
}

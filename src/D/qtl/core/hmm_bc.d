/**
 * hmm_bc
 **/

module qtl.core.hmm_bc;

import qtl.core.genotype;
import qtl.core.primitives;
import qtl.plugins.input.read_csv;
import std.math, std.stdio, std.path;
import std.exception;
import qtl.core.hmm_calcgenoprob;
import qtl.core.hmm_util;
import qtl.core.hmm_estmap;
import qtl.core.map_functions;

// things for the unit tests 
import qtl.plugins.input.read_csv;
import std.path;


// marginal genotype probability
double init(BC true_gen)
{
  if(true_gen == BC.A || true_gen == BC.H)
    return(-LN2);

  throw new Exception("true_gen not among the possible true genotypes");
}

unittest {
  writeln("Unit test " ~ __FILE__);
  writeln("    unit test init for BC");
  auto gen = [BC.A, BC.H];
  assert(abs(init(gen[0]) - log(0.5)) < 1e-12);
  assert(abs(init(gen[1]) - log(0.5)) < 1e-12);
}

// emission probability (marker genotype "penetrance")
double emit(Genotype!BC obs_gen, BC true_gen, double error_prob) 
{
  if(error_prob < 0.0 || error_prob > 1.0)
     throw new Exception("error_prob must be >= 0 and <= 1");
  
  switch(obs_gen.value) {
  case BC.NA: return(0.0);
  case BC.A: case BC.H:
    if(obs_gen.value == true_gen) {
      return(log(1.0-error_prob));
    } else {
      return(log(error_prob));
    }
  default: 
    throw new Exception("obs_gen.value not among the possible marker genotypes");
  }
}

unittest {
  writeln("    unit test emit for BC");
  auto gen = [BC.A, BC.H];
  auto ogen = [set_genotype!BC("-"), set_genotype!BC("A"), set_genotype!BC("H")];
  double err = 0.001;
  
  foreach(g; gen) {
    assert(abs( emit(ogen[0], g, err) ) < 1e-12);
  }
  assert(abs( emit(ogen[1], gen[0], err) - log(1-err) ) < 1e-12);
  assert(abs( emit(ogen[1], gen[1], err) - log(err) ) < 1e-12);
  assert(abs( emit(ogen[2], gen[0], err) - log(err) ) < 1e-12);
  assert(abs( emit(ogen[2], gen[1], err) - log(1-err) ) < 1e-12);
}

// transition probabilities 
double step(BC true_gen_left, BC true_gen_right, 
	      double rec_frac)
{
  if(true_gen_left != BC.A && true_gen_left != BC.H)
    throw new Exception("true_gen_left not among the possible true genotypes");
  if(true_gen_right != BC.A && true_gen_right != BC.H)
    throw new Exception("true_gen_right not among the possible true genotypes");
  if(rec_frac < 0.0 || rec_frac > 0.5)
     throw new Exception("rec_frac must be >= 0 and <= 1");

  if(true_gen_left == true_gen_right) {
    return(log(1.0-rec_frac));
  }
  else {
    return(log(rec_frac));
  }
}


unittest {
  writeln("    unit test step for BC");
  auto gen = [BC.A, BC.H];
  double rf = 0.01;
  assert(abs( step(gen[0], gen[0], rf) - log(1-rf) ) < 1e-12);
  assert(abs( step(gen[0], gen[1], rf) - log(rf) ) < 1e-12);
  assert(abs( step(gen[1], gen[0], rf) - log(rf) ) < 1e-12);
  assert(abs( step(gen[1], gen[1], rf) - log(1-rf) ) < 1e-12);
}

  
double nrec(BC true_gen_left, BC true_gen_right)
{
  if(true_gen_left != BC.A && true_gen_left != BC.H) 
    throw new Exception("true_gen_left not a valid genotype.");
  if(true_gen_right != BC.A && true_gen_right != BC.H) 
    throw new Exception("true_gen_right not a valid genotype.");

  if(true_gen_left == true_gen_right) return(0.0);
  else return(1.0);
}

unittest {
  writeln("    unit test nrec for BC");
  auto gen = [BC.A, BC.H];
  assert(nrec(gen[0], gen[0]) == 0.0);
  assert(nrec(gen[0], gen[1]) == 1.0);
  assert(nrec(gen[1], gen[0]) == 1.0);
  assert(nrec(gen[1], gen[1]) == 0.0);
}


// return vector of all possible true genotypes
BC[] allTrueGeno(BC one)
{
  return [BC.A, BC.H];
}


// return vector of all possible true genotypes, phase-known case
BC[] allTrueGenoPK(BC one)
{
  return [BC.A, BC.H];
}

unittest {
  writeln("    unit test allTrueGeno for BC");
  assert(allTrueGeno(BC.A) == [BC.A, BC.H]);
  assert(allTrueGenoPK(BC.A) == [BC.A, BC.H]);
}


// calcGenoprob for BC
mixin calcGenoprobCode!BC;

unittest {
  writeln("    unit test calcGenoprob for BC:");
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



  writeln("      - Run calcGenoprob for BC");
  auto genoprobs = calcGenoprob(chr_5_genotypes, rec_frac, 0.002);

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


// estmap for BC
mixin estmapCode!(BC, BC);

unittest {
  writeln("    unit test estmapBC:");
  alias std.path.join join;
  auto fn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","test","data","input","hyper_noX.csv");
  writeln("      - read CSV " ~ fn);
  auto data = new ReadSimpleCSV!BC(fn);
  
  Marker[] markers_on_chr_15;
  writeln("    Grab markers");
  foreach(marker; data.markers) {
    if(marker.chromosome.name=="15") {
      markers_on_chr_15 ~= marker;
    }
  }

  writeln("      - Subset genotype data");
  Genotype!BC[][] chr_15_genotypes;
  chr_15_genotypes.reserve(data.genotypes.length);
  foreach(i; 0..data.genotypes.length) {
    Genotype!BC[] an_individuals_genotype;
    an_individuals_genotype.reserve(markers_on_chr_15.length);
    foreach(j; 0..markers_on_chr_15.length) {
      an_individuals_genotype ~= data.genotypes[i][markers_on_chr_15[j].id];
    }
    chr_15_genotypes ~= an_individuals_genotype;
  }

  writeln("      - Get recombination fractions");
  double[] dist_cM;
  foreach(i; 1..markers_on_chr_15.length) {
    dist_cM ~= markers_on_chr_15[i].position - markers_on_chr_15[i-1].position;
  }
  auto rec_frac = mapFunction(dist_cM, "carter-falconer");

  auto rec_frac_rqtl = [0.0000000000010000000000000001818,
			0.0219999835094050709416446665045,
			0.0539985307617614693209695531095,
			0.0329998747687984747556377840283,
			0.0000000000010000000000000001818,
			0.0109999994856799232501032292930,
			0.0000000000010000000000000001818,
			0.1199204909742600105859722248169,
			0.2581768942841565772639000897470,
			0.0769913404562991288138107393024];

  foreach(i; 0..rec_frac.length) {
    assert(abs(rec_frac[i] - rec_frac_rqtl[i]) < 1e-13);
  }

  auto rec_frac_rev_rqtl = [0.020807743646802073778,
			    0.035733145049573165897,
			    0.037277422953319969134,
			    0.046231824599400450637,
			    0.027836290106153856183,
			    0.098343160668203172259,
			    0.030788500492084123344,
			    0.150350381248621794983,
			    0.162944847048261037825,
			    0.075138233383357094786];

  auto rec_frac_rev = estmap(chr_15_genotypes, rec_frac, 0.001, 100, 1e-6, true);

  foreach(i; 0..rec_frac.length) {
    assert(abs(rec_frac_rev[i] - rec_frac_rev_rqtl[i]) < 1e-5);
  }
}

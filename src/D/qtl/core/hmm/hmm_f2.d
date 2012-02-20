/**
 * hmm_f2
 **/

module qtl.core.hmm.hmm_f2;

import std.math, std.stdio, std.path;
import std.conv;

import qtl.core.deprecate.genotype_enum;
import qtl.core.primitives;
import qtl.plugins.deprecate.read_csv;
import qtl.core.map.genetic_map_functions;
import qtl.core.hmm.hmm_calcgenoprob;
import qtl.core.hmm.hmm_util;
import qtl.core.hmm.hmm_estmap;

// things for the unit tests 
import std.path;


// marginal genotype probability
double init(F2 true_gen)
{
  switch(true_gen) {
  case F2.H:
    return(-LN2);

  case F2.A: case F2.B: 
    return(-2.0*LN2);

  default:
    throw new Exception("true_gen not among the possible true genotypes");
  }
}

unittest {
  writeln("Unit test " ~ __FILE__);
  writeln("    unit test init for F2");
  auto gen = [F2.A, F2.H, F2.B];
  assert(abs(init(gen[0]) - log(0.25)) < 1e-12);
  assert(abs(init(gen[1]) - log(0.5)) < 1e-12);
  assert(abs(init(gen[2]) - log(0.25)) < 1e-12);
}


// emission probability (marker genotype "penetrance")
double emit(Genotype!F2 obs_gen, F2 true_gen, double error_prob) 
{
  if(true_gen != F2.A && true_gen != F2.H && true_gen != F2.B) 
    throw new Exception("true_gen not among the possible true genotypes");
  if(error_prob < 0.0 || error_prob > 1.0)
    throw new Exception("error_prob must be >= 0 and <= 1");
  
  switch(obs_gen.value) {
  case F2.NA: return(0.0);
    
  case F2.A: case F2.H: case F2.B:
    if(obs_gen.value == true_gen) {
      return(log(1.0-error_prob));
    } else {
      return(log(error_prob)-LN2);
    }

  case F2.HorB:
    if(true_gen == F2.A) {
      return(log(error_prob));
    } else {
      return(log(1.0-error_prob/2.0));
    }

  case F2.HorA:
    if(true_gen == F2.B) {
      return(log(error_prob));
    } else {
      return(log(1.0-error_prob/2.0));
    }

  default:
    throw new Exception("obs_gen.value not among possible observed marker genotypes.");
  }
}

unittest {
  writeln("    unit test emit for F2");
  auto gen = [F2.A, F2.H, F2.B];
  auto ogen = [set_genotype!F2("-"), 
	       set_genotype!F2("A"), set_genotype!F2("H"), set_genotype!F2("B"),
	       set_genotype!F2("D"), set_genotype!F2("C")];
  double err = 0.001;

  /* missing */
  foreach(g; gen) {
    assert(abs(emit(ogen[0], g, err) ) < 1e-12);
  }

  /* A */
  assert(abs(emit(ogen[1], gen[0], err) - log(1.0-err)) < 1e-12);
  assert(abs(emit(ogen[1], gen[1], err) - log(err/2.0)) < 1e-12);
  assert(abs(emit(ogen[1], gen[2], err) - log(err/2.0)) < 1e-12);

  /* H */
  assert(abs(emit(ogen[2], gen[0], err) - log(err/2.0)) < 1e-12);
  assert(abs(emit(ogen[2], gen[1], err) - log(1.0-err)) < 1e-12);
  assert(abs(emit(ogen[2], gen[2], err) - log(err/2.0)) < 1e-12);

  /* B */
  assert(abs(emit(ogen[3], gen[0], err) - log(err/2.0)) < 1e-12);
  assert(abs(emit(ogen[3], gen[1], err) - log(err/2.0)) < 1e-12);
  assert(abs(emit(ogen[3], gen[2], err) - log(1.0-err)) < 1e-12);

  /* AorH (aka D) */
  assert(abs(emit(ogen[4], gen[0], err) - log(1.0-err/2.0)) < 1e-12);
  assert(abs(emit(ogen[4], gen[1], err) - log(1.0-err/2.0)) < 1e-12);
  assert(abs(emit(ogen[4], gen[2], err) - log(err)) < 1e-12);

  /* HorB (aka C) */
  assert(abs(emit(ogen[5], gen[0], err) - log(err)) < 1e-12);
  assert(abs(emit(ogen[5], gen[1], err) - log(1.0-err/2.0)) < 1e-12);
  assert(abs(emit(ogen[5], gen[2], err) - log(1.0-err/2.0)) < 1e-12);
}


// transition probabilities 
double step(in F2 true_gen_left, in F2 true_gen_right, in double rec_frac) 
{
  if(true_gen_right != F2.A && true_gen_right != F2.H &&
     true_gen_right != F2.B) 
    throw new Exception("true_gen_right not among the possible true genotypes");
  if(rec_frac < 0.0 || rec_frac > 0.5)
    throw new Exception("rec_frac must be >= 0 and <= 0.5");

  switch(true_gen_left) {
    case F2.A:
      switch(true_gen_right) {
      case F2.A: return(2.0*log(1.0-rec_frac));
      case F2.H: return(LN2 + log(1.0-rec_frac) + log(rec_frac));
      case F2.B: return(2.0*log(rec_frac));
      default: break;
      }

    case F2.H:
      switch(true_gen_right) {
      case F2.A: case F2.B: return(log(rec_frac) + log(1.0-rec_frac));
      case F2.H: return(log((1.0-rec_frac)^^2 + rec_frac^^2));
      default: break;
      }

    case F2.B:
      switch(true_gen_right) {
      case F2.A: return(2.0*log(rec_frac));
      case F2.H: return(LN2 + log(1.0-rec_frac) + log(rec_frac));
      case F2.B: return(2.0*log(1.0-rec_frac));
      default: break;
      }

    default: 
      throw new Exception("true_gen_left not among the possible true genotypes");
  }
}

unittest {
  writeln("    unit test step for F2");
  auto gen = [F2.A, F2.H, F2.B];
  double rf = 0.01;

  assert( abs( step(gen[0], gen[0], rf) - log((1-rf)^^2) ) < 1e-12);
  assert( abs( step(gen[0], gen[1], rf) - log(2.0*rf*(1-rf)) ) < 1e-12);
  assert( abs( step(gen[0], gen[2], rf) - log(rf^^2) ) < 1e-12);

  assert( abs( step(gen[1], gen[0], rf) - log(rf*(1-rf)) ) < 1e-12);
  assert( abs( step(gen[1], gen[1], rf) - log((1-rf)^^2 + rf^^2) ) < 1e-12);
  assert( abs( step(gen[1], gen[2], rf) - log(rf*(1-rf)) ) < 1e-12);

  assert( abs( step(gen[2], gen[2], rf) - log((1-rf)^^2) ) < 1e-12);
  assert( abs( step(gen[2], gen[1], rf) - log(2.0*rf*(1-rf)) ) < 1e-12);
  assert( abs( step(gen[2], gen[0], rf) - log(rf^^2) ) < 1e-12);
}

  
double init(F2pk true_gen)
{
  if(true_gen != F2pk.AA && true_gen != F2pk.AB &&
     true_gen != F2pk.BA && true_gen != F2pk.BB)
    throw new Exception("true_gen not a valid genotype.");

  return(-2.0*LN2);  /* ln(0.25) */
}

unittest {
  writeln("    unit test init for F2 phase-known");
  auto gen = [F2pk.AA, F2pk.AB, F2pk.BA, F2pk.BB];
  foreach(g; gen) {
    assert(abs(init(g) - log(0.25)) < 1e-12);
  }
}




double emit(Genotype!F2 obs_gen, F2pk true_gen, double error_prob)
{
  if(true_gen != F2pk.AA && true_gen != F2pk.AB && 
     true_gen != F2pk.BA && true_gen != F2pk.BB)
    throw new Exception("true_gen not among the possible true genotypes");
  if(error_prob < 0.0 || error_prob > 1.0)
    throw new Exception("error_prob must be >= 0 and <= 1");
  
  switch(obs_gen.value) {
  case F2.NA: return(0.0);
    
  case F2.A: 
    if(true_gen==F2pk.AA) {
      return(log(1.0-error_prob));
    } else {
      return(log(error_prob)-LN2);
    }

  case F2.H:
    if(true_gen==F2pk.AB || true_gen==F2pk.BA) {
      return(log(1.0-error_prob));
    } else {
      return(log(error_prob)-LN2);
    }

  case F2.B:
    if(true_gen==F2pk.BB) {
      return(log(1.0-error_prob));
    } else {
      return(log(error_prob)-LN2);
    }

  case F2.HorB:
    if(true_gen == F2pk.AA) {
      return(log(error_prob));
    } else {
      return(log(1.0-error_prob/2.0));
    }

  case F2.HorA:
    if(true_gen == F2pk.BB) {
      return(log(error_prob));
    } else {
      return(log(1.0-error_prob/2.0));
    }

  default:
    throw new Exception("obs_gen.value not among possible observed marker genotypes.");
  }
}
    
  
unittest {
  writeln("    unit test emit for F2 phase-known");

  auto gen = [F2pk.AA, F2pk.AB, F2pk.BA, F2pk.BB];
  auto ogen = [set_genotype!F2("-"), 
	       set_genotype!F2("A"), set_genotype!F2("H"), set_genotype!F2("B"),
	       set_genotype!F2("D"), set_genotype!F2("C")];
  double err = 0.001;

  /* missing */
  foreach(g; gen) {
    assert(abs(emit(ogen[0], g, err) ) < 1e-12);
  }

  /* A */
  assert(abs(emit(ogen[1], gen[0], err) - log(1.0-err)) < 1e-12);
  assert(abs(emit(ogen[1], gen[1], err) - log(err/2.0)) < 1e-12);
  assert(abs(emit(ogen[1], gen[2], err) - log(err/2.0)) < 1e-12);
  assert(abs(emit(ogen[1], gen[3], err) - log(err/2.0)) < 1e-12);

  /* H */
  assert(abs(emit(ogen[2], gen[0], err) - log(err/2.0)) < 1e-12);
  assert(abs(emit(ogen[2], gen[1], err) - log(1.0-err)) < 1e-12);
  assert(abs(emit(ogen[2], gen[2], err) - log(1.0-err)) < 1e-12);
  assert(abs(emit(ogen[2], gen[3], err) - log(err/2.0)) < 1e-12);

  /* B */
  assert(abs(emit(ogen[3], gen[0], err) - log(err/2.0)) < 1e-12);
  assert(abs(emit(ogen[3], gen[1], err) - log(err/2.0)) < 1e-12);
  assert(abs(emit(ogen[3], gen[2], err) - log(err/2.0)) < 1e-12);
  assert(abs(emit(ogen[3], gen[3], err) - log(1.0-err)) < 1e-12);

  /* AorH (aka D) */
  assert(abs(emit(ogen[4], gen[0], err) - log(1.0-err/2.0)) < 1e-12);
  assert(abs(emit(ogen[4], gen[1], err) - log(1.0-err/2.0)) < 1e-12);
  assert(abs(emit(ogen[4], gen[2], err) - log(1.0-err/2.0)) < 1e-12);
  assert(abs(emit(ogen[4], gen[3], err) - log(err)) < 1e-12);

  /* HorB (aka C) */
  assert(abs(emit(ogen[5], gen[0], err) - log(err)) < 1e-12);
  assert(abs(emit(ogen[5], gen[1], err) - log(1.0-err/2.0)) < 1e-12);
  assert(abs(emit(ogen[5], gen[2], err) - log(1.0-err/2.0)) < 1e-12);
  assert(abs(emit(ogen[5], gen[3], err) - log(1.0-err/2.0)) < 1e-12);
}


double step(F2pk true_gen_left, F2pk true_gen_right, double rec_frac)
{
  if(true_gen_right != F2pk.AA && true_gen_right != F2pk.AB &&
     true_gen_right != F2pk.BA && true_gen_right != F2pk.BB)
    throw new Exception("true_gen_right not among the possible true genotypes");
  if(rec_frac < 0.0 || rec_frac > 0.5)
    throw new Exception("rec_frac must be >= 0 and <= 0.5");

  switch(true_gen_left) {
    case F2pk.AA:
      switch(true_gen_right) {
      case F2pk.AA: return(2.0*log(1.0-rec_frac));
      case F2pk.AB: case F2pk.BA: return(log(1.0-rec_frac) + log(rec_frac));
      case F2pk.BB: return(2.0*log(rec_frac));
      default: throw new Exception("true_gen_right not among possible observed marker genotypes.");
      }

    case F2pk.AB:
      switch(true_gen_right) {
      case F2pk.AA: case F2pk.BB: return(log(rec_frac) + log(1.0-rec_frac));
      case F2pk.AB: return(2.0*log(1.0-rec_frac));
      case F2pk.BA: return(2.0*log(rec_frac));
      default: throw new Exception("true_gen_right not among possible observed marker genotypes.");
      }

    case F2pk.BA:
      switch(true_gen_right) {
      case F2pk.AA: case F2pk.BB: return(log(rec_frac) + log(1.0-rec_frac));
      case F2pk.BA: return(2.0*log(1.0-rec_frac));
      case F2pk.AB: return(2.0*log(rec_frac));
      default: throw new Exception("true_gen_right not among possible observed marker genotypes.");
      }

    case F2pk.BB:
      switch(true_gen_right) {
      case F2pk.BB: return(2.0*log(1.0-rec_frac));
      case F2pk.AB: case F2pk.BA: return(log(1.0-rec_frac) + log(rec_frac));
      case F2pk.AA: return(2.0*log(rec_frac));
      default: throw new Exception("true_gen_right not among possible observed marker genotypes.");
      }

    default: 
      throw new Exception("true_gen_left not among the possible true genotypes");
  }
}

unittest {
  writeln("    unit test step for F2 phase-known");
  auto gen = [F2pk.AA, F2pk.AB, F2pk.BA, F2pk.BB];
  double rf = 0.01;

  assert( abs( step(gen[0], gen[0], rf) - log((1-rf)^^2) ) < 1e-12);
  assert( abs( step(gen[0], gen[1], rf) - log(rf*(1-rf)) ) < 1e-12);
  assert( abs( step(gen[0], gen[2], rf) - log(rf*(1-rf)) ) < 1e-12);
  assert( abs( step(gen[0], gen[3], rf) - log(rf^^2)     ) < 1e-12);

  assert( abs( step(gen[1], gen[0], rf) - log(rf*(1-rf)) ) < 1e-12);
  assert( abs( step(gen[1], gen[1], rf) - log((1-rf)^^2) ) < 1e-12);
  assert( abs( step(gen[1], gen[2], rf) - log(rf^^2)     ) < 1e-12);
  assert( abs( step(gen[1], gen[3], rf) - log(rf*(1-rf)) ) < 1e-12);

  assert( abs( step(gen[2], gen[3], rf) - log(rf*(1-rf)) ) < 1e-12);
  assert( abs( step(gen[2], gen[2], rf) - log((1-rf)^^2) ) < 1e-12);
  assert( abs( step(gen[2], gen[1], rf) - log(rf^^2)     ) < 1e-12);
  assert( abs( step(gen[2], gen[0], rf) - log(rf*(1-rf)) ) < 1e-12);

  assert( abs( step(gen[3], gen[3], rf) - log((1-rf)^^2) ) < 1e-12);
  assert( abs( step(gen[3], gen[2], rf) - log(rf*(1-rf)) ) < 1e-12);
  assert( abs( step(gen[3], gen[1], rf) - log(rf*(1-rf)) ) < 1e-12);
  assert( abs( step(gen[3], gen[0], rf) - log(rf^^2)     ) < 1e-12);
}



double nrec(F2pk true_gen_left, F2pk true_gen_right)
{
  if(true_gen_right != F2pk.AA && true_gen_right != F2pk.AB &&
     true_gen_right != F2pk.BA && true_gen_right != F2pk.BB)
    throw new Exception("true_gen_right not among the possible true genotypes");

  switch(true_gen_left) {
  case F2pk.AA:
    switch(true_gen_right) {
    case F2pk.AA: return(0.0);
    case F2pk.AB: case F2pk.BA: return(0.5);
    case F2pk.BB: return(1.0);
    default: throw new Exception("true_gen_right not among possible observed marker genotypes.");
    }
  case F2pk.AB:
    switch(true_gen_right) {
    case F2pk.AA: case F2pk.BB: return(0.5);
    case F2pk.AB: return(0.0);
    case F2pk.BA: return(1.0);
    default: throw new Exception("true_gen_right not among possible observed marker genotypes.");
    }
  case F2pk.BA:
    switch(true_gen_right) {
    case F2pk.AA: case F2pk.BB: return(0.5);
    case F2pk.BA: return(0.0);
    case F2pk.AB: return(1.0);
    default: throw new Exception("true_gen_right not among possible observed marker genotypes.");
    }
  case F2pk.BB:
    switch(true_gen_right) {
    case F2pk.BB: return(0.0);
    case F2pk.AB: case F2pk.BA: return(0.5);
    case F2pk.AA: return(1.0);
    default: throw new Exception("true_gen_right not among possible observed marker genotypes.");
    }
  default: 
    throw new Exception("true_gen_left not among the possible true genotypes");
  }

}

unittest {
  writeln("    unit test nrec for F2 phase-known");
  auto gen = [F2pk.AA, F2pk.AB, F2pk.BA, F2pk.BB];

  assert(nrec(gen[0], gen[0])==0.0);
  assert(nrec(gen[0], gen[1])==0.5);
  assert(nrec(gen[0], gen[2])==0.5);
  assert(nrec(gen[0], gen[3])==1.0);

  assert(nrec(gen[1], gen[0])==0.5);
  assert(nrec(gen[1], gen[1])==0.0);
  assert(nrec(gen[1], gen[2])==1.0);
  assert(nrec(gen[1], gen[3])==0.5);

  assert(nrec(gen[2], gen[3])==0.5);
  assert(nrec(gen[2], gen[2])==0.0);
  assert(nrec(gen[2], gen[1])==1.0);
  assert(nrec(gen[2], gen[0])==0.5);

  assert(nrec(gen[3], gen[3])==0.0);
  assert(nrec(gen[3], gen[2])==0.5);
  assert(nrec(gen[3], gen[1])==0.5);
  assert(nrec(gen[3], gen[0])==1.0);
}


// return vector of all possible true genotypes
F2[] allTrueGeno(F2 one)
{
  return [F2.A, F2.H, F2.B];
}


F2pk[] allTrueGenoPK(F2 one)
{
  return [F2pk.AA, F2pk.AB, F2pk.BA, F2pk.BB];
}

unittest {
  writeln("    unit test allTrueGeno for F2");
  assert(allTrueGeno(F2.A) == [F2.A, F2.H, F2.B]);
  assert(allTrueGenoPK(F2.A) == [F2pk.AA, F2pk.AB, F2pk.BA, F2pk.BB]);
}


// calcGenoprob for F2
// mixin calcGenoprobCode!F2;

unittest {
  writeln("    unit test calcGenoprob for F2:");
  alias std.path.buildPath buildPath;
  auto fn = dirName(__FILE__) ~ sep ~ buildPath("..","..","..","..","..","test","data","input","listeria.csv");
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
  auto rec_frac = dist_to_recfrac(dist_cM, GeneticMapFunc.Haldane);

  writeln("      - Run calcGenoprob for F2");
  auto genoprobs = calc_geno_prob(chr_4_genotypes, rec_frac, 0.002);

  writeln("      - Compare results to R/qtl");
  double[F2][size_t] genoprobs_from_rqtl;
  /* probs from R/qtl for individual 1 */
  genoprobs_from_rqtl[0] = [F2.A:0.99365258749878116, F2.H:0.005366774350785078, F2.B:0.00098063815043388934];
  genoprobs_from_rqtl[1] = [F2.A:0.01597337476839351, F2.H:0.984010456989931837, F2.B:0.00001616824167473149];
  genoprobs_from_rqtl[2] = [F2.A:0.98819056080900758, F2.H:0.010834274013193156, F2.B:0.00097516517779884938];
  genoprobs_from_rqtl[3] = [F2.A:0.00156449602454221, F2.H:0.998274388090304443, F2.B:0.00016111588515345984];
  
  foreach(i; 0..genoprobs[0].length) {
    foreach(j; [F2.A, F2.H, F2.B]) {
      assert(abs(genoprobs[0][cast()i][j] - genoprobs_from_rqtl[i][j]) < 1e-14);
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


// estmap
// mixin estmapCode!(F2, F2pk);

unittest {
  writeln("Unit test " ~ __FILE__);
  writeln("    unit test estmap for F2:");
  alias std.path.buildPath buildPath;
  auto fn = to!string(dirName(__FILE__) ~ sep ~ buildPath("..","..","..","..","..","test","data","input","listeria.csv"));
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
  auto rec_frac = dist_to_recfrac(dist_cM, GeneticMapFunc.Kosambi);

  auto rec_frac_rqtl = [0.18274786564786985044,
			0.15620001633845906341,
			0.28772927391795305452];

  foreach(i; 0..rec_frac.length) {
    assert(abs(rec_frac[i] - rec_frac_rqtl[i]) < 1e-14);
  }

  auto rec_frac_rev_rqtl = [0.18726831033621754719,
			    0.15929622543721855266,
			    0.29234793156548449788];

  auto rec_frac_rev = estmap(chr_4_genotypes, rec_frac, 0.002, 100, 1e-6, true);

  foreach(i; 0..rec_frac.length) {
    assert(abs(rec_frac_rev[i] - rec_frac_rev_rqtl[i]) < 1e-5);
  }
}


// forward and backward equations
// mixin forwardEquationsCode!(F2, F2);
// mixin backwardEquationsCode!(F2, F2);

// mixin forwardEquationsCode!(F2, F2pk);
// mixin backwardEquationsCode!(F2, F2pk);

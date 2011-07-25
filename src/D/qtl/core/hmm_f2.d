/**
 * hmm_f2
 **/

module qtl.core.hmm_f2;

import qtl.core.genotype;
import qtl.core.primitives;
import qtl.plugins.input.read_csv;
import std.math, std.stdio, std.path;

// marginal genotype probability
double initF2(F2 true_gen)
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
  writeln("    unit test initF2");
  auto gen = [F2.A, F2.H, F2.B];
  assert(abs(initF2(gen[0]) - log(0.25)) < 1e-12);
  assert(abs(initF2(gen[1]) - log(0.5)) < 1e-12);
  assert(abs(initF2(gen[2]) - log(0.25)) < 1e-12);
}


// emission probability (marker genotype "penetrance")
double emitF2(Genotype!F2 obs_gen, F2 true_gen, double error_prob) 
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
  writeln("    unit test emitF2");
  auto gen = [F2.A, F2.H, F2.B];
  auto ogen = [set_genotype!F2("-"), 
	       set_genotype!F2("A"), set_genotype!F2("H"), set_genotype!F2("B"),
	       set_genotype!F2("D"), set_genotype!F2("C")];
  double err = 0.001;

  /* missing */
  foreach(g; gen) {
    assert(abs(emitF2(ogen[0], g, err) ) < 1e-12);
  }

  /* A */
  assert(abs(emitF2(ogen[1], gen[0], err) - log(1.0-err)) < 1e-12);
  assert(abs(emitF2(ogen[1], gen[1], err) - log(err/2.0)) < 1e-12);
  assert(abs(emitF2(ogen[1], gen[2], err) - log(err/2.0)) < 1e-12);

  /* H */
  assert(abs(emitF2(ogen[2], gen[0], err) - log(err/2.0)) < 1e-12);
  assert(abs(emitF2(ogen[2], gen[1], err) - log(1.0-err)) < 1e-12);
  assert(abs(emitF2(ogen[2], gen[2], err) - log(err/2.0)) < 1e-12);

  /* B */
  assert(abs(emitF2(ogen[3], gen[0], err) - log(err/2.0)) < 1e-12);
  assert(abs(emitF2(ogen[3], gen[1], err) - log(err/2.0)) < 1e-12);
  assert(abs(emitF2(ogen[3], gen[2], err) - log(1.0-err)) < 1e-12);

  /* AorH (aka D) */
  assert(abs(emitF2(ogen[4], gen[0], err) - log(1.0-err/2.0)) < 1e-12);
  assert(abs(emitF2(ogen[4], gen[1], err) - log(1.0-err/2.0)) < 1e-12);
  assert(abs(emitF2(ogen[4], gen[2], err) - log(err)) < 1e-12);

  /* HorB (aka C) */
  assert(abs(emitF2(ogen[5], gen[0], err) - log(err)) < 1e-12);
  assert(abs(emitF2(ogen[5], gen[1], err) - log(1.0-err/2.0)) < 1e-12);
  assert(abs(emitF2(ogen[5], gen[2], err) - log(1.0-err/2.0)) < 1e-12);
}


// transition probabilities 
double stepF2(F2 true_gen_left, F2 true_gen_right, 
	      double rec_frac) 
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
  writeln("    unit test stepF2");
  auto gen = [F2.A, F2.H, F2.B];
  double rf = 0.01;

  assert( abs( stepF2(gen[0], gen[0], rf) - log((1-rf)^^2) ) < 1e-12);
  assert( abs( stepF2(gen[0], gen[1], rf) - log(2.0*rf*(1-rf)) ) < 1e-12);
  assert( abs( stepF2(gen[0], gen[2], rf) - log(rf^^2) ) < 1e-12);

  assert( abs( stepF2(gen[1], gen[0], rf) - log(rf*(1-rf)) ) < 1e-12);
  assert( abs( stepF2(gen[1], gen[1], rf) - log((1-rf)^^2 + rf^^2) ) < 1e-12);
  assert( abs( stepF2(gen[1], gen[2], rf) - log(rf*(1-rf)) ) < 1e-12);

  assert( abs( stepF2(gen[2], gen[2], rf) - log((1-rf)^^2) ) < 1e-12);
  assert( abs( stepF2(gen[2], gen[1], rf) - log(2.0*rf*(1-rf)) ) < 1e-12);
  assert( abs( stepF2(gen[2], gen[0], rf) - log(rf^^2) ) < 1e-12);
}

  
double initF2pk(F2pk true_gen)
{
  if(true_gen != F2pk.AA && true_gen != F2pk.AB &&
     true_gen != F2pk.BA && true_gen != F2pk.BB)
    throw new Exception("true_gen not a valid genotype.");

  return(-2.0*LN2);  /* ln(0.25) */
}

unittest {
  writeln("    unit test initF2pk");
  auto gen = [F2pk.AA, F2pk.AB, F2pk.BA, F2pk.BB];
  foreach(g; gen) {
    assert(abs(initF2pk(g) - log(0.25)) < 1e-12);
  }
}




double emitF2pk(Genotype!F2 obs_gen, F2pk true_gen, double error_prob)
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
  writeln("    unit test emitF2pk");

  auto gen = [F2pk.AA, F2pk.AB, F2pk.BA, F2pk.BB];
  auto ogen = [set_genotype!F2("-"), 
	       set_genotype!F2("A"), set_genotype!F2("H"), set_genotype!F2("B"),
	       set_genotype!F2("D"), set_genotype!F2("C")];
  double err = 0.001;

  /* missing */
  foreach(g; gen) {
    assert(abs(emitF2pk(ogen[0], g, err) ) < 1e-12);
  }

  /* A */
  assert(abs(emitF2pk(ogen[1], gen[0], err) - log(1.0-err)) < 1e-12);
  assert(abs(emitF2pk(ogen[1], gen[1], err) - log(err/2.0)) < 1e-12);
  assert(abs(emitF2pk(ogen[1], gen[2], err) - log(err/2.0)) < 1e-12);
  assert(abs(emitF2pk(ogen[1], gen[3], err) - log(err/2.0)) < 1e-12);

  /* H */
  assert(abs(emitF2pk(ogen[2], gen[0], err) - log(err/2.0)) < 1e-12);
  assert(abs(emitF2pk(ogen[2], gen[1], err) - log(1.0-err)) < 1e-12);
  assert(abs(emitF2pk(ogen[2], gen[2], err) - log(1.0-err)) < 1e-12);
  assert(abs(emitF2pk(ogen[2], gen[3], err) - log(err/2.0)) < 1e-12);

  /* B */
  assert(abs(emitF2pk(ogen[3], gen[0], err) - log(err/2.0)) < 1e-12);
  assert(abs(emitF2pk(ogen[3], gen[1], err) - log(err/2.0)) < 1e-12);
  assert(abs(emitF2pk(ogen[3], gen[2], err) - log(err/2.0)) < 1e-12);
  assert(abs(emitF2pk(ogen[3], gen[3], err) - log(1.0-err)) < 1e-12);

  /* AorH (aka D) */
  assert(abs(emitF2pk(ogen[4], gen[0], err) - log(1.0-err/2.0)) < 1e-12);
  assert(abs(emitF2pk(ogen[4], gen[1], err) - log(1.0-err/2.0)) < 1e-12);
  assert(abs(emitF2pk(ogen[4], gen[2], err) - log(1.0-err/2.0)) < 1e-12);
  assert(abs(emitF2pk(ogen[4], gen[3], err) - log(err)) < 1e-12);

  /* HorB (aka C) */
  assert(abs(emitF2pk(ogen[5], gen[0], err) - log(err)) < 1e-12);
  assert(abs(emitF2pk(ogen[5], gen[1], err) - log(1.0-err/2.0)) < 1e-12);
  assert(abs(emitF2pk(ogen[5], gen[2], err) - log(1.0-err/2.0)) < 1e-12);
  assert(abs(emitF2pk(ogen[5], gen[3], err) - log(1.0-err/2.0)) < 1e-12);
}


double stepF2pk(F2pk true_gen_left, F2pk true_gen_right, double rec_frac)
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
  writeln("    unit test stepF2pk");
  auto gen = [F2pk.AA, F2pk.AB, F2pk.BA, F2pk.BB];
  double rf = 0.01;

  assert( abs( stepF2pk(gen[0], gen[0], rf) - log((1-rf)^^2) ) < 1e-12);
  assert( abs( stepF2pk(gen[0], gen[1], rf) - log(rf*(1-rf)) ) < 1e-12);
  assert( abs( stepF2pk(gen[0], gen[2], rf) - log(rf*(1-rf)) ) < 1e-12);
  assert( abs( stepF2pk(gen[0], gen[3], rf) - log(rf^^2)     ) < 1e-12);

  assert( abs( stepF2pk(gen[1], gen[0], rf) - log(rf*(1-rf)) ) < 1e-12);
  assert( abs( stepF2pk(gen[1], gen[1], rf) - log((1-rf)^^2) ) < 1e-12);
  assert( abs( stepF2pk(gen[1], gen[2], rf) - log(rf^^2)     ) < 1e-12);
  assert( abs( stepF2pk(gen[1], gen[3], rf) - log(rf*(1-rf)) ) < 1e-12);

  assert( abs( stepF2pk(gen[2], gen[3], rf) - log(rf*(1-rf)) ) < 1e-12);
  assert( abs( stepF2pk(gen[2], gen[2], rf) - log((1-rf)^^2) ) < 1e-12);
  assert( abs( stepF2pk(gen[2], gen[1], rf) - log(rf^^2)     ) < 1e-12);
  assert( abs( stepF2pk(gen[2], gen[0], rf) - log(rf*(1-rf)) ) < 1e-12);

  assert( abs( stepF2pk(gen[3], gen[3], rf) - log((1-rf)^^2) ) < 1e-12);
  assert( abs( stepF2pk(gen[3], gen[2], rf) - log(rf*(1-rf)) ) < 1e-12);
  assert( abs( stepF2pk(gen[3], gen[1], rf) - log(rf*(1-rf)) ) < 1e-12);
  assert( abs( stepF2pk(gen[3], gen[0], rf) - log(rf^^2)     ) < 1e-12);
}



double nrecF2pk(F2pk true_gen_left, F2pk true_gen_right)
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
  writeln("    unit test nrecF2pk");
  auto gen = [F2pk.AA, F2pk.AB, F2pk.BA, F2pk.BB];

  assert(nrecF2pk(gen[0], gen[0])==0.0);
  assert(nrecF2pk(gen[0], gen[1])==0.5);
  assert(nrecF2pk(gen[0], gen[2])==0.5);
  assert(nrecF2pk(gen[0], gen[3])==1.0);

  assert(nrecF2pk(gen[1], gen[0])==0.5);
  assert(nrecF2pk(gen[1], gen[1])==0.0);
  assert(nrecF2pk(gen[1], gen[2])==1.0);
  assert(nrecF2pk(gen[1], gen[3])==0.5);

  assert(nrecF2pk(gen[2], gen[3])==0.5);
  assert(nrecF2pk(gen[2], gen[2])==0.0);
  assert(nrecF2pk(gen[2], gen[1])==1.0);
  assert(nrecF2pk(gen[2], gen[0])==0.5);

  assert(nrecF2pk(gen[3], gen[3])==0.0);
  assert(nrecF2pk(gen[3], gen[2])==0.5);
  assert(nrecF2pk(gen[3], gen[1])==0.5);
  assert(nrecF2pk(gen[3], gen[0])==1.0);
}

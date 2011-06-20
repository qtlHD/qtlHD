/**
 * hmm_f2
 **/

module qtl.core.hmm_f2;

import qtl.core.genotype;
import qtl.core.primitives;
import qtl.plugins.input.read_csv;
import std.math, std.stdio, std.path;

// marginal genotype probability
double initF2(Genotype!F2 true_gen)
in {
  assert(true_gen.value == F2.A || 
	 true_gen.value == F2.H || 
	 true_gen.value == F2.B);
}
body {
  switch(true_gen.value) {
    case F2.H:
      return(-LN2);
    case F2.A: case F2.B: 
      return(-2.0*LN2);
  }
  return(0.0); /* shouldn't get here */
}

unittest {
  writeln("Unit test " ~ __FILE__);
  writeln("    unit test initFC");
  auto gen = [set_genotype!F2("A"), set_genotype!F2("H"), set_genotype!F2("B")];
  assert(abs(initF2(gen[0]) - log(0.25)) < 1e-12);
  assert(abs(initF2(gen[1]) - log(0.5)) < 1e-12);
  assert(abs(initF2(gen[2]) - log(0.25)) < 1e-12);
}


// emission probability (marker genotype "penetrance")
double emitF2(Genotype!F2 obs_gen, Genotype!F2 true_gen, double error_prob) 
in {
  assert(true_gen.value == F2.A || 
	 true_gen.value == F2.H || 
	 true_gen.value == F2.B);
  assert(error_prob >= 0 && error_prob <= 1.0);
}
body {
  switch(obs_gen.value) {
  case F2.NA: return(0.0);
  case F2.A: case F2.H: case F2.B:
    if(obs_gen.value == true_gen.value) {
      return(log(1.0-error_prob));
    } else {
      return(log(error_prob)-LN2);
    }
  case F2.HorB:
    if(true_gen.value == F2.A) {
      return(log(error_prob));
    } else {
      return(log(1.0-error_prob/2.0));
    }
  case F2.HorA:
    if(true_gen.value == F2.B) {
      return(log(error_prob));
    } else {
      return(log(1.0-error_prob/2.0));
    }
  }
  return(0.0); /* shouldn't get here */
}

unittest {
  writeln("    unit test emitFC");
  auto gen = [set_genotype!F2("A"), set_genotype!F2("H"), set_genotype!F2("B")];
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
double stepF2(Genotype!F2 true_gen_left, Genotype!F2 true_gen_right, 
	      double rec_frac) 
in {
  assert(true_gen_left.value == F2.A || 
	 true_gen_left.value == F2.H || 
	 true_gen_left.value == F2.B);
  assert(true_gen_right.value == F2.A || 
	 true_gen_right.value == F2.H || 
	 true_gen_right.value == F2.B);
  assert(rec_frac >= 0 && rec_frac <= 0.5);
}
body {
  switch(true_gen_left.value) {
    case F2.A:
      switch(true_gen_right.value) {
      case F2.A: return(2.0*log(1.0-rec_frac));
      case F2.H: return(LN2 + log(1.0-rec_frac) + log(rec_frac));
      case F2.B: return(2.0*log(rec_frac));
      }
    case F2.H:
      switch(true_gen_right.value) {
      case F2.A: case F2.B: return(log(rec_frac) + log(1.0-rec_frac));
      case F2.H: return(log((1.0-rec_frac)^^2 + rec_frac^^2));
      }
    case F2.B:
      switch(true_gen_right.value) {
      case F2.A: return(2.0*log(rec_frac));
      case F2.H: return(LN2 + log(1.0-rec_frac) + log(rec_frac));
      case F2.B: return(2.0*log(1.0-rec_frac));
      }
  }
  return(log(-1.0)); /* shouldn't get here */
}

unittest {
  writeln("    unit test stepFC");
  auto gen = [set_genotype!F2("A"), set_genotype!F2("H"), set_genotype!F2("B")];
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

  

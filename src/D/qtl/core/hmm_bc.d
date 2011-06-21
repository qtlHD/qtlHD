/**
 * hmm_bc
 **/

module qtl.core.hmm_bc;

import qtl.core.genotype;
import qtl.core.primitives;
import qtl.plugins.input.read_csv;
import std.math, std.stdio, std.path;

// marginal genotype probability
double initBC(BC true_gen)
in {
  assert(true_gen == BC.A || true_gen == BC.H,
	 "true_gen not among the possible true genotypes");
}
body {
  return(-LN2);
}

unittest {
  writeln("Unit test " ~ __FILE__);
  writeln("    unit test initBC");
  auto gen = [BC.A, BC.H];
  assert(abs(initBC(gen[0]) - log(0.5)) < 1e-12);
  assert(abs(initBC(gen[1]) - log(0.5)) < 1e-12);
}

// emission probability (marker genotype "penetrance")
double emitBC(Genotype!BC obs_gen, BC true_gen, double error_prob) 
in {
  assert(true_gen == BC.A || true_gen == BC.H,
	 "true_gen not among the possible true genotypes");
  assert(error_prob >= 0 && error_prob <= 1.0,
	 "error_prob must be >= 0 and <= 1");
}
body {
  switch(obs_gen.value) {
    case BC.NA: return(0.0);
    case BC.A: case BC.H:
      if(obs_gen.value == true_gen) {
        return(log(1.0-error_prob));
      } else {
        return(log(error_prob));
      }
  }
  return(0.0); /* shouldn't get here */
}

unittest {
  writeln("    unit test emitBC");
  auto gen = [BC.A, BC.H];
  auto ogen = [set_genotype!BC("-"), set_genotype!BC("A"), set_genotype!BC("H")];
  double err = 0.001;

  foreach(g; gen) {
    assert(abs( emitBC(ogen[0], g, err) ) < 1e-12);
  }
  assert(abs( emitBC(ogen[1], gen[0], err) - log(1-err) ) < 1e-12);
  assert(abs( emitBC(ogen[1], gen[1], err) - log(err) ) < 1e-12);
  assert(abs( emitBC(ogen[2], gen[0], err) - log(err) ) < 1e-12);
  assert(abs( emitBC(ogen[2], gen[1], err) - log(1-err) ) < 1e-12);
}

// transition probabilities 
double stepBC(BC true_gen_left, BC true_gen_right, 
	      double rec_frac)
in {
  assert(true_gen_left == BC.A || true_gen_left == BC.H,
	 "true_gen_left not among the possible true genotypes");
  assert(true_gen_right == BC.A || true_gen_right == BC.H,
	 "true_gen_right not among the possible true genotypes");
  assert(rec_frac >= 0.0 && rec_frac <= 0.5,
	 "rec_frac must be >= 0 and <= 1");
}
body {
  if(true_gen_left == true_gen_right) {
    return(log(1.0-rec_frac));
  }
  else {
    return(log(rec_frac));
  }
}


unittest {
  writeln("    unit test stepBC");
  auto gen = [BC.A, BC.H];
  double rf = 0.01;
  assert(abs( stepBC(gen[0], gen[0], rf) - log(1-rf) ) < 1e-12);
  assert(abs( stepBC(gen[0], gen[1], rf) - log(rf) ) < 1e-12);
  assert(abs( stepBC(gen[1], gen[0], rf) - log(rf) ) < 1e-12);
  assert(abs( stepBC(gen[1], gen[1], rf) - log(1-rf) ) < 1e-12);
}

  

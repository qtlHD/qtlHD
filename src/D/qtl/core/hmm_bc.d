/**
 * hmm_bc
 **/

module qtl.core.hmm_bc;

import qtl.core.genotype;
import qtl.core.primitives;
import qtl.plugins.input.read_csv;
import std.math, std.stdio, std.path;
import std.exception;


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

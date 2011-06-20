
module qtl.core.hmm_bc;
import qtl.core.genotype;

enum BCtrue { A=BC.A, H=BC.H };


// marginal genotype probability
double initBC(Genotype!BCtrue true_gen)
{
  return(-LN2);
}

// emission probability (marker genotype "penetrance")
double emitBC(Genotype!BC obs_gen, Genotype!BCtrue true_gen, double error_prob) 
{
  switch(obs_gen.value) {
    case BC.NA: return(0.0);
    case BC.A: case BC.H:
      if(obs_gen.value == true_gen.value) {
        return(log(1.0-error_prob));
      } else {
        return(log(error_prob));
      }
  }
  return(0.0); /* shouldn't get here */
}

// transition probabilities 
double stepBC(Genotype!BCtrue true_gen_left, Genotype!BCtrue true_gen_right, 
	      double rec_frac)
{
  if(true_gen_left.value == true_gen_right.value) {
    return(log(1.0-rec_frac));
  }
  else {
    return(log(rec_frac));
  }
}


unittest {
  writeln("Unit test " ~ __FILE__);
  alias std.path.join join;
}

  

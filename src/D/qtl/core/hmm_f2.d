
module qtl.core.hmm_f2;
import qtl.core.genotype;

enum F2true { A=F2.A, H=F2.H, B=F2.B };


// marginal genotype probability
double initF2(Genotype!F2true true_gen)
{
  switch(true_gen.value) {
    case F2true.H:
      return(-LN2);
    case F2true.A: case F2true.B: 
      return(-2.0*LN2);
  }
  return(0.0); /* shouldn't get here */
}

// emission probability (marker genotype "penetrance")
double emitF2(Genotype!F2 obs_gen, Genotype!F2true true_gen, double error_prob) 
{
  switch(obs_gen.value) {
  case F2.NA: return(0.0);
  case F2.A: case F2.H: case F2.B:
    if(obs_gen.value == true_gen.value) {
      return(log(1.0-error_prob));
    } else {
      return(log(error_prob)-LN2);
    }
  case F2.HorB:
    if(true_gen.value != F2.A) {
      return(log(1.0-error_prob/2.0));
    } else {
      return(log(error_prob)-LN2);
    }
  case F2.HorA:
    if(true_gen.value != F2.B) {
      return(log(1.0-error_prob/2.0));
    } else {
      return(log(error_prob)-LN2);
    }
  }
  return(0.0); /* shouldn't get here */
}

// transition probabilities 
double stepF2(Genotype!F2true true_gen_left, Genotype!F2true true_gen_right, 
	      double rec_frac) {
  switch(true_gen_left.value) {
    case F2true.A:
      switch(true_gen_right.value) {
      case F2true.A: return(2.0*log(1.0-rec_frac));
      case F2true.H: return(LN2 + log(1.0-rec_frac) + log(rec_frac));
      case F2true.B: return(2.0*log(rec_frac));
      }
    case F2true.H:
      switch(true_gen_right.value) {
      case F2true.A: case F2true.B: return(log(rec_frac) + log(1.0-rec_frac));
      case F2true.H: return(log((1.0-rec_frac)^^2 + rec_frac^^2));
      }
    case F2true.B:
      switch(true_gen_right.value) {
      case F2true.A: return(2.0*log(rec_frac));
      case F2true.H: return(LN2 + log(1.0-rec_frac) + log(rec_frac));
      case F2true.B: return(2.0*log(1.0-rec_frac));
      }
  }
  return(log(-1.0)); /* shouldn't get here */
}


unittest {
  writeln("Unit test " ~ __FILE__);
  alias std.path.join join;
  auto fn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","test","data","input","listeria.csv");
  writeln("  - read CSV " ~ fn);
  auto data = new ReadSimpleCSV!F2(fn);
}

  

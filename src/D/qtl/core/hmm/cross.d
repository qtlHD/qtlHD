/**
 * HMM backcross functions
 **/

module qtl.core.hmm.cross;

import std.math, std.stdio, std.path;
import std.exception;
import std.conv;

import qtl.core.primitives;
import qtl.core.genotype;
import qtl.core.hmm.calcgenoprob;
import qtl.core.hmm.estmap;
import qtl.core.map.genetic_map_functions;

// class to contain HMM-related functions
class Cross {
  string cross_type;
  
  TrueGenotype[] all_true_geno, all_true_geno_phaseknown;
  
  abstract double init(in TrueGenotype truegen);
  abstract double init_phaseknown(in TrueGenotype truegen);

  abstract double step(in TrueGenotype truegen_left, in TrueGenotype truegen_right, in double rec_frac);
  abstract double step_phaseknown(in TrueGenotype truegen_left, in TrueGenotype truegen_right, in double rec_frac);

  abstract double emit(in GenotypeCombinator obsgen, in TrueGenotype truegen, in double error_prob);

  abstract double nrec(in TrueGenotype truegen_left, in TrueGenotype truegen_right);

}

// the switch
Cross form_cross(string which_cross) {
  Cross cross;

  switch(which_cross) {
  case "BC": cross = new BC; break;
  case "F2": cross = new F2; break;
  default: throw new Exception("cross type " ~ which_cross ~ " not supported.");
  }
  return cross;
}


// BC (backcross)
class BC : Cross {
  this() {
    cross_type = "BC";

    // vector of true genotypes
    auto g1 = new TrueGenotype(0,0);
    auto g2 = new TrueGenotype(1,0);
    all_true_geno = [g1, g2];
    all_true_geno_phaseknown = all_true_geno;
  }

  // ln Pr(true genotype)
  double init(in TrueGenotype truegen)
  {
    if(truegen != all_true_geno[0] &&
       truegen != all_true_geno[1])
      throw new Exception("truegen must be 0,0 or 1,0");
    
    return(-LN2);
  }

  double init_phaseknown(in TrueGenotype truegen) { return init(truegen); }

  // ln Pr(genotype at right marker | genotype at left marker)
  double step(in TrueGenotype truegen_left, in TrueGenotype truegen_right, in double rec_frac)
  {
    alias all_true_geno atg;

    if(truegen_left == atg[0]) {
      if(truegen_right == atg[0]) // A -> A
        return(log(1.0-rec_frac));
      if(truegen_right == atg[1]) // A -> H
        return(log(rec_frac));
    }

    if(truegen_left == atg[1]) {
      if(truegen_right == atg[1]) // H -> H
        return(log(1.0-rec_frac));
      if(truegen_right == atg[0]) // H -> A
        return(log(rec_frac));
    }

    throw new Exception("inputs not among the possible true genotypes");
  }

  double step_phaseknown(in TrueGenotype truegen_left, in TrueGenotype truegen_right,
                         in double rec_frac) {
    return step(truegen_left, truegen_right, rec_frac);
  }

  // ln Pr(observed genotype | true genotype)
  double emit(in GenotypeCombinator obsgen, in TrueGenotype truegen, in double error_prob)
  {
    if(obsgen.list.length==0) // missing value
      return(0.0); // log(1.0)

    if(obsgen.match(truegen)) // compatible with truegen?
      return(log(1.0-error_prob));
    else
      return(log(error_prob));
  }

  // No. recombination events
  double nrec(in TrueGenotype truegen_left, in TrueGenotype truegen_right)
  {
    alias all_true_geno atg;

    if(truegen_left == atg[0]) {
      if(truegen_right == atg[0]) // A -> A
        return(0.0);
      if(truegen_right == atg[1]) // A -> H
        return(1.0);
    }

    if(truegen_left == atg[1]) {
      if(truegen_right == atg[1]) // H -> H
        return(0.0);
      if(truegen_right == atg[0]) // H -> A
        return(1.0);
    }

    throw new Exception("inputs not among the possible true genotypes");
  }
}

unittest {
  writeln("Unit test " ~ __FILE__ ~ " BC");

  // unit test all_true_geno for BC
  auto bc = form_cross("BC");
  auto atg = bc.all_true_geno;
  assert(atg.length == 2);
  assert(atg[0] == new TrueGenotype(0,0));
  assert(atg[1] == new TrueGenotype(1,0));

  // unit test init for BC
  foreach(gv; atg)
    assert(to!string(bc.init(gv)) == to!string(log(0.5)));

  // unit test step for BC
  double rf = 0.01;

  assert( to!string( bc.step(atg[0], atg[0], rf) ) ==
          to!string( log((1-rf)) ));
  assert( to!string( bc.step(atg[0], atg[1], rf) ) ==
          to!string( log(rf) ));

  assert( to!string( bc.step(atg[1], atg[0], rf) ) ==
          to!string( log(rf) ));
  assert( to!string( bc.step(atg[1], atg[1], rf) ) ==
          to!string( log((1-rf)) ));

  // unit test emit for BC
  auto NA = new GenotypeCombinator("-");
  auto A = new GenotypeCombinator("A");
  A ~= new TrueGenotype(0,0);
  auto H = new GenotypeCombinator("H");
  H ~= new TrueGenotype(1,0);

  double error_prob = 0.01;

  assert(bc.emit(NA, atg[0], error_prob) == 0);
  assert(bc.emit(NA, atg[1], error_prob) == 0);
  assert(to!string(bc.emit(A, atg[0], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(bc.emit(A, atg[1], error_prob)) == to!string(log(error_prob)));
  assert(to!string(bc.emit(H, atg[1], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(bc.emit(H, atg[0], error_prob)) == to!string(log(error_prob)));

  // unit test nrec for BC
  assert(bc.nrec(atg[0], atg[0]) == 0.0);
  assert(bc.nrec(atg[0], atg[1]) == 1.0);

  assert(bc.nrec(atg[1], atg[0]) == 1.0);
  assert(bc.nrec(atg[1], atg[1]) == 0.0);
}


// F2 (intercross)
class F2 : Cross {

  this() {
    cross_type = "F2";

    auto g1 = new TrueGenotype(0,0);
    auto g2 = new TrueGenotype(1,0);
    auto g3 = new TrueGenotype(0,1);
    auto g4 = new TrueGenotype(1,1);
    all_true_geno = [g1, g2, g4];
    all_true_geno_phaseknown = [g1, g2, g3, g4];
  }

  // ln Pr(true genotype)
  double init(in TrueGenotype truegen)
  {
    alias all_true_geno atg;

    if(truegen==atg[0] || truegen==atg[2]) // AA or BB
      return(-2.0*LN2);

    if(truegen==atg[1]) // AB
      return(-LN2);

    throw new Exception("truegen not among possible true genotypes");
  }

  // ln Pr(true genotype) for est_map
  double init_phaseknown(in TrueGenotype truegen)
  {
    alias all_true_geno_phaseknown atg;

    if(truegen==atg[0] || truegen==atg[1] ||
       truegen==atg[2] || truegen==atg[3])
      return(-2.0*LN2);

    throw new Exception("truegen not among possible true genotypes");
  }

  // ln Pr(genotype at right marker | genotype at left marker)
  double step(in TrueGenotype truegen_left, in TrueGenotype truegen_right, in double rec_frac)
  {
    alias all_true_geno atg;

    if(truegen_left==atg[0]) {
      if(truegen_right==atg[0]) // A -> A
        return(2.0*log(1.0-rec_frac));
      if(truegen_right==atg[1]) // A -> H
        return(LN2 + log(1.0-rec_frac) + log(rec_frac));
      if(truegen_right==atg[2]) // A -> B
        return(2.0*log(rec_frac));
    }
    else if(truegen_left == atg[1]) {
      if(truegen_right==atg[0] || // H -> A
         truegen_right==atg[2]) // H -> B
        return(log(rec_frac) + log(1.0-rec_frac));
      if(truegen_right==atg[1]) // H -> H
        return(log((1.0-rec_frac)^^2 + rec_frac^^2));
    }
    else if(truegen_left == atg[2]) {
      if(truegen_right==atg[0]) // B -> A
        return(2.0*log(rec_frac));
      if(truegen_right==atg[1]) // B -> H
        return(LN2 + log(1.0-rec_frac) + log(rec_frac));
      if(truegen_right==atg[2]) // B -> B
        return(2.0*log(1.0-rec_frac));
    }
    
    throw new Exception("inputs not among the possible true genotypes");
  }
  
  // like step, but with phase-known genotypes (for use with est_map)
  double step_phaseknown(in TrueGenotype truegen_left, in TrueGenotype truegen_right, in double rec_frac)
  {
    alias all_true_geno_phaseknown atgpk;

    if(truegen_left==atgpk[0]) {
      if(truegen_right==atgpk[0]) // AA -> AA
        return(2.0*log(1.0-rec_frac));
      if(truegen_right==atgpk[1] || truegen_right==atgpk[2]) // AA -> AB; AA -> BA
        return(log(1.0-rec_frac) + log(rec_frac));
      if(truegen_right==atgpk[3]) // AA -> BB
        return(2.0*log(rec_frac));
    }
    else if(truegen_left == atgpk[1]) {
      if(truegen_right==atgpk[0] || // AB -> AA
         truegen_right==atgpk[3]) // AB -> BB
        return(log(rec_frac) + log(1.0-rec_frac));
      if(truegen_right==atgpk[1]) // AB -> AB
        return(log((1.0-rec_frac)^^2));
      if(truegen_right==atgpk[2]) // AB -> BA
        return(log(rec_frac^^2));
    }
    else if(truegen_left == atgpk[2]) {
      if(truegen_right==atgpk[0] || // BA -> AA
         truegen_right==atgpk[3]) // BA -> BB
        return(log(rec_frac) + log(1.0-rec_frac));
      if(truegen_right==atgpk[2]) // BA -> BA
        return(log((1.0-rec_frac)^^2));
      if(truegen_right==atgpk[1]) // BA -> AB
        return(log(rec_frac^^2));
    }
    else if(truegen_left == atgpk[3]) {
      if(truegen_right==atgpk[0]) // BB -> AA
        return(2.0*log(rec_frac));
      if(truegen_right==atgpk[1] || // BB -> AB
         truegen_right==atgpk[2])   // BB -> BA
        return(log(1.0-rec_frac) + log(rec_frac));
      if(truegen_right==atgpk[3]) // BB -> BB
        return(2.0*log(1.0-rec_frac));
    }
    
    throw new Exception("inputs not among the possible true genotypes");
  }

  // ln Pr(observed genotype | true genotype)
  double emit(in GenotypeCombinator obsgen, in TrueGenotype truegen, in double error_prob)
  {
    auto n_obsgen = obsgen.list.length;

    if(n_obsgen==0) // missing value
      return(0.0); // log(1.0)

    if(obsgen.match(truegen)) { // compatible with truegen
      if(n_obsgen>2) // AorH and HorB cases
        return(log(1.0-error_prob/2.0));
      else // A, H, B cases
        return(log(1.0-error_prob));
    }
    else {
      if(n_obsgen>2) // AorH and HorB cases
        return(log(error_prob));
      else // A, H, B cases
        return(log(error_prob/2.0));
    }
  }

  // proportion of recombination events
  double nrec(in TrueGenotype truegen_left, in TrueGenotype truegen_right)
  {
    alias all_true_geno_phaseknown atgpk;

    if(truegen_left==atgpk[0]) {
      if(truegen_right==atgpk[0]) // AA -> AA
        return(0.0);
      if(truegen_right==atgpk[1] || truegen_right==atgpk[2]) // AA -> AB; AA -> BA
        return(0.5);
      if(truegen_right==atgpk[3]) // AA -> BB
        return(1.0);
    }
    else if(truegen_left == atgpk[1]) {
      if(truegen_right==atgpk[0] || // AB -> AA
         truegen_right==atgpk[3]) // AB -> BB
        return(0.5);
      if(truegen_right==atgpk[1]) // AB -> AB
        return(0.0);
      if(truegen_right==atgpk[2]) // AB -> BA
        return(1.0);
    }
    else if(truegen_left == atgpk[2]) {
      if(truegen_right==atgpk[0] || // BA -> AA
         truegen_right==atgpk[3]) // BA -> BB
        return(0.5);
      if(truegen_right==atgpk[2]) // BA -> BA
        return(0.0);
      if(truegen_right==atgpk[1]) // BA -> AB
        return(1.0);
    }
    else if(truegen_left == atgpk[3]) {
      if(truegen_right==atgpk[0]) // BB -> AA
        return(1.0);
      if(truegen_right==atgpk[1] || // BB -> AB
         truegen_right==atgpk[2])   // BB -> BA
        return(0.5);
      if(truegen_right==atgpk[3]) // BB -> BB
        return(0.0);
    }
    
    throw new Exception("inputs not among the possible true genotypes");
  }

}

unittest {
  writeln("Unit test " ~ __FILE__ ~ " F2");

  // unit test allTrueGeno for F2
  auto f2 = form_cross("F2");
  auto g = f2.all_true_geno;
  auto gPK = f2.all_true_geno_phaseknown;

  assert(g.length == 3);
  assert(g[0] == new TrueGenotype(0,0));
  assert(g[1] == new TrueGenotype(1,0));
  assert(g[2] == new TrueGenotype(1,1));

  // unit test allTrueGeno for F2 (phase-known)
  assert(gPK.length == 4);
  assert(gPK[0] == new TrueGenotype(0,0));
  assert(gPK[1] == new TrueGenotype(1,0));
  assert(gPK[2] == new TrueGenotype(0,1));
  assert(gPK[3] == new TrueGenotype(1,1));

  // unit test init for F2
  double val[] = [log(0.25), log(0.5), log(0.25)];

  foreach(i, gv; g)
    assert(to!string(f2.init(gv)) == to!string(val[i]));

  // unit test init for F2 (phase-known)
  foreach(gv; gPK)
    assert(to!string(f2.init_phaseknown(gv)) == to!string(val[0]));

  //  unit test step for F2
  double rf = 0.01;

  assert( to!string( f2.step(g[0], g[0], rf) ) ==
          to!string( log((1-rf)^^2) ));
  assert( to!string( f2.step(g[0], g[1], rf) ) ==
          to!string( log(2.0*rf*(1-rf)) ));
  assert( to!string( f2.step(g[0], g[2], rf) ) ==
          to!string( log(rf^^2) ));

  assert( to!string( f2.step(g[1], g[0], rf) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( f2.step(g[1], g[1], rf) ) ==
          to!string( log(rf^^2 + (1-rf)^^2) ));
  assert( to!string( f2.step(g[1], g[2], rf) ) ==
          to!string( log(rf*(1-rf)) ));

  assert( to!string( f2.step(g[2], g[0], rf) ) ==
          to!string( log(rf^^2) ));
  assert( to!string( f2.step(g[2], g[1], rf) ) ==
          to!string( log(2.0*rf*(1-rf)) ));
  assert( to!string( f2.step(g[2], g[2], rf) ) ==
          to!string( log((1-rf)^^2) ));


  // unit test step for F2 (phase-known)
  assert( to!string( f2.step_phaseknown(gPK[0], gPK[0], rf) ) ==
          to!string( log((1-rf)^^2) ));
  assert( to!string( f2.step_phaseknown(gPK[0], gPK[1], rf) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( f2.step_phaseknown(gPK[0], gPK[2], rf) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( f2.step_phaseknown(gPK[0], gPK[3], rf) ) ==
          to!string( log(rf^^2) ));

  assert( to!string( f2.step_phaseknown(gPK[1], gPK[0], rf) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( f2.step_phaseknown(gPK[1], gPK[1], rf) ) ==
          to!string( log((1-rf)^^2) ));
  assert( to!string( f2.step_phaseknown(gPK[1], gPK[2], rf) ) ==
          to!string( log(rf^^2) ));
  assert( to!string( f2.step_phaseknown(gPK[1], gPK[3], rf) ) ==
          to!string( log(rf*(1.0-rf)) ));

  assert( to!string( f2.step_phaseknown(gPK[2], gPK[0], rf) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( f2.step_phaseknown(gPK[2], gPK[2], rf) ) ==
          to!string( log((1-rf)^^2) ));
  assert( to!string( f2.step_phaseknown(gPK[2], gPK[1], rf) ) ==
          to!string( log(rf^^2) ));
  assert( to!string( f2.step_phaseknown(gPK[2], gPK[3], rf) ) ==
          to!string( log(rf*(1.0-rf)) ));

  assert( to!string( f2.step_phaseknown(gPK[3], gPK[3], rf) ) ==
          to!string( log((1-rf)^^2) ));
  assert( to!string( f2.step_phaseknown(gPK[3], gPK[1], rf) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( f2.step_phaseknown(gPK[3], gPK[2], rf) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( f2.step_phaseknown(gPK[3], gPK[0], rf) ) ==
          to!string( log(rf^^2) ));

  // unit test emit for F2
  auto NA = new GenotypeCombinator("-");
  auto A = new GenotypeCombinator("A");
  A ~= new TrueGenotype(0,0);
  auto H = new GenotypeCombinator("H");
  H ~= new TrueGenotype(1,0);
  H ~= new TrueGenotype(0,1);
  auto B = new GenotypeCombinator("B");
  B ~= new TrueGenotype(1,1);
  auto AorH = new GenotypeCombinator("AorH");
  AorH ~= new TrueGenotype(0,0);
  AorH ~= new TrueGenotype(1,0);
  AorH ~= new TrueGenotype(0,1);
  auto HorB = new GenotypeCombinator("HorB");
  HorB ~= new TrueGenotype(1,0);
  HorB ~= new TrueGenotype(0,1);
  HorB ~= new TrueGenotype(1,1);

  double error_prob = 0.01;

  assert(f2.emit(NA, g[0], error_prob) == 0);
  assert(f2.emit(NA, g[1], error_prob) == 0);
  assert(f2.emit(NA, g[2], error_prob) == 0);

  assert(to!string(f2.emit(A, g[0], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2.emit(A, g[1], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2.emit(A, g[2], error_prob)) == to!string(log(error_prob/2.0)));

  assert(to!string(f2.emit(H, g[1], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2.emit(H, g[0], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2.emit(H, g[2], error_prob)) == to!string(log(error_prob/2.0)));

  assert(to!string(f2.emit(B, g[2], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2.emit(B, g[0], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2.emit(B, g[1], error_prob)) == to!string(log(error_prob/2.0)));

  assert(to!string(f2.emit(AorH, g[0], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2.emit(AorH, g[1], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2.emit(AorH, g[2], error_prob)) == to!string(log(error_prob)));

  assert(to!string(f2.emit(HorB, g[0], error_prob)) == to!string(log(error_prob)));
  assert(to!string(f2.emit(HorB, g[1], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2.emit(HorB, g[2], error_prob)) == to!string(log(1.0-error_prob/2.0)));

  assert(f2.emit(NA, gPK[0], error_prob) == 0);
  assert(f2.emit(NA, gPK[1], error_prob) == 0);
  assert(f2.emit(NA, gPK[2], error_prob) == 0);
  assert(f2.emit(NA, gPK[3], error_prob) == 0);

  assert(to!string(f2.emit(A, gPK[0], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2.emit(A, gPK[1], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2.emit(A, gPK[2], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2.emit(A, gPK[3], error_prob)) == to!string(log(error_prob/2.0)));

  assert(to!string(f2.emit(H, gPK[1], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2.emit(H, gPK[2], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2.emit(H, gPK[0], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2.emit(H, gPK[3], error_prob)) == to!string(log(error_prob/2.0)));

  assert(to!string(f2.emit(B, gPK[3], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2.emit(B, gPK[0], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2.emit(B, gPK[1], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2.emit(B, gPK[2], error_prob)) == to!string(log(error_prob/2.0)));

  assert(to!string(f2.emit(AorH, gPK[0], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2.emit(AorH, gPK[1], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2.emit(AorH, gPK[2], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2.emit(AorH, gPK[3], error_prob)) == to!string(log(error_prob)));

  assert(to!string(f2.emit(HorB, gPK[0], error_prob)) == to!string(log(error_prob)));
  assert(to!string(f2.emit(HorB, gPK[1], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2.emit(HorB, gPK[2], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2.emit(HorB, gPK[3], error_prob)) == to!string(log(1.0-error_prob/2.0)));

  // unit test nrec for F2 (phase-known)
  assert( f2.nrec(gPK[0], gPK[0]) == 0.0 );
  assert( f2.nrec(gPK[0], gPK[1]) == 0.5 );
  assert( f2.nrec(gPK[0], gPK[2]) == 0.5 );
  assert( f2.nrec(gPK[0], gPK[3]) == 1.0 );

  assert( f2.nrec(gPK[1], gPK[0]) == 0.5 );
  assert( f2.nrec(gPK[1], gPK[1]) == 0.0 );
  assert( f2.nrec(gPK[1], gPK[2]) == 1.0 );
  assert( f2.nrec(gPK[1], gPK[3]) == 0.5 );

  assert( f2.nrec(gPK[2], gPK[0]) == 0.5 );
  assert( f2.nrec(gPK[2], gPK[2]) == 0.0 );
  assert( f2.nrec(gPK[2], gPK[1]) == 1.0 );
  assert( f2.nrec(gPK[2], gPK[3]) == 0.5 );

  assert( f2.nrec(gPK[3], gPK[3]) == 0.0 );
  assert( f2.nrec(gPK[3], gPK[1]) == 0.5 );
  assert( f2.nrec(gPK[3], gPK[2]) == 0.5 );
  assert( f2.nrec(gPK[3], gPK[0]) == 1.0 );
}

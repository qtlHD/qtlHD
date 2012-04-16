/**
 * HMM intercross functions
 **/

module qtl.core.hmm.f2;

import std.math, std.stdio, std.path;
import std.exception;
import std.conv;

import qtl.core.primitives;
import qtl.core.genotype;
import qtl.core.hmm.calcgenoprob;
import qtl.core.hmm.estmap;
import qtl.core.map.genetic_map_functions;

// things for the unit tests
import qtl.plugins.qtab.read_qtab;

TrueGenotype[] allTrueGeno_F2()
{
  auto g1 = new TrueGenotype(0,0);
  auto g2 = new TrueGenotype(1,0);
  auto g3 = new TrueGenotype(1,1);
  return [g1, g2, g3];
}

// phase-known version (for est_map)
TrueGenotype[] allTrueGeno_F2PK()
{
  auto g1 = new TrueGenotype(0,0);
  auto g2 = new TrueGenotype(1,0);
  auto g3 = new TrueGenotype(0,1);
  auto g4 = new TrueGenotype(1,1);
  return [g1, g2, g3, g4];
}


unittest {
  writeln("Unit test " ~ __FILE__);

  writeln("    unit test allTrueGeno for F2");
  auto g = allTrueGeno_F2();

  assert(g.length == 3);
  assert(g[0] == new TrueGenotype(0,0));
  assert(g[1] == new TrueGenotype(1,0));
  assert(g[2] == new TrueGenotype(1,1));


  writeln("    unit test allTrueGeno for F2 (phase-known)");
  auto gPK = allTrueGeno_F2PK();

  assert(gPK.length == 4);
  assert(gPK[0] == new TrueGenotype(0,0));
  assert(gPK[1] == new TrueGenotype(1,0));
  assert(gPK[2] == new TrueGenotype(0,1));
  assert(gPK[3] == new TrueGenotype(1,1));
}

// ln Pr(true genotype)
double init_F2(in TrueGenotype truegen)
{
  auto g = allTrueGeno_F2();

  if(truegen==g[0] || truegen==g[2]) // AA or BB
    return(-2.0*LN2);

  if(truegen==g[1]) // AB
    return(-LN2);

  throw new Exception("truegen not among possible true genotypes");
}

// ln Pr(true genotype) for est_map
double init_F2PK(in TrueGenotype truegen)
{
  auto g = allTrueGeno_F2PK();

  if(truegen==g[0] || truegen==g[1] ||
     truegen==g[2] || truegen==g[3])
    return(-2.0*LN2);

  throw new Exception("truegen not among possible true genotypes");
}

unittest {
  writeln("    unit test init for F2");
  auto g = allTrueGeno_F2();
  double val[] = [log(0.25), log(0.5), log(0.25)];

  foreach(i, gv; g)
    assert(to!string(init_F2(gv)) == to!string(val[i]));

  writeln("    unit test init for F2 (phase-known)");
  auto gPK = allTrueGeno_F2PK();
  foreach(gv; gPK)
    assert(to!string(init_F2PK(gv)) == to!string(val[0]));
}

// ln Pr(genotype at right marker | genotype at left marker)
double step_F2(in TrueGenotype truegen_left, in TrueGenotype truegen_right, in double rec_frac)
{
  auto atg = allTrueGeno_F2();
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

// like step_F2, but with phase-known genotypes (for use with est_map)
double step_F2PK(in TrueGenotype truegen_left, in TrueGenotype truegen_right, in double rec_frac)
{
  auto atg = allTrueGeno_F2PK();
  if(truegen_left==atg[0]) {
    if(truegen_right==atg[0]) // AA -> AA
      return(2.0*log(1.0-rec_frac));
    if(truegen_right==atg[1] || truegen_right==atg[2]) // AA -> AB; AA -> BA
      return(log(1.0-rec_frac) + log(rec_frac));
    if(truegen_right==atg[3]) // AA -> BB
      return(2.0*log(rec_frac));
  }
  else if(truegen_left == atg[1]) {
    if(truegen_right==atg[0] || // AB -> AA
       truegen_right==atg[3]) // AB -> BB
      return(log(rec_frac) + log(1.0-rec_frac));
    if(truegen_right==atg[1]) // AB -> AB
      return(log((1.0-rec_frac)^^2));
    if(truegen_right==atg[2]) // AB -> BA
      return(log(rec_frac^^2));
  }
  else if(truegen_left == atg[2]) {
    if(truegen_right==atg[0] || // BA -> AA
       truegen_right==atg[3]) // BA -> BB
      return(log(rec_frac) + log(1.0-rec_frac));
    if(truegen_right==atg[2]) // BA -> BA
      return(log((1.0-rec_frac)^^2));
    if(truegen_right==atg[1]) // BA -> AB
      return(log(rec_frac^^2));
  }
  else if(truegen_left == atg[3]) {
    if(truegen_right==atg[0]) // BB -> AA
      return(2.0*log(rec_frac));
    if(truegen_right==atg[1] || // BB -> AB
       truegen_right==atg[2])   // BB -> BA
      return(log(1.0-rec_frac) + log(rec_frac));
    if(truegen_right==atg[3]) // BB -> BB
      return(2.0*log(1.0-rec_frac));
  }

  throw new Exception("inputs not among the possible true genotypes");
}


unittest {
  writeln("    unit test step for F2");
  auto g = allTrueGeno_F2();
  double rf = 0.01;

  assert( to!string( step_F2(g[0], g[0], rf) ) ==
          to!string( log((1-rf)^^2) ));
  assert( to!string( step_F2(g[0], g[1], rf) ) ==
          to!string( log(2.0*rf*(1-rf)) ));
  assert( to!string( step_F2(g[0], g[2], rf) ) ==
          to!string( log(rf^^2) ));

  assert( to!string( step_F2(g[1], g[0], rf) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( step_F2(g[1], g[1], rf) ) ==
          to!string( log(rf^^2 + (1-rf)^^2) ));
  assert( to!string( step_F2(g[1], g[2], rf) ) ==
          to!string( log(rf*(1-rf)) ));

  assert( to!string( step_F2(g[2], g[0], rf) ) ==
          to!string( log(rf^^2) ));
  assert( to!string( step_F2(g[2], g[1], rf) ) ==
          to!string( log(2.0*rf*(1-rf)) ));
  assert( to!string( step_F2(g[2], g[2], rf) ) ==
          to!string( log((1-rf)^^2) ));


  writeln("    unit test step for F2 (phase-known)");
  auto gPK = allTrueGeno_F2PK();

  assert( to!string( step_F2PK(gPK[0], gPK[0], rf) ) ==
          to!string( log((1-rf)^^2) ));
  assert( to!string( step_F2PK(gPK[0], gPK[1], rf) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( step_F2PK(gPK[0], gPK[2], rf) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( step_F2PK(gPK[0], gPK[3], rf) ) ==
          to!string( log(rf^^2) ));

  assert( to!string( step_F2PK(gPK[1], gPK[0], rf) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( step_F2PK(gPK[1], gPK[1], rf) ) ==
          to!string( log((1-rf)^^2) ));
  assert( to!string( step_F2PK(gPK[1], gPK[2], rf) ) ==
          to!string( log(rf^^2) ));
  assert( to!string( step_F2PK(gPK[1], gPK[3], rf) ) ==
          to!string( log(rf*(1.0-rf)) ));

  assert( to!string( step_F2PK(gPK[2], gPK[0], rf) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( step_F2PK(gPK[2], gPK[2], rf) ) ==
          to!string( log((1-rf)^^2) ));
  assert( to!string( step_F2PK(gPK[2], gPK[1], rf) ) ==
          to!string( log(rf^^2) ));
  assert( to!string( step_F2PK(gPK[2], gPK[3], rf) ) ==
          to!string( log(rf*(1.0-rf)) ));

  assert( to!string( step_F2PK(gPK[3], gPK[3], rf) ) ==
          to!string( log((1-rf)^^2) ));
  assert( to!string( step_F2PK(gPK[3], gPK[1], rf) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( step_F2PK(gPK[3], gPK[2], rf) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( step_F2PK(gPK[3], gPK[0], rf) ) ==
          to!string( log(rf^^2) ));
}


// ln Pr(observed genotype | true genotype)
double emit_F2(in GenotypeCombinator obsgen, in TrueGenotype truegen, in double error_prob)
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

alias emit_F2 emit_F2PK;

unittest {
  writeln("    unit test emit for F2");
  auto atg = allTrueGeno_F2();
  auto atgPK = allTrueGeno_F2PK();

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

  assert(emit_F2(NA, atg[0], error_prob) == 0);
  assert(emit_F2(NA, atg[1], error_prob) == 0);
  assert(emit_F2(NA, atg[2], error_prob) == 0);

  assert(to!string(emit_F2(A, atg[0], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(emit_F2(A, atg[1], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(emit_F2(A, atg[2], error_prob)) == to!string(log(error_prob/2.0)));

  assert(to!string(emit_F2(H, atg[1], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(emit_F2(H, atg[0], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(emit_F2(H, atg[2], error_prob)) == to!string(log(error_prob/2.0)));

  assert(to!string(emit_F2(B, atg[2], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(emit_F2(B, atg[0], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(emit_F2(B, atg[1], error_prob)) == to!string(log(error_prob/2.0)));

  assert(to!string(emit_F2(AorH, atg[0], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(emit_F2(AorH, atg[1], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(emit_F2(AorH, atg[2], error_prob)) == to!string(log(error_prob)));

  assert(to!string(emit_F2(HorB, atg[0], error_prob)) == to!string(log(error_prob)));
  assert(to!string(emit_F2(HorB, atg[1], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(emit_F2(HorB, atg[2], error_prob)) == to!string(log(1.0-error_prob/2.0)));

  assert(emit_F2PK(NA, atgPK[0], error_prob) == 0);
  assert(emit_F2PK(NA, atgPK[1], error_prob) == 0);
  assert(emit_F2PK(NA, atgPK[2], error_prob) == 0);
  assert(emit_F2PK(NA, atgPK[3], error_prob) == 0);

  assert(to!string(emit_F2PK(A, atgPK[0], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(emit_F2PK(A, atgPK[1], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(emit_F2PK(A, atgPK[2], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(emit_F2PK(A, atgPK[3], error_prob)) == to!string(log(error_prob/2.0)));

  assert(to!string(emit_F2PK(H, atgPK[1], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(emit_F2PK(H, atgPK[2], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(emit_F2PK(H, atgPK[0], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(emit_F2PK(H, atgPK[3], error_prob)) == to!string(log(error_prob/2.0)));

  assert(to!string(emit_F2PK(B, atgPK[3], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(emit_F2PK(B, atgPK[0], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(emit_F2PK(B, atgPK[1], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(emit_F2PK(B, atgPK[2], error_prob)) == to!string(log(error_prob/2.0)));

  assert(to!string(emit_F2PK(AorH, atgPK[0], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(emit_F2PK(AorH, atgPK[1], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(emit_F2PK(AorH, atgPK[2], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(emit_F2PK(AorH, atgPK[3], error_prob)) == to!string(log(error_prob)));

  assert(to!string(emit_F2PK(HorB, atgPK[0], error_prob)) == to!string(log(error_prob)));
  assert(to!string(emit_F2PK(HorB, atgPK[1], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(emit_F2PK(HorB, atgPK[2], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(emit_F2PK(HorB, atgPK[3], error_prob)) == to!string(log(1.0-error_prob/2.0)));
}

// proportion of recombination events
double nrec_F2PK(in TrueGenotype truegen_left, in TrueGenotype truegen_right)
{
  auto atg = allTrueGeno_F2PK();
  if(truegen_left==atg[0]) {
    if(truegen_right==atg[0]) // AA -> AA
      return(0.0);
    if(truegen_right==atg[1] || truegen_right==atg[2]) // AA -> AB; AA -> BA
      return(0.5);
    if(truegen_right==atg[3]) // AA -> BB
      return(1.0);
  }
  else if(truegen_left == atg[1]) {
    if(truegen_right==atg[0] || // AB -> AA
       truegen_right==atg[3]) // AB -> BB
      return(0.5);
    if(truegen_right==atg[1]) // AB -> AB
      return(0.0);
    if(truegen_right==atg[2]) // AB -> BA
      return(1.0);
  }
  else if(truegen_left == atg[2]) {
    if(truegen_right==atg[0] || // BA -> AA
       truegen_right==atg[3]) // BA -> BB
      return(0.5);
    if(truegen_right==atg[2]) // BA -> BA
      return(0.0);
    if(truegen_right==atg[1]) // BA -> AB
      return(1.0);
  }
  else if(truegen_left == atg[3]) {
    if(truegen_right==atg[0]) // BB -> AA
      return(1.0);
    if(truegen_right==atg[1] || // BB -> AB
       truegen_right==atg[2])   // BB -> BA
      return(0.5);
    if(truegen_right==atg[3]) // BB -> BB
      return(0.0);
  }

  throw new Exception("inputs not among the possible true genotypes");
}


unittest {
  writeln("    unit test nrec for F2 (phase-known)");
  auto gPK = allTrueGeno_F2PK();

  assert( nrec_F2PK(gPK[0], gPK[0]) == 0.0 );
  assert( nrec_F2PK(gPK[0], gPK[1]) == 0.5 );
  assert( nrec_F2PK(gPK[0], gPK[2]) == 0.5 );
  assert( nrec_F2PK(gPK[0], gPK[3]) == 1.0 );

  assert( nrec_F2PK(gPK[1], gPK[0]) == 0.5 );
  assert( nrec_F2PK(gPK[1], gPK[1]) == 0.0 );
  assert( nrec_F2PK(gPK[1], gPK[2]) == 1.0 );
  assert( nrec_F2PK(gPK[1], gPK[3]) == 0.5 );

  assert( nrec_F2PK(gPK[2], gPK[0]) == 0.5 );
  assert( nrec_F2PK(gPK[2], gPK[2]) == 0.0 );
  assert( nrec_F2PK(gPK[2], gPK[1]) == 1.0 );
  assert( nrec_F2PK(gPK[2], gPK[3]) == 0.5 );

  assert( nrec_F2PK(gPK[3], gPK[3]) == 0.0 );
  assert( nrec_F2PK(gPK[3], gPK[1]) == 0.5 );
  assert( nrec_F2PK(gPK[3], gPK[2]) == 0.5 );
  assert( nrec_F2PK(gPK[3], gPK[0]) == 1.0 );
}

double[][][] calc_geno_prob_F2(in GenotypeCombinator[][] genotypes,
                               in Marker[] marker_map,
                               in double[] rec_frac,
                               in double error_prob)
{
  auto all_true_geno = allTrueGeno_F2();

  return( calc_geno_prob!(init_F2, emit_F2, step_F2)(genotypes, all_true_geno,
                                                     marker_map, rec_frac,
                                                     error_prob) );
}


double[] estmap_F2(in GenotypeCombinator[][] genotypes, 
                   in Marker[] marker_map,
                   in double[] rec_frac, 
                   in double error_prob, 
                   in int max_iterations,
                   in double tol, 
                   in bool verbose)
{
  auto all_true_geno = allTrueGeno_F2PK();

  return( estmap!(init_F2PK, emit_F2PK, step_F2PK, nrec_F2PK)(genotypes, all_true_geno,
                                                              marker_map, rec_frac,
                                                              error_prob, max_iterations,
                                                              tol, verbose) );
}

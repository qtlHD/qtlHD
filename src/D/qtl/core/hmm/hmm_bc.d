/**
 * hmm_bc
 **/

module qtl.core.hmm.hmm_bc;

import std.math, std.stdio, std.path;
import std.exception;
import std.conv;

import qtl.core.primitives;
import qtl.core.genotype;
import qtl.core.hmm.hmm_calcgenoprob;
import qtl.core.hmm.hmm_estmap;
import qtl.core.map.genetic_map_functions;

// things for the unit tests
import qtl.plugins.qtab.read_qtab;

// vector of true genotypes
TrueGenotype[] allTrueGeno_BC()
{
  auto g1 = new TrueGenotype(0,0);
  auto g2 = new TrueGenotype(1,0);
  return [g1, g2];
}

unittest {
  writeln("Unit test " ~ __FILE__);

  writeln("    unit test allTrueGeno for BC");
  auto g = allTrueGeno_BC();

  assert(g.length == 2);
  assert(g[0] == new TrueGenotype(0,0));
  assert(g[1] == new TrueGenotype(1,0));
}

// ln Pr(true genotype)
double init_BC(in TrueGenotype truegen)
{
  if(truegen != new TrueGenotype(0,0) &&
     truegen != new TrueGenotype(1,0))
    throw new Exception("truegen must be 0,0 or 1,0");

  return(-LN2);
}

unittest {
  writeln("    unit test init for BC");

  auto g = allTrueGeno_BC();
  foreach(gv; g)
    assert(to!string(init_BC(gv)) == to!string(log(0.5)));
}

// ln Pr(genotype at right marker | genotype at left marker)
double step_BC(in TrueGenotype truegen_left, in TrueGenotype truegen_right, in double rec_frac)
{
  auto atg = allTrueGeno_BC();

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

unittest {
  writeln("    unit test step for BC");
  auto g = allTrueGeno_BC();
  double rf = 0.01;

  assert( to!string( step_BC(g[0], g[0], rf) ) ==
          to!string( log((1-rf)) ));
  assert( to!string( step_BC(g[0], g[1], rf) ) ==
          to!string( log(rf) ));

  assert( to!string( step_BC(g[1], g[0], rf) ) ==
          to!string( log(rf) ));
  assert( to!string( step_BC(g[1], g[1], rf) ) ==
          to!string( log((1-rf)) ));
}

// ln Pr(observed genotype | true genotype)
double emit_BC(in GenotypeCombinator obsgen, in TrueGenotype truegen, in double error_prob)
{
  if(obsgen.list.length==0) // missing value
    return(0.0); // log(1.0)

  if(obsgen.match(truegen)) // compatible with truegen?
    return(log(1.0-error_prob));
  else
    return(log(error_prob));
}

unittest {
  writeln("    unit test emit for BC");
  auto atg = allTrueGeno_BC();

  auto NA = new GenotypeCombinator("-");
  auto A = new GenotypeCombinator("A");
  A ~= new TrueGenotype(0,0);
  auto H = new GenotypeCombinator("H");
  H ~= new TrueGenotype(1,0);

  double error_prob = 0.01;

  assert(emit_BC(NA, atg[0], error_prob) == 0);
  assert(emit_BC(NA, atg[1], error_prob) == 0);
  assert(to!string(emit_BC(A, atg[0], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(emit_BC(A, atg[1], error_prob)) == to!string(log(error_prob)));
  assert(to!string(emit_BC(H, atg[1], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(emit_BC(H, atg[0], error_prob)) == to!string(log(error_prob)));
}

// No. recombination events
double nrec_BC(in TrueGenotype truegen_left, in TrueGenotype truegen_right)
{
  auto atg = allTrueGeno_BC();

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

unittest {
  writeln("    unit test nrec for BC");
  auto g = allTrueGeno_BC();

  assert(nrec_BC(g[0], g[0]) == 0.0);
  assert(nrec_BC(g[0], g[1]) == 1.0);

  assert(nrec_BC(g[1], g[0]) == 1.0);
  assert(nrec_BC(g[1], g[1]) == 0.0);
}

Probability[][][] calc_geno_prob_BC(in GenotypeCombinator[][] genotypes,
                               in Marker[] marker_map,
                               in double[] rec_frac,
                               in double error_prob)
{
  auto all_true_geno = allTrueGeno_BC();

  return( calc_geno_prob!(init_BC, emit_BC, step_BC)(genotypes, all_true_geno,
                                                     marker_map, rec_frac,
                                                     error_prob) );
}


double[] estmap_BC(in GenotypeCombinator[][] genotypes, 
                   in Marker[] marker_map,
                   in double[] rec_frac, 
                   in double error_prob, 
                   in int max_iterations,
                   in double tol, 
                   in bool verbose)
{
  auto all_true_geno = allTrueGeno_BC();

  return( estmap!(init_BC, emit_BC, step_BC, nrec_BC)(genotypes, all_true_geno,
                                                      marker_map, rec_frac,
                                                      error_prob, max_iterations,
                                                      tol, verbose) );
}

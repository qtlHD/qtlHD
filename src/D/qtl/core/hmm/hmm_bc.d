/**
 * hmm_bc
 **/

module qtl.core.hmm.hmm_bc;

import std.math, std.stdio, std.path;
import std.exception;
import std.conv;

import qtl.core.primitives;
import qtl.core.genotype;
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
double emit_BC(Gref obsgen, TrueGenotype truegen, double error_prob)
{
  return 0.0;
}

unittest {
  writeln("    unit test emit for BC");
  auto atg = allTrueGeno_BC();
  ObservedGenotypes aog[3];
}
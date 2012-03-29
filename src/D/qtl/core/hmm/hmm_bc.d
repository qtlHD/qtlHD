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

  auto g = allTrueGeno_BC();

  assert(g.length == 2);
  assert(g[0] == new TrueGenotype(0,0));
  assert(g[1] == new TrueGenotype(1,0));
}

// Pr(true genotype)
double init_BC(TrueGenotype truegen)
{
  if(truegen != new TrueGenotype(0,0) &&
     truegen != new TrueGenotype(1,0))
    throw new Exception("truegen must be 0,0 or 1,0");

  return(-LN2);
}

unittest {
  auto g = allTrueGeno_BC();
  foreach(gv; g)
    assert(to!string(init_BC(gv)) == to!string(log(0.5)));
}

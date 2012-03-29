/**
 * hmm_f2
 **/

module qtl.core.hmm.hmm_f2;

import std.math, std.stdio, std.path;
import std.exception;
import std.conv;

import qtl.core.primitives;
import qtl.core.genotype;
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

// phase-known version (for est.map)
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

  auto g = allTrueGeno_F2();
  auto gPK = allTrueGeno_F2PK();

  assert(g.length == 3);
  assert(g[0] == new TrueGenotype(0,0));
  assert(g[1] == new TrueGenotype(1,0));
  assert(g[2] == new TrueGenotype(1,1));

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

unittest {
  auto g = allTrueGeno_F2();
  double val[] = [log(0.25), log(0.5), log(0.25)];

  foreach(i, gv; g)
    assert(to!string(init_F2(gv)) == to!string(val[i]));
}


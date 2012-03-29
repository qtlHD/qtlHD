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
  auto gPK = allTrueGeno_F2PK();
  double val[] = [log(0.25), log(0.5), log(0.25)];

  foreach(i, gv; g)
    assert(to!string(init_F2(gv)) == to!string(val[i]));

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
}

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

TrueGenotype[] allTrueGenoBC()
{
  auto g1 = new TrueGenotype("0,0");
  auto g2 = new TrueGenotype("1,0");
  return [g1, g2];
}

unittest {
  writeln("Unit test " ~ __FILE__);
  
  auto g = allTrueGenoBC();
  writeln(g.length);
  foreach(gv; g) writeln(gv);
}
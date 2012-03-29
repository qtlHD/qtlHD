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

TrueGenotype[] allTrueGenoF2()
{
  auto g1 = new TrueGenotype(0,0);
  auto g2 = new TrueGenotype(1,0);
  auto g3 = new TrueGenotype(1,1);
  return [g1, g2, g3];
}

unittest {
  writeln("Unit test " ~ __FILE__);
  
  auto g = allTrueGenoF2();
  writeln(g.length);
  foreach(gv; g) writeln(gv);
}
/**
 * Test mqm routines, using hyper_noX (CSV) set
 */

module test.mqm.test_mqmscan;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.core.marker;
import qtl.core.map;
import qtl.core.make_map;
import qtl.plugins.input.read_csv;
import qtl.core.scanone_hk;
import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.algorithm;

alias std.algorithm.find find;

static bool VERBOSE = false;

unittest{
  writeln("Unit test " ~ __FILE__);
  alias std.path.join join;
  auto fn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","test","data","input","hyper_noX.csv");
  if(VERBOSE) writeln("  - reading CSV " ~ fn);
}
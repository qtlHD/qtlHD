/**
 * Test mqm routines, using hyper_noX (CSV) set
 */

module test.mqm.test_mqmscan;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.core.individual;
import example.genotype_examples;
import qtl.core.mqm.mqm_types;
import qtl.core.marker;
import qtl.core.map.map;
import qtl.core.mqm.matrix;
import qtl.core.map.make_map;
import qtl.plugins.csv.read_csv;
import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.algorithm;

alias std.algorithm.find find;

unittest{
  writeln("Unit test " ~ __FILE__);
  alias std.path.buildPath buildPath;
  auto fn = dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","test","data","input","hyper_noX.csv");
  if(VERBOSE) writeln("  - reading CSV " ~ fn);
  auto indata = new ReadSimpleCSV!(F2,ObservedF2)(fn);

}

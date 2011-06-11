/**
 * Chromosome related
 */

module qtl.core.chromosome;

import std.container;
import qtl.core.primitives;

import std.stdio;
unittest {
  writeln("Unit test " ~ __FILE__);
  // For chromosome, for now, simply use an associative array
  Chromosome[string] chromosomes;
  immutable one = "1";
  auto c1 = new Chromosome(one);
  chromosomes[one] = c1;
  assert(c1.id == 1);
  assert(c1.name == "1");
  immutable X = "X";
  auto cx = new Chromosome(X);
  chromosomes[X] = cx;
  assert(cx.id == 0);
  assert(cx.name == "X");
}

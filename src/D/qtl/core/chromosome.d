/**
 * Chromosome related
 */

module qtl.core.chromosome;

import std.container;
import std.conv;
import qtl.core.primitives;

import std.stdio;

/**
 * Create a new Chromosome object, basing the id on the content of name
 */

Chromosome get_chromosome_with_id(string name) {
  uint id;
  if (name == "X") id = 0;  
  else             id = to!int(name);

  return new Chromosome(name,id);
}

unittest {
  writeln("Unit test " ~ __FILE__);
  // For chromosome, for now, simply use an associative array
  Chromosome[string] chromosomes;
  immutable one = "1";
  auto c1 = new Chromosome(one,1);
  chromosomes[one] = c1;
  assert(c1.id == 1);
  assert(c1.name == "1");
  immutable X = "X";
  auto cx = new Chromosome(X,0);
  chromosomes[X] = cx;
  assert(cx.id == 0);
  assert(cx.name == "X");
}

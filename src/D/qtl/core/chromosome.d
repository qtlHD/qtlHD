/**
 * Chromosome related
 */

module qtl.core.chromosome;

import std.container;
import std.conv;
import qtl.core.primitives;

import std.stdio;

/**
 * Create a new Chromosome object, basing the id on the content of name.
 * Currently string 'X' returns a SexChromosome.
 * (this may change)
 */

Chromosome get_chromosome_with_id(string name) {
  uint id;
  if (name == "X") id = 0;  
  else             id = to!int(name);
  return get_chromosome(name,id);
}

/**
 * Create new Chromosome object. Currently, if id is zero it 
 * returns a SexChromosome, otherwise an Autosome.
 * (this may change)
 */

Chromosome get_chromosome(string name, uint id) {
  if (id==0) 
    return new SexChromosome(name);
  else
    return new Autosome(name,id);
}

/*
 * Test for chromosome sex - FIXME should be done on type instead
 */

bool is_sex(Chromosome chromosome) {
  return (chromosome.id == 0);
}

unittest {
  writeln("Unit test " ~ __FILE__);
  // For chromosome, for now, simply use an associative array
  Chromosome[string] chromosomes;
  immutable one = "1";
  auto c1 = get_chromosome_with_id(one);
  auto c12 = get_chromosome(one,2);
  chromosomes[one] = c1;
  assert(c1.id == 1);
  assert(c1.name == "1");
  assert(c12.id == 2);
  immutable X = "X";
  auto cx = get_chromosome_with_id(X);
  chromosomes[X] = cx;
  assert(cx.id == 0);
  assert(cx.name == "X");

  // test for sex
  assert(!is_sex(c1));
  assert(is_sex(cx),typeof(cx).stringof ~ to!string(cx.id));
}

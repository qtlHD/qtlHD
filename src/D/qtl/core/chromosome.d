/**
 * Chromosome related functions
 */

module qtl.core.chromosome;

import std.container;
import std.conv;
import std.variant;
import qtl.core.primitives;

import std.stdio;
// import std.container;
import std.typecons; 


/**
 * Create a new Chromosome object, basing the id on the content of name.
 * Currently string 'X' returns a SexChromosome.
 * (this may change)
 */

Chromosome get_chromosome_with_id(string name) {
  uint id;
  if (name == "X") id = 0;  
  else             id = to!int(name);
  return get_chromosome(name,id,(name == "X"));
}

/**
 * Create new Chromosome object. Currently, if is_sex is true it 
 * returns a SexChromosome, otherwise an Autosome.
 */

Chromosome get_chromosome(string name, uint id, bool is_sex=false) {
  if (is_sex) 
    return new SexChromosome(name);
  else
    return new Autosome(name,id);
}

/**
 * Test for chromosome sex - based on object type
 */

static bool is_sex(Chromosome chromosome) {
  return (typeid(chromosome) == typeid(SexChromosome));
}

/**
 * Take a full list of markers and return a Tuple list of
 * (chromosomes,Markers), i.e. Chromosomes with their markers attached.
 */


Tuple!(Chromosome,Ms)[] chromosome_markers(Ms)(in Ms markers) {
  Ms[string] alist;
  foreach(m ; markers.list) {
    if ((m.chromosome.name) !in alist) {
      alist[m.chromosome.name] = new Ms();
    }
    alist[m.chromosome.name].list ~= new Marker(m); // FIXME: what about other types?
  }
  // convert to ret type
  Tuple!(Chromosome, Ms)[] list; 
  foreach( cname, ms ; alist) {
     list ~= Tuple!(Chromosome, Ms)(ms.list[0].chromosome,ms);
  }
  return list;
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
  assert(cx.id == ID_UNKNOWN);
  assert(cx.name == "X");

  // test for sex
  assert(is_sex(cx),typeof(cx).stringof ~ to!string(cx.id) ~ to!string(is_sex(cx)));
  assert(!is_sex(c1));

  auto markers = new Markers!Marker();
  auto m1 = new Marker(c1,10.0);
  auto m2 = new Marker(c1,20.0);
  auto m3 = new Marker(c1,30.0);
  auto m4 = new Marker(cx,11.0);
  auto m5 = new Marker(cx,14.0);
  markers.list ~= m1;
  markers.list ~= m2;
  markers.list ~= m3;
  markers.list ~= m4;
  markers.list ~= m5;
  // fetch markers by Chromosome
  auto tlist = chromosome_markers(markers);
  foreach(c, ms ; tlist) {
    writeln(c,ms);
  }
  assert(tlist[0][0].name == '1');
}

/**
 * Chromosome related functions
 */

module qtl.core.chromosome;

import std.container;
import std.conv;
import std.variant;
import std.random;
import qtl.core.primitives;
import std.algorithm; // to use sort()

import std.stdio;
// import std.container;
import std.typecons;


/**
 * Create a new Chromosome object, basing the id on the content of name.
 * Currently string 'X' returns an Xchromosome.
 * (this may change)
 */

Chromosome get_chromosome_with_id(string name) {
  uint id;
  if (name == "X") id = 0;
  else             id = to!int(name);
  return get_chromosome(name,id,(name == "X"));
}

/**
 * Create new Chromosome object. Currently, if is_X_chr is true it
 * returns a Xchromosome, otherwise an Autosome.
 */

Chromosome get_chromosome(string name, uint id, bool is_X_chr=false) {
  if (is_X_chr)
    return new Xchromosome(name);
  else
    return new Autosome(name,id);
}

/**
 * Test for chromosome X - based on object type
 */

static bool is_X_chr(Chromosome chromosome) {
  return (typeid(chromosome) == typeid(Xchromosome));
}

/**
 * Take a full list of markers and return a Tuple list of
 * (chromosomes,Markers), i.e. Chromosomes with their markers (unsorted). 
 *
 * Note conversion does not affect derived types, such as PseudoMarker.
 */

Tuple!(Chromosome,Ms)[] get_markers_by_chromosome(Ms)(in Ms markers) {
  Ms[Chromosome] c_markers;   // store markers temporarily by chromosome
  foreach(m ; markers.list) {
    c_markers[m.chromosome] ~= cast(Marker)m;
  }
  // ---- convert to ret type
  Tuple!(Chromosome, Ms)[] list;
  foreach( cname, ms ; c_markers) {
    list ~= tuple(ms[0].chromosome,ms);
  }
  return list;
}

static auto VERBOSE = false;

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

  // test for X chromosome
  assert(is_X_chr(cx),typeof(cx).stringof ~ to!string(cx.id) ~ to!string(is_X_chr(cx)));
  assert(!is_X_chr(c1));

  /*
  auto markers = new Markers!Marker();
  auto m1 = new Marker(c1,10.0);
  auto m2 = new Marker(c1,20.0);
  auto m3 = new Marker(c1,30.0);
  auto m4 = new Marker(cx,11.0);
  auto m5 = new Marker(cx,14.0);
  auto m6 = new PseudoMarker(cx,12.0);
  markers.list ~= m1;
  markers.list ~= m2;
  markers.list ~= m3;
  markers.list ~= m4;
  markers.list ~= m5;
  markers.list ~= m6;
  // fetch markers by Chromosome
  auto c_mslist = get_markers_by_chromosome(markers);
  assert(c_mslist.length == 2);
  if (VERBOSE) {
    foreach(c, ms ; c_mslist) {
       writeln(c,ms);
       foreach (m; ms[1].list) {
         writeln(typeid(m));
       }
    }
  }
  auto chromosome1 = c_mslist[0][0];
  assert(chromosome1.name == "X");
  */
}


// comparison function for sorting
static bool compare_chromosomes_by_marker_id(Tup)(in Tup a, in Tup b)
{
  return a[1][0].id <= b[1][0].id;
}


// sort the chromosomes in the results of get_markers_by_chromosomes
Tuple!(Chromosome,Ms)[] sort_chromosomes_by_marker_id(Ms)(Tuple!(Chromosome,Ms)[] chr_tuples)
{
  // need a copy of chr_tuples, so it doesn't get sorted
  auto output = new Tuple!(Chromosome,Ms)[chr_tuples.length];
  foreach(i, m; chr_tuples)
    output[i] = m;

  sort!(compare_chromosomes_by_marker_id)(output);
  return(output);
}

unittest {
  writeln("test sort_chromosome_by_marker_id");
  Marker markers[];
  int i, j;
  uint k;
  for(j=0, k=0; j<5; j++) {
    for(i=0; i<5; i++, k++) {
      auto chr = new Chromosome(to!string(j+1));
      auto marker = new Marker(chr, i*5.0, "marker" ~ to!string(k+1), k);
      markers ~= marker;
    }
  }

  // markers split by chromosome and then sorted
  auto markers_by_chr = get_markers_by_chromosome(markers);
  auto markers_by_chr_sorted = sort_chromosomes_by_marker_id(markers_by_chr);
}

/**
 * ChromosomeMap combines Chromosome and Marker list.
 */

/**
 * Disabled fow now
class ChromosomeMap(T) {
  Chromosome chromosome;
  // Markers!(MarkerRef!T) markers;
}
*/

/**
 * The FullMap has an ordered chromosome map.

class FullMap(T) {
  SList!(ChromosomeMap!T) chromosome_map;
}
unittest {
  // e.g. for F2 auto map = new FullMap!F2();
  // this should also compile:
  auto map = new FullMap!uint();
  foreach ( c ; map.chromosome_map ) {
    auto markers = c.markers;
    foreach ( m ; markers.list ) {
    }
  }
}

*/


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


Tuple!(Chromosome,Ms)[] get_markers_by_chromosome(Ms)(in Ms markers) {
  Ms[string] hash_cs;   // store chromosomes temporarily in associative array hash_cs
  auto ms = markers.list;
  foreach(m ; ms) {
    // note the cast does not affect derived types, such as PseudoMarker
    static if (is(Ms : Object)) {
      if ((m.chromosome.name) !in hash_cs)
        hash_cs[m.chromosome.name] = new Ms();
      hash_cs[m.chromosome.name].list ~= cast(Marker)m;
    }
    else {
      hash_cs[m.chromosome.name] ~= cast(Marker)m;
    }
  }
  // convert to ret type
  Tuple!(Chromosome, Ms)[] list;
  foreach( cname, ms ; hash_cs) {
    list ~= Tuple!(Chromosome, Ms)(ms[0].chromosome,ms);
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

  // test for sex
  assert(is_sex(cx),typeof(cx).stringof ~ to!string(cx.id) ~ to!string(is_sex(cx)));
  assert(!is_sex(c1));

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

  writeln("  Markers split by chromosome:");
  auto markers_by_chr = get_markers_by_chromosome(markers);
  foreach(chr; markers_by_chr)
    writefln("%2s (%2d): %-9s (%3d)", chr[0].name, chr[1].length, chr[1][0].name, chr[1][0].id);

  writeln("\n  Chromosomes sorted:");
  auto markers_by_chr_sorted = sort_chromosomes_by_marker_id(markers_by_chr);
  foreach(chr; markers_by_chr_sorted)
    writefln("%2s (%2d): %-9s (%3d)", chr[0].name, chr[1].length, chr[1][0].name, chr[1][0].id);
  writeln();

  writeln("Has the original marker_by_chr changed?");
  foreach(chr; markers_by_chr)
    writefln("%2s (%2d): %-9s (%3d)", chr[0].name, chr[1].length, chr[1][0].name, chr[1][0].id);
  writeln();

  // it's the same stuff; just pointed to in a different order
  for(i=0; i<5; i++)
    for(j=0; j<5; j++)
      if(markers_by_chr[i] is markers_by_chr_sorted[j])
	writeln(i, " is ", j);
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


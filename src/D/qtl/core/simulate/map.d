/*
 Code to generate a marker map
*/

module qtl.core.simulate.map;

import qtl.core.primitives;
import qtl.core.chromosome;

import math.distributions.random;
import std.algorithm;
import std.stdio;
import std.string;
import std.conv;
import std.random;
import std.math;

// generate map of equally-spaced markers for one chromosome
// if first_marker_number=5, markers will be named like "m5", "m6", ...
Marker[] generate_map_eqspacing_onechr(in double chrlen, in uint n_markers, in uint first_marker_number, 
                                       Chromosome chromosome)
in {
  assert(chrlen >= 0, "chrlen must be >= 0");
  assert(n_markers > 0, "n_markers must be > 0");
  assert(chrlen > 0 || n_markers == 1, "If n_markers > 1, chrlen must be > 0");
}
body {
  Marker[] map;

  foreach(i; 0..n_markers) {
    double location = i*chrlen/(n_markers-1);

    Marker m = new Marker(location, "m" ~ to!string(first_marker_number + i), first_marker_number+i);
    m.chromosome = chromosome;
    map ~= m;
  }
  return map;
}

unittest {
  writeln("Unit test " ~ __FILE__);

  writeln("Equal spacing:");
  auto chromosome = get_chromosome_with_id("3");
  auto marker_map = generate_map_eqspacing_onechr(150.0, 4, 21, chromosome);

  foreach(m; marker_map) {
    writefln("%s %3s %3d %6.2f", m.chromosome.name, m.name, m.id, m.get_position);
  }
  writeln;
}


// generate map of randomly spaced markers for one chromosome
Marker[] generate_map_random_onechr(in double chrlen, in uint n_markers, in uint first_marker_number,
                                    Chromosome chromosome, in bool anchor_tel, ref Random gen)
in {
  assert(chrlen >= 0, "chrlen must be >= 0");
  assert(n_markers > 0, "n_markers must be > 0");
  assert(chrlen > 0 || n_markers == 1, "If n_markers > 1, chrlen must be > 0");
  assert(n_markers >= 2 || !anchor_tel, "If anchor_tel, n_markers must be >= 2");
}
body {
  Marker[] map;
  Marker m;

  uint n_markers_local = n_markers;
  uint first_marker_number_local = first_marker_number;

  if(anchor_tel) { // place markers at 0.0 and chrlen
    m = new Marker(0.0, "m" ~ to!string(first_marker_number_local), first_marker_number_local);
    m.chromosome = chromosome;
    map ~= m;

    int marker_number = first_marker_number_local + n_markers_local - 1;
    m = new Marker(chrlen, "m" ~ to!string(marker_number), marker_number);
    m.chromosome = chromosome;
    map ~= m;

    if(n_markers_local == 2) return(map);

    n_markers_local -= 2;
    first_marker_number_local++;
  }

  // simulate marker locations
  auto locations = new double[](n_markers_local);
  foreach(i; 0..n_markers_local) {
    locations[i] = uniform(0.0, chrlen, gen);
  }
  sort(locations);

  foreach(i; 0..n_markers_local) {
    m = new Marker(locations[i], "m" ~ to!string(first_marker_number_local+i), first_marker_number_local+i);
    m.chromosome = chromosome;
    map ~= m;
  }
  
  sort(map);
  return(map);
}

unittest {
  Random gen;
  gen.seed(unpredictableSeed);

  writeln("Random spacing; anchored telomeres:");
  auto chromosome = get_chromosome_with_id("7");
  auto marker_map = generate_map_random_onechr(150.0, 16, 31, chromosome, true, gen);

  foreach(m; marker_map) {
    writefln("%s %3s %3d %6.2f", m.chromosome.name, m.name, m.id, m.get_position);
  }
  writeln;

  writeln("Random spacing; telomeres not anchored:");
  marker_map = generate_map_random_onechr(200.0, 10, 51, chromosome, false, gen);

  foreach(m; marker_map) {
    writefln("%s %3s %3d %6.2f", m.chromosome.name, m.name, m.id, m.get_position);
  }
  writeln;
}

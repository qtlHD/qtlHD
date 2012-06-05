/**
 * Find max LOD score from scanone
 */

module qtl.core.scanone.peaks;

import std.stdio;
import std.algorithm;
import std.container;
import std.variant;
import std.conv;
import std.typecons;
import std.random;

import qtl.core.chromosome;
import qtl.core.primitives;

// find maximum lod score and the position at which it occurred
// if multiple positions share maximum, return a random one
Tuple!(double, Marker) get_peak_scanone(double[] lod, Marker[] map, ref Random gen)
in {
  assert(lod.length > 0, "lod should have length > 0");
  assert(lod.length == map.length, "lod and map should have the same length");
}
body {
  if(lod.length==1)
    return Tuple!(double, Marker)(lod[0], map[0]);

  double maxlod = lod[0];
  size_t[] index_maxlod = [0];

  foreach(i; 1..lod.length) {
    if(lod[i] == maxlod)
      index_maxlod ~= i;
    if(lod[i] > maxlod) {
      index_maxlod = [i];
      maxlod = lod[i];
    }
  }

  if(index_maxlod.length > 1)
    index_maxlod[0] = index_maxlod[uniform(0, index_maxlod.length, gen)];

  return Tuple!(double, Marker)(lod[index_maxlod[0]], map[index_maxlod[0]]);
}

// same with self-initiated random seed
Tuple!(double, Marker) get_peak_scanone(double[] lod, Marker[] map)
{
  Random gen;
  gen.seed(unpredictableSeed);
  return get_peak_scanone(lod, map, gen);
}


// same for doubly-index LOD scores
Tuple!(double, Marker)[] get_peak_scanone(double[][] lod, Marker[] map, ref Random gen)
in {
  assert(lod.length == map.length, "lod and map should have the same length");
}
body {
  Tuple!(double, Marker)[] peaks;

  foreach(i; 0..lod[0].length) {
    auto thislod = new double[](lod.length);
    foreach(j; 0..lod.length)
      thislod[j] = lod[j][i];

    peaks ~= get_peak_scanone(thislod, map, gen);
  }

  return peaks;
}

// same with self-initiated random seed
Tuple!(double, Marker)[] get_peak_scanone(double[][] lod, Marker[] map)
{
  Random gen;
  gen.seed(unpredictableSeed);
  return get_peak_scanone(lod, map, gen);
}

unittest {
  writeln("Unit test " ~ __FILE__);

  import qtl.core.simulate.map;
  import qtl.core.chromosome;

  // singly indexed lod
  double[] lod = [1.0, 3.0, 2.3, 5.6, 1.8];

  auto chromosome = get_chromosome_with_id("1");
  auto map = generate_map_eqspacing_onechr(100.0, cast(uint)lod.length, 1, chromosome);

  auto peak = get_peak_scanone(lod, map);
  assert(peak[0] == 5.6 && peak[1].name == "m4");

  // doubly indexed lod
  double[][] lod_2d = [[4.0, 5.0], [5.0, 2.6], [4.9, 0.3], [4.9, 5.3], [2.8, 2.8]];
  auto peak_2d = get_peak_scanone(lod_2d, map);
  assert(peak_2d[0][0] == 5.0 && peak_2d[0][1].name == "m2");
  assert(peak_2d[1][0] == 5.3 && peak_2d[1][1].name == "m4");

  // multiple positions sharing maximum LOD
  lod = [1.0, 5.6, 2.3, 5.6, 1.8];

  Random gen;
  gen.seed(51118652);
  peak = get_peak_scanone(lod, map, gen);
  assert(peak[0] == 5.6 && peak[1].name == "m4");

  gen.seed(20276563);
  peak = get_peak_scanone(lod, map, gen);
  assert(peak[0] == 5.6 && peak[1].name == "m2");
}

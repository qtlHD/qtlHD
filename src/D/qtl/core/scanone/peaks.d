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

import qtl.core.chromosome;
import qtl.core.primitives;

// find maximum lod score and the position at which it occurred
Tuple!(double, Marker) get_peak_scanone(double[] lod, Marker[] map)
in {
  assert(lod.length == map.length, "lod and map should have the same length");
}
body {
  auto index = new size_t[](lod.length);

  // FIX ME: if multiple locations having the maximum LOD, pick random one?
  makeIndex!("b < a")(lod, index);

  return Tuple!(double, Marker)(lod[index[0]], map[index[0]]);
}

// same for doubly-index LOD scores
Tuple!(double, Marker)[] get_peak_scanone(double[][] lod, Marker[] map)
in {
  assert(lod.length == map.length, "lod and map should have the same length");
}
body {
  Tuple!(double, Marker)[] peaks;

  foreach(i; 0..lod[0].length) {
    auto thislod = new double[](lod.length);
    foreach(j; 0..lod.length)
      thislod[j] = lod[j][i];

    peaks ~= get_peak_scanone(thislod, map);
  }

  return peaks;
}

unittest {
  writeln("Unit test " ~ __FILE__);

  import qtl.core.simulate.map;
  import qtl.core.chromosome;

  double[] lod = [1.0, 3.0, 2.3, 5.6, 1.8];

  auto chromosome = get_chromosome_with_id("1");
  auto map = generate_map_eqspacing_onechr(100.0, cast(uint)lod.length, 1, chromosome);

  auto peak = get_peak_scanone(lod, map);
  assert(peak[0] == 5.6 && peak[1].name == "m4");

  double[][] lod_2d = [[4.0, 5.0], [5.0, 2.6], [4.9, 0.3], [4.9, 5.3], [2.8, 2.8]];
  auto peak_2d = get_peak_scanone(lod_2d, map);
  assert(peak_2d[0][0] == 5.0 && peak_2d[0][1].name == "m2");
  assert(peak_2d[1][0] == 5.3 && peak_2d[1][1].name == "m4");
}
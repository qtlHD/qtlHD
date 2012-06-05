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

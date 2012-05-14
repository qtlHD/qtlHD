/**
 * Functions for the map 
 */

module qtl.core.map.map;

import qtl.core.primitives;
import qtl.core.map.genetic_map_functions;

import std.algorithm;
import std.math;
import std.array;
import std.stdio;

/**
 * Find first (most-left) marker - do not assume markers is sorted
 */

Marker map_first(Ms)(in Ms markers) {
  return reduce!("(a.get_position()<b.get_position()?a:b)")(markers);
}

/**
 * Find last (most-right) marker - do not assume markers is sorted
 */

Marker map_last(Ms)(in Ms markers) {
  return reduce!("(a.get_position()>b.get_position()?a:b)")(markers);
}

/**
 * Calculate map size of a list of Markers - do not assume
 * the markers are ordered.
 */

double map_size(Ms)(in Ms markers) {
  return map_last(markers).get_position - map_first(markers).get_position;
}

/** 
 * Calculate recombination fractions between markers
 */

double[] recombination_fractions(OrderedMs)(OrderedMs markers, in GeneticMapFunc which_genmapfunc = GeneticMapFunc.Haldane) {
  // ps = positions(name)
  // ds = []
  // ps.each_cons(2) { | a | ds.push a[1]-a[0] }
  double[] distances;
  distances.reserve(markers.length);
  Marker leftmarker = null;
  foreach(m; markers) {
    if (leftmarker) {
      distances ~= m.get_position - leftmarker.get_position;
    }
    leftmarker = m;
  }

  return array(dist_to_recfrac(distances, which_genmapfunc)); // Haldane
}

// unittests in test_scanone.d

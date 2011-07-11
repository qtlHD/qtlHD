/**
 * Functions for the map 
 */

module qtl.core.map;

import qtl.core.primitives;

import std.algorithm;
import std.math;
import std.stdio;

/**
 * Find first (most-left) marker - do not assume markers is sorted
 */

Marker map_first(Ms)(Ms markers) {
  return reduce!("(a.get_position()<b.get_position()?a:b)")(markers);
}

/**
 * Find last (most-right) marker - do not assume markers is sorted
 */

Marker map_last(Ms)(Ms markers) {
  return reduce!("(a.get_position()>b.get_position()?a:b)")(markers);
}

/**
 * Calculate map size of a list of Markers - do not assume
 * the markers are ordered.
 */

double map_size(Ms)(Ms markers) {
  return map_last(markers).get_position - map_first(markers).get_position;
}

/** 
 * Calculate recombination between markers
 */

double[] recombination_fractions(Ms)(Ms markers) {
 return null;
}

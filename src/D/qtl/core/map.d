/**
 * Functions for the map 
 */

module qtl.core.marker;

import std.algorithm;
import std.math;
import std.stdio;

/**
 * Calculate map size of a list of Markers - do not assume
 * the markers are ordered.
 */

double map_size(Ms)(Ms markers) {
  auto min = reduce!("(a.get_position()<b.get_position()?a:b)")
    (markers);
  auto max = reduce!("(a.get_position()>b.get_position()?a:b)")
    (markers);
  return max.get_position - min.get_position;
}


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

import std.conv;
import std.container;
import std.typecons;


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

double[] recombination_fractions(OrderedMs)(in OrderedMs markers, in GeneticMapFunc which_genmapfunc = GeneticMapFunc.Haldane) {
  auto distances = new double[markers.length-1];
  foreach(i; 0..(markers.length-1)) {
    distances[i] = markers[i+1].get_position - markers[i].get_position;
  }

  return dist_to_recfrac(distances, which_genmapfunc);
}

// Same, taking tuple
Tuple!(double[])[] recombination_fractions(OrderedMs)(in Tuple!(Chromosome, OrderedMs[])[] markers_by_chr, 
                                                        in GeneticMapFunc which_genmapfunc = GeneticMapFunc.Haldane)
{
  Tuple!(double[])[] list;

  foreach(chr; markers_by_chr) {
    auto recfrac = recombination_fractions(chr[1], which_genmapfunc);
    list ~= Tuple!(double[])(recfrac);
  }

  return list;
}



// unittests in test_scanone.d

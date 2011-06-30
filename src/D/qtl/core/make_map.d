/**
 * Module for making a new map with inserted inter-marker locations
 */

module qtl.core.make_map;

import std.container;
import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.exception;
import std.algorithm;
alias std.algorithm.find find;

import qtl.core.primitives;
import qtl.core.genotype;

/**
 * The make_fixed_map function creates a new map for a chromosome, and 
 * inserts Pseudo markers at fixed positions.
 *
 * In R/qtl one function existed for fixed_distance_map, 
 * fixed_distance_map_sex and add_marker_if_single.
 *
 * markers:        List of (unsorted) markers (one chromosome)
 * step:           Step size (in cM)
 * off_end:        Max distance (in cM) past the terminal marker
 *
 * Returns:
 *
 *   Marker container Ms, including pseudomarkers. Note: this list
 *   is unsorted.
 *
 * (status: under review)
 */

Ms add_stepped_markers_autosome(Ms)(in Ms markers, Position step=1.0, Position off_end=0.0) {
  enforce(step>0);
  enforce(off_end>=0);
  // With R/qtl, if step is zero, the purpose was to add a marker if there is
  // only one. Now, use the function add_marker_if_single instead. Adding
  // markers (at off_end) with step=0 is not supported here. That should also
  // be a separate function.
  auto new_markers = new Ms(markers);
  auto sorted_markers = markers.sorted();
  if (markers.list.length > 1) {
    auto list = sorted_markers.list;
    auto minpos = list[0].get_position();
    auto maxpos = list[$-1].get_position();
    for (auto npos = minpos ; npos < maxpos; npos += step) {
      auto pm = new PseudoMarker(npos);
      if (countUntil(sorted_markers.list, pm) == -1)
        // marker does not exist: add pseudo marker 
        new_markers.add(pm);
    }
    // always add one marker beyond each end (no matter
    // the fixed position distance)
    new_markers.add(new PseudoMarker(minpos - off_end));
    new_markers.add(new PseudoMarker(maxpos + off_end));
    // FIXME remove Pseudo markers too close to other markers
    // (was this in R/qtl?)
  }
  return new_markers;
}

/**
 * NYI FIXME - implementation of sex chromosome will be done later
 */

Ms add_stepped_markers_sex(Ms)(in Ms markers, Position step=1.0, Position off_end=0.0)
{
  throw new Exception("Function not implemented");
}

/**
 * Add a marker if there is only one. Ms is a list of markers, and 
 * this function works as long as there markers.list and a
 * markers.add function are defined in the Ms type(!)
 */

Ms add_one_if_single_marker(Ms)(in Ms markers, Position step_right=1.0) {
  enforce(step_right>0);
  auto new_markers = new Ms(markers);
  if (markers.list.length == 1) {
    // create a single new marker to the right
    auto marker = new_markers.list[0];
    auto npos = marker.get_position() + step_right;
    auto pm = new PseudoMarker(npos);
    new_markers.add(pm);
  }
  return new_markers;
}

/**
 * Add a range of markers if there is only one in the list
 *
 * off_end:   The map is filled with stepped markers from
 *            marker.position-off_end to marker.position+off_end. 
 */

Ms add_stepped_if_single_marker(Ms)(in Ms markers, Position step=1.0, Position off_end=0.0) {
  enforce(step>0);
  enforce(off_end>=0);
  auto new_markers = new Ms(markers);
  if (markers.list.length == 1) {
    auto marker = new_markers.list[0];
    auto position = marker.get_position();
    for (auto npos = position-off_end ; npos < position+off_end ; npos += step) {
      auto pm = new PseudoMarker(npos);
      new_markers.add(pm);
    }
  }
  return new_markers;
}

unittest {
  writeln("Unit test " ~ __FILE__);
  auto markers = new Markers!(MarkerRef!F2)();
  markers.list ~= new MarkerRef!F2(10.0);
  assert(markers.list.length == 1);
  // call function. Note we don't have to give the type!
  auto new_markers = add_one_if_single_marker(markers,2.0);
  // make sure the original list did not change...
  assert(markers.list.length == 1, "Length is " ~ to!string(markers.list.length));
  // now test the new list
  assert(new_markers.list.length == 2, "Length is " ~ to!string(new_markers.list.length));
  assert(new_markers.list[1].marker.position == 12.0);
  assert(new_markers.list[1].marker.name == "loc12", new_markers.list[1].marker.name);
  // It should work for any list of markers
  auto markers1 = new Markers!(Marker)();
  markers1.list ~= new Marker(10.0);
  assert(markers1.list.length == 1);
  // call function. Note we don't have to give the type!
  auto new_markers1 = add_one_if_single_marker(markers1,2.0);
  // make sure the original list did not change...
  assert(markers1.list.length == 1, "Length is " ~ to!string(markers.list.length));
  // now test the new list
  assert(new_markers1.list.length == 2, "Length is " ~ to!string(new_markers.list.length));
  assert(new_markers1.list[1].position == 12.0);
  assert(new_markers1.list[1].name == "loc12", new_markers1.list[1].name);

  auto new_markers2 = add_stepped_if_single_marker(markers,1.0,5.0);
  assert(new_markers2.list.length == 11, "Length is " ~ to!string(new_markers2.list.length));
  // --- test autosome pseudomarker insertion
  auto markers2 = new Markers!(Marker)();
  // start with three markers
  markers2.list ~= new Marker(10.0);
  markers2.list ~= new Marker(20.0);
  markers2.list ~= new Marker(30.0);
  auto res2 = add_stepped_markers_autosome(markers2,1.0,1.0);
  auto list = res2.list; // new marker list
  // assert there are no duplicates
  assert(list.length == 23, to!string(list.length)); 
  auto uniq_list = uniq!"a.get_position() == b.get_position()"(list);
  auto pos_list = map!"a.get_position()"(uniq_list);
  assert(equal(pos_list,[10, 20, 30, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29, 9, 31]), to!string(pos_list));
  auto ts = map!"to!double(a)"(pos_list);
  // nudge mar Result into a double[]
  double ds[];
  foreach (u ; pos_list) { ds ~= u; }  // this can be done better, I am sure
  assert(ds[0] == 10);
  assert(ds.length == 23); // tis proof 
}


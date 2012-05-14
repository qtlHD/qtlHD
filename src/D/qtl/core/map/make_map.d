/**
 * Module for making a new map with inserted inter-marker locations
 */

module qtl.core.map.make_map;

import std.container;
import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.exception;
import std.algorithm;
import std.math;
alias std.algorithm.find find;

import qtl.core.marker;
import qtl.core.primitives;
import qtl.core.chromosome;

/**
 * In R/qtl one function existed for create.map (in util.R), now split into
 * fixed_distance_map, fixed_distance_map_sex and add_marker_if_single.
 */

 /**
 * The make_fixed_map function creates a new map for a chromosome, and
 * inserts Pseudo markers at fixed positions.
 *
 * markerlist:     List of (unsorted) markers (one chromosome); i.e., Marker[]
 * step:           Step size (in cM)
 * off_end:        Max distance (in cM) past the terminal marker
 *
 * Returns:
 *
 *   Marker[] container Ms, including pseudomarkers. Note: this list
 *   is unsorted.
 *
 * (status: under review)
 */

Ms add_stepped_markers_autosome(Ms)(in Ms markerlist, in Position step=1.0, in Position off_end=0.0) {
  enforce(step>0);
  enforce(off_end>=0);
  // With R/qtl, if step is zero, the purpose was to add a marker if there is
  // only one. Now, use the function add_marker_if_single instead. Adding
  // markers (at off_end) with step=0 is not supported here. That should also
  // be a separate function. Variable step size is, again, another function.

  // Note: it wasn't so much the "purpose" for the case of one marker; that was a hack
  //       to deal with bugs that came up when there was just one marker/pseudomarker
  //       on a chromosome
  auto new_markerlist = markerlist.dup();

  sort(new_markerlist);
  auto minpos = new_markerlist[0].get_position();
  auto maxpos = new_markerlist[$-1].get_position();

  if (new_markerlist.length > 1) {
    for (auto npos = minpos+step ; npos < maxpos; npos += step) {
      auto pm = new PseudoMarker(new_markerlist[0].chromosome, npos);
      if (countUntil(new_markerlist, pm) == -1)
        // marker does not exist: add pseudo marker
        new_markerlist ~= pm;
    }
    // FIXME remove Pseudo markers too close to other markers
    // (was this in R/qtl? No, so maybe I leave this)
  }

  if(off_end >= step) {
    for(auto npos=minpos-step; npos >= minpos - off_end; npos -= step)
      new_markerlist ~= new PseudoMarker(new_markerlist[0].chromosome, npos);
    for(auto npos=maxpos+step; npos <= maxpos + off_end; npos += step)
      new_markerlist ~= new PseudoMarker(new_markerlist[0].chromosome, npos);
  }

  // sort the result
  sort(new_markerlist);

  return new_markerlist;
}

// like add_stepped_markers_autosome, but add minimal number of pseudomarkers so that gaps < step
Ms add_minimal_markers_autosome(Ms)(in Ms markerlist, in Position step=1.0, in Position off_end=0.0) {
  enforce(step>0);
  enforce(off_end>=0);

  auto new_markerlist = markerlist.dup();

  sort(new_markerlist);

  auto minpos = new_markerlist[0].get_position();
  auto maxpos = new_markerlist[$-1].get_position();

  if(off_end > 0) {
    new_markerlist ~= new PseudoMarker(new_markerlist[0].chromosome, minpos-off_end);
    new_markerlist ~= new PseudoMarker(new_markerlist[0].chromosome, maxpos+off_end);

    // update minpos and maxpos and sorted list
    minpos -= step;
    maxpos += step;
    sort(new_markerlist);
  }

  if (new_markerlist.length > 1) {

    auto n_markers = new_markerlist.length;

    for(auto left=0; left < n_markers-1; left++) {
      auto leftpos = new_markerlist[left].get_position();
      auto rightpos = new_markerlist[left+1].get_position();
      auto dist =  rightpos - leftpos;

      if(dist > step) {
	auto n_pseudomarkers = ceil(dist/step)-1;
	auto dist_to_step = dist/(n_pseudomarkers+1);

	for(auto pmarpos=leftpos+dist_to_step; pmarpos < rightpos-dist_to_step/2; pmarpos += dist_to_step)
	  new_markerlist ~= new PseudoMarker(new_markerlist[0].chromosome, pmarpos);

      }
    }

    sort(new_markerlist);
  }

  return new_markerlist;
}


/**
 * NYI FIXME - implementation of sex chromosome will be done later
 */

Ms add_stepped_markers_sex(Ms)(in Ms markers, in Position step=1.0, in Position off_end=0.0)
{
  throw new Exception("Function not implemented");
}

/**
 * Add a marker if there is only one. Ms is a list of markers, and
 * this function works as long as there markers.list and a
 * markers.add function are defined in the Ms type(!)
 */

Ms add_one_if_single_marker(Ms)(in Ms markers, in Position step_right=1.0) {
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
 *            marker.get_position()-off_end to marker.get_position()+off_end.
 */

Ms add_stepped_if_single_marker(Ms)(in Ms markers, in Position step=1.0, in Position off_end=0.0) {
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
  auto markers = new Markers!Marker();
  markers.list ~= new Marker(10.0);
  assert(markers.list.length == 1);
  // call function. Note we don't have to give the type!
  auto new_markers = add_one_if_single_marker(markers,2.0);
  // make sure the original list did not change...
  assert(markers.list.length == 1, "Length is " ~ to!string(markers.list.length));
  // now test the new list
  assert(new_markers.list.length == 2, "Length is " ~ to!string(new_markers.list.length));
  assert(new_markers.list[1].get_position() == 12.0);
  assert(new_markers.list[1].name == "loc12", new_markers.list[1].name);
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
  assert(new_markers1.list[1].get_position() == 12.0);
  assert(new_markers1.list[1].name == "loc12", new_markers1.list[1].name);

  auto new_markers2 = add_stepped_if_single_marker(markers,1.0,5.0);
  assert(new_markers2.list.length == 11, "Length is " ~ to!string(new_markers2.list.length));

  writeln("test autosome pseudomarker insertion:");
  // --- test autosome pseudomarker insertion
  auto markers2 = new Markers!(Marker)();
  // start with three markers
  markers2.list ~= new Marker(10.0);
  markers2.list ~= new Marker(20.0);
  markers2.list ~= new Marker(30.0);
  auto res2 = add_stepped_markers_autosome(markers2.list,1.0,1.0);
  auto list = res2; // new marker list
  // assert there are no duplicates
  assert(list.length == 23, to!string(list.length));
  auto uniq_list = uniq!"a.get_position() == b.get_position()"(list);
  auto pos_list = map!"a.get_position()"(uniq_list);
  assert(equal(pos_list,[9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31]), to!string(pos_list));
  auto ds = std.array.array(map!"to!double(a)"(pos_list));
  assert(ds[1] == 10);
  assert(ds.length == 23); // tis proof

  writeln("test with off_end > step:");
  // test with off_end > step
  auto markers3 = new Markers!(Marker)();
  markers3.list ~= new Marker(10.0);
  markers3.list ~= new Marker(20.0);
  markers3.list ~= new Marker(30.0);
  auto res3 = add_stepped_markers_autosome(markers3.list, 5.0, 7.5);
  auto list3 = res3;
  assert(list3.length == 7, to!string(list3.length));
  auto uniq_list3 = uniq!"a.get_position() == b.get_position()"(list3);
  auto pos_list3 = map!"a.get_position()"(uniq_list3);
  assert(equal(pos_list3, [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0]), to!string(pos_list3));

  writeln("test add_minimal_markers_autosome:");
  // test add_minimal_markers_autosome
  auto markers4 = new Markers!(Marker)();
  markers4.list ~= new Marker(10.0);
  markers4.list ~= new Marker(17.0);
  markers4.list ~= new Marker(27.0);
  auto res4 = add_minimal_markers_autosome(markers4.list, 2.0, 7.5);
  auto list4 = res4;
  assert(list4.length == 18, to!string(list4.length));
  auto uniq_list4 = uniq!"a.get_position() == b.get_position()"(list4);
  auto pos_list4 = map!"a.get_position()"(uniq_list4);
  assert(equal(pos_list4, [2.5, 4.375, 6.25, 8.125, 10.0, 11.75, 13.5, 15.25, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0,
			   28.875, 30.75, 32.625, 34.5]),
	 to!string(pos_list4));

  writeln("test add_minimal_markers_autosome with off_end=0:");
  // test add_minimal_markers_autosome with off_end=0
  auto markers5 = new Markers!(Marker)();
  markers5.list ~= new Marker(5.0);
  markers5.list ~= new Marker(8.3);
  markers5.list ~= new Marker(9.8);
  auto res5 = add_minimal_markers_autosome(markers5.list, 1.0, 0);
  auto list5 = res5;
  assert(list5.length == 7, to!string(list5.length));
  auto uniq_list5 = uniq!"a.get_position() == b.get_position()"(list5);
  auto pos_list5 = map!"a.get_position()"(uniq_list5);
  assert(equal(to!string(pos_list5), to!string([5.0, 5.825, 6.65, 7.475, 8.3, 9.05, 9.8])),
	 to!string(pos_list5));

  writeln("end of unit tests in make_map.d");
}


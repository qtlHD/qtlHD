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

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.deprecate.genotype_enum;

/**
 * In R/qtl one function existed for create.map (in util.R), now split into
 * fixed_distance_map, fixed_distance_map_sex and add_marker_if_single.
 */

 /**
 * The make_fixed_map function creates a new map for a chromosome, and 
 * inserts Pseudo markers at fixed positions.
 *
 * markerlist:     List of (unsorted) markers (one chromosome)
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

Ms add_stepped_markers_autosome(Ms)(in Ms markerlist, Position step=1.0, Position off_end=0.0) {
  enforce(step>0);
  enforce(off_end>=0);
  // With R/qtl, if step is zero, the purpose was to add a marker if there is
  // only one. Now, use the function add_marker_if_single instead. Adding
  // markers (at off_end) with step=0 is not supported here. That should also
  // be a separate function. Variable step size is, again, another function.

  auto new_markerlist = cast(Ms)markerlist.xdup;
  writeln("**",typeid(new_markerlist));
  // **const(qtl.core.primitives.Marker)[]
  /*
  Ms new_markerlist;
  foreach(m; markerlist) {
    new_markerlist ~= cast(Marker)m;
  }
  */

/*
  sort_markers_by_position(new_markerlist);
  auto minpos = new_markerlist[0].get_position();
  auto maxpos = new_markerlist[$-1].get_position();

  if (new_markerlist.length > 1) {
    for (auto npos = minpos+step ; npos < maxpos; npos += step) {
      auto pm = new PseudoMarker(npos);
      if (countUntil(new_markerlist, pm) == -1)
        // marker does not exist: add pseudo marker 
        new_markerlist ~= pm;
    }
    // FIXME remove Pseudo markers too close to other markers
    // (was this in R/qtl? No, so maybe I leave this)
  }

  if(off_end >= step) {  
    for(auto npos=minpos-step; npos >= minpos - off_end; npos -= step) 
      new_markerlist ~= new PseudoMarker(npos);
    for(auto npos=maxpos+step; npos <= maxpos + off_end; npos += step) 
      new_markerlist ~= new PseudoMarker(npos);
  }

  // sort the result
  sort_markers_by_position(new_markerlist);

  */
  return new_markerlist;
}

// like add_stepped_markers_autosome, but add minimal number of pseudomarkers so that gaps < step
Ms add_minimal_markers_autosome(Ms)(in Ms markers, Position step=1.0, Position off_end=0.0) {
  enforce(step>0);
  enforce(off_end>=0);

  auto new_markers = new Ms(markers);
  auto sorted_markers = new_markers.sorted();

  auto list = sorted_markers.list;
  auto minpos = list[0].get_position();
  auto maxpos = list[$-1].get_position();

  if(off_end > 0) {  
    new_markers.add(new PseudoMarker(minpos-off_end));
    new_markers.add(new PseudoMarker(maxpos+off_end));

    // update minpos and maxpos and sorted list
    minpos -= step;
    maxpos += step;
    sorted_markers = new_markers.sorted();
    list = sorted_markers.list;
  }

  if (list.length > 1) {
    for(auto left=0; left < list.length-1; left++) {
      auto leftpos = list[left].get_position();
      auto rightpos = list[left+1].get_position();
      auto dist =  rightpos - leftpos;

      if(dist > step) {
	auto n_pseudomarkers = ceil(dist/step)-1;
	auto dist_to_step = dist/(n_pseudomarkers+1);

	for(auto pmarpos=leftpos+dist_to_step; pmarpos < rightpos; pmarpos += dist_to_step) {
	  new_markers.add(new PseudoMarker(pmarpos));
	}
      }
    }

    new_markers = new_markers.sorted();
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

  // test add_minimal_markers_autosome
  auto markers4 = new Markers!(Marker)();
  markers4.list ~= new Marker(10.0);
  markers4.list ~= new Marker(17.0);
  markers4.list ~= new Marker(27.0);
  auto res4 = add_minimal_markers_autosome(markers4, 2.0, 7.5);
  auto list4 = res4.list;
  assert(list4.length == 18, to!string(list4.length));
  auto uniq_list4 = uniq!"a.get_position() == b.get_position()"(list4);
  auto pos_list4 = map!"a.get_position()"(uniq_list4);
  assert(equal(pos_list4, [2.5, 4.375, 6.25, 8.125, 10.0, 11.75, 13.5, 15.25, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0,
			   28.875, 30.75, 32.625, 34.5]), 
	 to!string(pos_list4));

  // test add_minimal_markers_autosome with off_end=0
  auto markers5 = new Markers!(Marker)();
  markers5.list ~= new Marker(5.0);
  markers5.list ~= new Marker(8.3);
  markers5.list ~= new Marker(9.8);
  auto res5 = add_minimal_markers_autosome(markers5, 1.0, 0);
  auto list5 = res5.list;
  assert(list5.length == 7, to!string(list5.length));
  auto uniq_list5 = uniq!"a.get_position() == b.get_position()"(list5);
  auto pos_list5 = map!"a.get_position()"(uniq_list5);
  assert(equal(to!string(pos_list5), to!string([5.0, 5.825, 6.65, 7.475, 8.3, 9.05, 9.8])), 
	 to!string(pos_list5));
}


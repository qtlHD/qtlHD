/**
 * Marker related functions
 */


module qtl.core.marker;

import std.container;
import std.conv;
import std.variant;
import std.random;
import qtl.core.primitives;
import std.algorithm; // to use sort()

import std.stdio;
// import std.container;
import std.typecons;

/**
 * The Markers container keeps track of Markers, throwing in some useful
 * funtionality. Note there is no guarantee the list is ordered.
 *
 * Normally use the default Marker[] instead.
 */

class Markers(M) {
  mixin ActList!M; 
  this() {}
  this(in Markers!M markers) {
    // list = markers.list.dup;  // make sure to clone all data
    foreach(m ; markers.list) {
      list ~= cast(M) m;
    }
  }
  this(M[] markers) {
    list = markers.dup;
  }
  void add(in Marker m) {
    list ~= new M(m);
  }

  const auto sorted() { 
    auto ms = new Markers(this); // make a copy
    sort(ms.list); // sorts in place
    return ms;
  }
  const size_t length() { return list.length; }
  M opIndex(int i) { return list[i]; }
  /// find by name
  M find(string name) {
    foreach(m; list) { if (m.name == name) return m; }
    return null;
  }
  /// find by ID
  M find(uint id) {
    foreach(m; list) { if (m.id == id) return m; }
    return null;
  }
}

@property M[] list(M)(Markers!M ms) { return ms.list; }

/**
 * Convenience function for duplication Marker[]
 */

Marker[] dup(in Marker[] ms) 
{ 
  Marker[] new_ms;
  foreach(m; ms) {
    new_ms ~= cast(Marker)m;
  }
  return new_ms;
}

unittest {
  writeln("Unit test " ~ __FILE__);
  Marker m1 = new Marker(4.6,"m1",1);
  Marker m2 = new Marker(4.8,"m2",2);
  m2.chromosome = new Autosome("1",1);
  m2.attrib_list = new Attribute[1];
  auto mref1 = new MarkerRef!uint(m1);
  auto mref2 = new MarkerRef!uint(m2);
  PseudoMarker pm1 = new PseudoMarker(4.7,"pm1",3);
  auto pmref1 = new MarkerRef!uint(pm1);
  auto markers = new Markers!(MarkerRef!uint)();
  markers.list ~= mref1;  // 1, 4.6
  markers.list ~= mref2;  // 2, 4.8
  markers.list ~= pmref1; // 3, 4.7
  // test list
  assert(list(markers) == markers.list);
  // find by index
  assert(markers[0].name == "m1");
  // find by name
  assert(markers.find("m1").name == "m1");
  uint[] result;
  foreach ( m ; markers ) {
    result ~= m.marker.id;
  }
  assert(result==cast(uint[])[1,2,3]);
  // reverse sort markers
  auto list2 = sort!("a.get_position() > b.get_position()")(markers.list);
  foreach ( m ; list2 ) {
    result ~= m.marker.id; // add reverse sorted list (2,3,1)

  }
  assert(result==cast(uint[])[1,2,3,2,3,1],to!string(result));
  auto ms = markers.sorted(); 
  foreach ( m ; ms ) {
    result ~= m.marker.id; // add sorted marker list (1,3,2)
  }
  assert(result==cast(uint[])[1,2,3,2,3,1,1,3,2],to!string(result));
  auto ulist2 = uniq!("a.get_position() == b.get_position")(ms.list);
  auto pos_list = map!"a.get_position()"(ulist2);
  assert(equal(pos_list,[4.6, 4.7, 4.8]),to!string(pos_list));
}
